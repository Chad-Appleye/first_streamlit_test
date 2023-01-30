import sys
import time
import subprocess as sp
import argparse
import multiprocessing as mp
from collections import defaultdict
import pysam
import json
import os
import io

parser = argparse.ArgumentParser(description="Drive segmented diagnostic databases")
parser.add_argument("input_binner", type=str, help="input binner file from txtll output, with reads named <seq id>.<pos>.<read id>")
parser.add_argument("config_file", type=str, help="config file name")
parser.add_argument("namespace", type=str, help="specify namespace.  can be one of [16S, 18S, viral, bg]")
parser.add_argument("parent_taxid", type=int, help="parent taxid used to subset references, likely parent of bacteria (2), fungi (4751), or virus (10239)")
parser.add_argument("n_processes", type=int, help="number of processes to run concurrently")
parser.add_argument("output_directory", type=str, help="output directory for simulations, will create seqid specific outputs under this")
parser.add_argument("--v4_binner", action='store_true', help="set if the binner output is from v4", default=False)
args = parser.parse_args()

"""
json config file

top_level_dir:          path to top level folder that contains all the relevant databases, resources path are relative to this point
ncbi_tax_path:          path to directory that contains nodes.dmp, names.dmp, merged.dmp file from ncbi
ncbi_tax_utils:         path to directory that contains ncbi_taxonomy_utils.py
binner_dna_executable:  complete path to binner executable
binner_dna_db:          complete path to binner database.  If 'ramdisk' is in the path, then it will not load db into ram
binner_rna_executable:  complete path to binner executable
binner_rna_db:          complete path to binner database.  If 'ramdisk' is in the path, then it will not load db into ram
ncbi_dna_map_file:      path to file containing mapping of taxids from binner output to ncbi taxid and category.  <taxid> <ncbi taxid> <category> -- necessary because remappings occur in the top level
taxid_reporting_names:  path to json file containing taxid, reporting name, organism class, subclass, medical relevance, etc
rna_leaf_ref_files:     complete path to file with <reference fasta><tab><taxid>
bact_16S_tgs:           complete path to bacterial 16S .tgs file
fp_18S_tgs:             complete path to fungal / parasite .tgs file
viral_tgs:              complete path to viral .tgs file
viral_leaf_files:       complete path to file with <reference fasta><tab><taxid>
bact_gen_tgs:           complete path to bacterial genomes .tgs file
bact_gen_leaf_files     complete path to file with <reference fasta><tab><taxid>
viral_prot_db:          complete path to viral protein database
txdx_sm_binary:         complete path to txdx_sm
species_read_cutoff     number of reads required to use a species in alignment
bact_16S_txdx_params:   bacterial 16S alignment parameters for txdx_sm as dict:  {c:<int>,x:<int>,z:<int>,m:<float>,k:<int>}
bact_16S_txdx_db:       pre-built 16S txdx database
bact_16S_txdx_fa:       file with all 16S sequences in pre-built 16S database, not to be built on the fly
bact_gen_txdx_params:   bacterial genome alignment parameters for txdx_sm as dict:  {c:<int>,x:<int>,z:<int>,m:<float>,k:<int>}
fp_18S_txdx_params:     fungal / parasite alignment parameters for txdx_sm as dict:  {c:<int>,x:<int>,z:<int>,m:<float>,k:<int>}
18S_txdx_fa:            fungal / parasite 18S reference sequence fasta will ALL sequences, for txdx_sm
18S_txdx_db:            fungal / parasite 18S database pre-built with 18S_txdx_fa using txdx_sm
viral_txdx_params:      viral alignment parameters for txdx_sm as dict:  {c:<int>,x:<int>,z:<int>,m:<float>,k:<int>}
viral_txdx_db:          pre-built viral txdx database
viral_txdx_fa:          file with all viral sequences, not to be constructed on the fly
fungal_gen_ba_db:       complete path to fungal genome bit array database
fungal_gen_meta:        complete path to fungal genome meta info file 
parasite_gen_ba_db:     complete path to parasite genome bit array database
parasite_gen_meta:      complete path to parasite genome meta info file
bit_search_executable:  complete path to bit search executable
binner_ncpus:           number of cpus to use while binning
align_ncpus:            number of cpus for txdx_sm and ba_search to use, set to around 10
process_binner_ncpus:   number of independent processes to use while parsing binner output
dna_ncpus:              number of independent processes to use while running txdx_sm and ba_search (for dna data)
rna_ncpus:              number of independent processes to use while running txdx_sm (for rna data)
summary_info:           information for summary file including the version
"""

config_args = json.load(open(args.config_file, 'r'))
if "binner_ncpus" not in config_args:
    config_args["binner_ncpus"] = os.cpu_count()  # if binner ncpus is not set, use all cpus on machine

# change working directory to top level that contains all the database paths
os.chdir(config_args["top_level_dir"])

taxid_reporting_name_info = json.load(open(config_args["taxid_reporting_names"], 'r'))

"""
load ncbi taxonomy
"""
sys.path.append(config_args["ncbi_tax_utils"])
from ncbi_taxonomy_utils import *
ncbi_taxonomy_path = config_args["ncbi_tax_path"].rstrip("/")
ncbi_tax = ncbi_taxonomy("%s/merged.dmp" % ncbi_taxonomy_path, "%s/nodes.dmp" % ncbi_taxonomy_path, "%s/names.dmp" % ncbi_taxonomy_path)

"""
load fasta sequences by leaf taxid
"""

leaf_taxids = {}  # {taxid}->StringIO of reference sequences
nuc_file_taxids = config_args["rna_leaf_ref_files"]
if args.namespace == "viral":  # need to use different leaf taxid file since the viral database was updated and not the binner
    nuc_file_taxids = config_args["viral_leaf_files"]
elif args.namespace == "bg":  # need to use different leaf taxid file with the bacterial genomes
    nuc_file_taxids = config_args["bact_gen_leaf_files"]

for line in open(nuc_file_taxids, 'r'):
    data = line.strip().split("\t")
    taxid = ncbi_tax.merged.get(int(data[1]), int(data[1]))
    if taxid not in ncbi_tax.nodes_rel:
        sys.stderr.write("WARNING: taxid %s from reference files not found in taxonomy and will be skipped.\n" % data[1])
        sys.stderr.flush()
        continue
    if args.parent_taxid in ncbi_tax.get_path(taxid):
        for seq in pysam.FastxFile(data[0]):
            leaf_taxids.setdefault(ncbi_tax.merged.get(int(data[1]), int(data[1])), io.StringIO()).write(">%s\n%s\n" % (seq.name, seq.sequence))


def gen_binner_output_file(output_dir, accession, seq_type):
    return "%s/%s.%s.binner.out" % (output_dir, accession, seq_type)


def gen_sample_composition_file(output_dir, accession, seq_type):
    return "%s/%s.%s.sample_composition.out" % (output_dir, accession, seq_type)


def gen_binner_class_fasta_file(output_dir, accession, seq_type, org_class):
    return "%s/%s.%s.%s" % (output_dir, accession, seq_type, org_class)


def gen_done_file(tmp_input):
    return "%s/%s.%s.done" % (tmp_input[2], tmp_input[3], tmp_input[1])


cached_genus_taxids = {}  #{taxid}->genus taxid
def get_genus_taxid(taxid):
    genus_taxid = cached_genus_taxids.get(taxid, None)
    if genus_taxid is None:
        lnr = ncbi_tax.get_lineage_names_ranks(taxid, canonical=True)
        genus_taxid = lnr["genus"][1]
        cached_genus_taxids[taxid] = genus_taxid
    return genus_taxid


cached_taxid_paths = {}  # {taxid}->[path]
def get_taxid_path(taxid):
    path = cached_taxid_paths.get(taxid, None)
    if path is None:
        path = ncbi_tax.get_path(taxid)
        cached_taxid_paths[taxid] = path
    return path


def run_rna_alignment(taxid_propagated_counts, align_reads_fa, dxsm_c=2, dxsm_x=5, dxsm_z=1000, dxsm_m=0.02, dxsm_kl=25):
    """
    :param taxid_propagated_counts:     Dictionary that saves the number of read counts by taxid
    :param align_reads_fa:              Fasta file of reads to be aligned (mapped to this organism class)
    :param dxsm_c                       txdx_sm cutoff for number of kmers from read that have to match reference
    :param dxsm_x:                      txdx_sm cutoff for number of reads going to reference in order to be reported
    :param dxsm_z:                      txdx_sm cutoff for number of references a read can be shared with
    :param dxsm_m:                      txdx_sm cutoff for percent mismatch tolerated in alignment
    :param dxsm_kl:                     txdx_sm kmer length
    :return:
    """

    # collect leaf taxids based on read counts in classifier.  not done in viral or 16S namespace, since all references are used
    cutoff = config_args["species_read_cutoff"]
    if args.namespace != "viral" and args.namespace != "16S" and args.namespace != "18S":
        taxids_to_use = set()
        for ltx in leaf_taxids:
            lineage_name_rank = ncbi_tax.get_lineage_names_ranks(ltx, canonical=True)
            species_taxid = lineage_name_rank.get("species")[1]
            genus_taxid = lineage_name_rank.get("genus")[1]
            total_reads = taxid_propagated_counts.get(ltx, 0)  # use leaf taxid count if species taxid doesnt exist
            if species_taxid > 0:
                tmp_reads = taxid_propagated_counts.get(species_taxid, 0)
                if tmp_reads > total_reads:  # only use if count is greater than previously
                    total_reads = tmp_reads
            if total_reads <= cutoff and (args.namespace == "16S" or args.namespace == "18S"):  # for 16S / 18S check genus level if cutoff not met
                total_reads = taxid_propagated_counts.get(genus_taxid, 0)
            if total_reads > cutoff:  # include leaf taxid if read count threshold is met
                taxids_to_use.add(ltx)

        sys.stderr.write("TOTAL NUMBER OF TAXIDS ABOVE THRESHOLD: %d\n" % len(taxids_to_use))
        sys.stderr.flush()

    # collect relevant references
    dxsm_ref_fasta = align_reads_fa + ".dxsm.ref.fa"
    if args.namespace != "viral" and args.namespace != "16S" and args.namespace != "18S":  # not use for viruses or 16S since all references will be used with viruses
        ref_fa_out = open(dxsm_ref_fasta, 'w')
        for ltx in taxids_to_use:
            seqs = leaf_taxids.get(ltx, None)
            ref_fa_out.write("%s" % seqs.getvalue())
        ref_fa_out.close()

    dxsm_out = align_reads_fa + ".dxsm.out"
    #run txdx_sm
    if args.namespace == "viral":  #viral_txdx_db viral_txdx_fa
        command = "%s -i '%s' -r '%s' -o '%s' -l %s -k %d -c %d -x %d -z %d -m %f -n %s" % (
            config_args["txdx_sm_binary"],
            align_reads_fa,
            config_args["viral_txdx_fa"],
            dxsm_out,
            config_args["viral_txdx_db"],
            dxsm_kl, dxsm_c, dxsm_x, dxsm_z, dxsm_m, config_args["align_ncpus"])
    elif args.namespace == "16S":  #bact_16S_txdx_db bact_16S_txdx_fa
        command = "%s -i '%s' -r '%s' -o '%s' -l %s -k %d -c %d -x %d -z %d -m %f -n %s" % (
            config_args["txdx_sm_binary"],
            align_reads_fa,
            config_args["bact_16S_txdx_fa"],
            dxsm_out,
            config_args["bact_16S_txdx_db"],
            dxsm_kl, dxsm_c, dxsm_x, dxsm_z, dxsm_m, config_args["align_ncpus"])
    elif args.namespace == "18S":
        command = "%s -i '%s' -r '%s' -o '%s' -l %s -k %d -c %d -x %d -z %d -m %f -n %s" % (
            config_args["txdx_sm_binary"],
            align_reads_fa,
            config_args["18S_txdx_fa"],
            dxsm_out,
            config_args["18S_txdx_db"],
            dxsm_kl, dxsm_c, dxsm_x, dxsm_z, dxsm_m, config_args["align_ncpus"])
    else:
        command = "%s -i '%s' -r '%s' -o '%s' -k %d -c %d -x %d -z %d -m %f -n %s" % (config_args["txdx_sm_binary"], align_reads_fa, dxsm_ref_fasta, dxsm_out, dxsm_kl, dxsm_c, dxsm_x, dxsm_z, dxsm_m, config_args["align_ncpus"])
    sys.stderr.write("%s\n" % command)
    sys.stderr.flush()
    sp.call(command, shell=True)


def get_alignment_parameters():
    c = None
    x = None
    z = None
    m = None
    kl = None
    if args.namespace == "viral":
        c = config_args["viral_txdx_params"]["c"]
        x = config_args["viral_txdx_params"]["x"]
        z = config_args["viral_txdx_params"]["z"]
        m = config_args["viral_txdx_params"]["m"]
        kl = config_args["viral_txdx_params"]["k"]
    elif args.namespace == "16S":
        c = config_args["bact_16S_txdx_params"]["c"]
        x = config_args["bact_16S_txdx_params"]["x"]
        z = config_args["bact_16S_txdx_params"]["z"]
        m = config_args["bact_16S_txdx_params"]["m"]
        kl = config_args["bact_16S_txdx_params"]["k"]
    elif args.namespace == "18S":
        c = config_args["fp_18S_txdx_params"]["c"]
        x = config_args["fp_18S_txdx_params"]["x"]
        z = config_args["fp_18S_txdx_params"]["z"]
        m = config_args["fp_18S_txdx_params"]["m"]
        kl = config_args["fp_18S_txdx_params"]["k"]
    elif args.namespace == "bg":
        c = config_args["bact_gen_txdx_params"]["c"]
        x = config_args["bact_gen_txdx_params"]["x"]
        z = config_args["bact_gen_txdx_params"]["z"]
        m = config_args["bact_gen_txdx_params"]["m"]
        kl = config_args["bact_gen_txdx_params"]["k"]

    return c, x, z, m, kl


def process_binner_block(block_info):
    custom_tax_map = None
    if args.namespace == 'bg':
        # load mapping from custom taxonomy to ncbi taxonomy
        custom_tax_map = {}  # {custom taxid}->ncbi_taxid
        for line in open(config_args["ncbi_dna_map_file"], 'r'):
            data = line.strip().split("\t")
            custom_tax_map[data[0]] = int(data[1])
        custom_tax_map['0'] = 0
    offset = block_info[0]
    n_lines = block_info[1]
    seqid = block_info[2]
    seqid_output_dir = args.output_directory.rstrip("/") + "/%s" % seqid
    os.mkdir(seqid_output_dir)
    read_file_name = seqid_output_dir + "/sim_reads.fa"
    out_fa = open(read_file_name, 'w')
    f = open(args.input_binner, 'r')
    f.seek(offset)
    taxid_propagated_count = defaultdict(int)
    for i in range(n_lines):
        line = f.readline()
        data = line.strip().split("\t")
        taxid = int(data[1])
        if custom_tax_map is not None:
            taxid = custom_tax_map[data[1]]  # mapping to the ncbi taxonomy
        if taxid == 1:
            if args.v4_binner is False:
                out_fa.write(">%s\t%d\n%s\n" % (data[0], taxid, data[8]))
            else:
                out_fa.write(">%s\t%d\n%s\n" % (data[0], taxid, data[6]))
        else:
            for tx in get_taxid_path(taxid):
                taxid_propagated_count[tx] += 1
                if args.parent_taxid == tx:
                    if args.v4_binner is False:
                        out_fa.write(">%s\t%d\n%s\n" % (data[0], taxid, data[8]))
                    else:
                        out_fa.write(">%s\t%d\n%s\n" % (data[0], taxid, data[6]))
    f.close()
    out_fa.close()
    c, x, z, m, kl = get_alignment_parameters()
    run_rna_alignment(taxid_propagated_count, read_file_name, dxsm_c=c, dxsm_x=x, dxsm_z=z, dxsm_m=m, dxsm_kl=kl)


"""
Parse binner output
sequence_name   max_taxid       max_window      max_score       n_ties  reading_frame   reading_frame_length    taxid_string    query_sequence
"""
binner_output_indices = []  # [(offset, n lines, seqid)]
cur_seqid = None
cur_offset = 0
block_byte_size = 0
n_lines = 0
f = open(args.input_binner, 'r')
cur_offset += len(f.readline())  # skip header
for line in f:
    data = line.strip().split("\t")
    tmp_seqid = data[0].split(".")[0]
    if tmp_seqid == cur_seqid or cur_seqid is None:
        cur_seqid = tmp_seqid
        block_byte_size += len(line)
        n_lines += 1
    else:
        binner_output_indices.append((cur_offset, n_lines, cur_seqid))
        cur_offset += block_byte_size
        n_lines = 1  # set to 1 to include the current line
        cur_seqid = tmp_seqid
        block_byte_size = len(line)
if n_lines > 0:  # add last block of data
    binner_output_indices.append((cur_offset, n_lines, cur_seqid))
f.close()

p = mp.Pool(args.n_processes)
results = p.imap(process_binner_block, binner_output_indices, 1)
for res in results:
    pass

