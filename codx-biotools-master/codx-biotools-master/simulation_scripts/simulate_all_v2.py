import sys
import time
import subprocess as sp
import argparse
import multiprocessing as mp
from collections import defaultdict
import pysam
import json
import os
import re
import itertools
from sim_reads import *


parser = argparse.ArgumentParser(description="Drive simulations of organisms in input file")
parser.add_argument("input_file", type=str, help="input file where each line is <namespace 16S, 18S, viral, bg, fg, pg><tab><reporting id><tab><reporting name><tab><taxid><tab><rna / dna>")
parser.add_argument("config_file", type=str, help="json file that contains information needed for simulations")
parser.add_argument("output_dir", type=str, help="output directory for all the simulations")
parser.add_argument("--pipeline_tag", type=str, help="tag to include in .summary.done files", default=".")
args = parser.parse_args()

"""
simulation config:
16S_fasta               16S fasta to be used in simulation
16S_tgs                 16S tgs that maps to 16S fasta
16S_sim_depth           depth of coverage to simulate 16S sequences
18S_fasta               18S fasta to be used in simulation
18S_tgs                 18S tgs that maps to 18S fasta
18S_sim_depth           depth of coverage to simulate 18S sequences
viral_fasta             viral fasta to be used in simulation
viral_tgs               viral tgs that maps to viral fasta
viral_sim_depth         depth of coverage to simulate viral sequences
bg_genome_paths         bacterial genome paths file  <genome path><tab><id>
bg_minhash_meta         bacterial genome minhash meta info file
bg_sim_depth            depth of coverage to simulate bacterial genomes
fg_genome_paths         fungal genome paths file <genome path><tab><id>
fg_minhash_meta         fungal genome minhash meta info file
fg_sim_depth            depth of coverage to simulate fungal genomes
pg_genome_paths         parasite genome paths file <genome path><tab><id>
pg_minhash_meta         parasite genome minhash meta info file
pg_sim_depth            depth of coverage to simulate parasite genomes
python_interpreter      python interpreter to use when calling pipeline
pipeline_script         complete path to exbox_pipeline_v2.py
pipeline_config_file    complete path to pipeline config
taxid_reporting_names   complete path to json file with reporting ids, names
ncbi_tax_utils          complete path to folder containing ncbi_taxonomy_utils.py
ncbi_tax_path           complete path to folder containing merged.dmp, nodes.dmp, names.dmp
read_length             read length to simulate
"""

"""
load config files
"""
config = json.load(open(args.config_file, 'r'))
config['ncpus_sim'] = 10

# load taxid reporting name info
taxid_reporting_name_info = {}  # {"taxid"}->{reporting_id:, reporting_name:, compound_id:, class_type:, subclass:}
f = open(config["taxid_reporting_names"])
f.readline()  # skip header
for line in f:
    data = line.strip().split("\t")
    for taxid in data[2].split(","):
        taxid_reporting_name_info[taxid] = {
            "reporting_id": data[1],
            "reporting_name": data[0],
            "compound_id": data[3],
            "class_type": data[4],
            "subclass": [x for x in data[5].split(",")],
            "nucleic_acid": data[6]
        }
f.close()

"""
load ncbi taxonomy
"""
sys.path.append(config["ncbi_tax_utils"])
from ncbi_taxonomy_utils import *
ncbi_taxonomy_path = config["ncbi_tax_path"].rstrip("/")
ncbi_tax = ncbi_taxonomy("%s/merged.dmp" % ncbi_taxonomy_path, "%s/nodes.dmp" % ncbi_taxonomy_path, "%s/names.dmp" % ncbi_taxonomy_path)


def simulate_genomes(info):
    sim_taxid = info['taxid']
    genome_id_info = info['genome_id_info']
    pipeline_inputs = []
    for genome_id, g_info in genome_id_info.items():
        sim_taxid_output_dir = "%s/%s/%d" % (info['output_dir'], info['namespace'], sim_taxid)
        try:
            os.mkdir(sim_taxid_output_dir)
        except FileExistsError:
            pass
        # create simulated reads
        sim_file = "%s/%d.sim.fa" % (sim_taxid_output_dir, genome_id)
        out = open(sim_file, 'w')
        for seq in pysam.FastxFile(g_info[1]):
            simulate_reads_tiled_by_depth("%s-%s" % (g_info[0], seq.name), seq.sequence, info['depth'], info['read_length'], file_handle=out)
        out.close()
        pipeline_inputs.append("%s\t%s\t%s\t%d.sim\t%s" % (sim_file, info['nucl_type'], sim_taxid_output_dir, genome_id, args.pipeline_tag))

    return pipeline_inputs


def combine_simulate_segments(info):
    sim_taxid = info['taxid']
    geneid_info = info['geneid_info']
    geneid_list = []
    seqids_to_use = set()
    for geneid, seqids_accessions in geneid_info.items():
        geneid_list.append(seqids_accessions)
        seqids_to_use.update([x[0] for x in seqids_accessions])

    seqid_seq = {}
    for seq in pysam.FastxFile(info["fasta"]):
        tmp_seqid = int(seq.name)
        if tmp_seqid in seqids_to_use:
            seqid_seq[tmp_seqid] = seq.sequence

    sim_taxid_output_dir = "%s/%s/%d" % (info["output_dir"], info["namespace"], sim_taxid)
    try:
        os.mkdir(sim_taxid_output_dir)
    except FileExistsError:
        pass

    pipeline_inputs = []
    new_seqids = set('need a place holder to start iteration')  # reset at every iteration, describes which new sequences have been used
    used_seqids = set()  # saves which sequences have already been included in a simulation
    sim_index = 0
    while len(new_seqids) > 0:
        new_seqids = set()
        sim_seqids = []
        for geneid, seqids_accessions in geneid_info.items():
            found_seqid = False
            for seq_accession in seqids_accessions:
                if seq_accession[0] not in used_seqids:
                    sim_seqids.append(seq_accession)
                    new_seqids.add(seq_accession[0])
                    found_seqid = True
                    break
            if found_seqid is False:  # just add the first if all have already been used
                sim_seqids.append(seqids_accessions[0])
        if len(new_seqids) > 0:  # only simulate if at least 1 new sequence is present
            sim_file = "%s/%d.sim.fa" % (sim_taxid_output_dir, sim_index)
            out = open(sim_file, 'w')
            for seq in sim_seqids:
                out.write("%s" % ("".join(simulate_reads_tiled_by_depth("%s-%d" % (seq[1], seq[0]), seqid_seq[seq[0]], info['depth'], info['read_length']))))
            out.close()
            pipeline_inputs.append("%s\t%s\t%s\t%d.sim\t%s" % (sim_file, info['nucl_type'], sim_taxid_output_dir, sim_index, args.pipeline_tag))
        sim_index += 1
        used_seqids.update(new_seqids)

    return pipeline_inputs


def drive_sims(namespace, org_info, output_dir, config, overall_out_fh):
    if namespace in ['16S', '18S', 'viral']:
        """
        1) simulate reads
        2) run pipeline
        """
        # simulate reads
        sim_taxid_seqids = {}  # {taxid}->{geneid}->[seqids]
        taxid_org_info = {}  # {taxid}->{'reporting_id':, 'reporting_name':, 'taxid':, 'nucl_type'}
        for line in open(config["%s_tgs" % namespace], 'r'):
            data = line.strip().split("\t")
            taxid = int(data[0])
            geneid = int(data[1])
            seqid = int(data[2])
            json_obj = json.loads(data[3])
            for tx in ncbi_tax.get_path(taxid):
                reporting_info = taxid_reporting_name_info.get("%d" % tx, None)
                if reporting_info is not None and reporting_info["reporting_id"] in org_info:
                    sim_taxid_seqids.setdefault(taxid, {}).setdefault(geneid, []).append((seqid, json_obj["accession"]))
                    oi = org_info[reporting_info["reporting_id"]]
                    taxid_org_info[taxid] = {'reporting_id': reporting_info["reporting_id"], 'reporting_name': reporting_info["reporting_name"], 'taxid': taxid, 'nucl_type': oi["nucl_type"]}
                    break
        # flatten sequence information into lists for multiprocessing simulations
        simulation_info = []
        for taxid, geneid_info in sim_taxid_seqids.items():
            info_dict = {}
            info_dict['geneid_info'] = geneid_info
            info_dict['taxid'] = taxid
            info_dict['output_dir'] = output_dir
            info_dict['fasta'] = config["%s_fasta" % namespace]
            info_dict['read_length'] = config["read_length"]
            info_dict['namespace'] = namespace
            info_dict['depth'] = config["%s_sim_depth" % namespace]
            info_dict.update(taxid_org_info[taxid])
            simulation_info.append(info_dict)

        pipeline_input_lines = []
        p = mp.Pool(config['ncpus_sim'])
        results = p.imap(combine_simulate_segments, simulation_info)
        for res in results:
            pipeline_input_lines.extend(res)
        p.close()
        p.join()

        # run pipeline, process 200 at a time
        for tmp_input in itertools.zip_longest(*[iter(pipeline_input_lines)]*200):
            pipeline_input_file = "%s/%s.input.txt" % (output_dir, namespace)
            out = open(pipeline_input_file, 'w')
            out.write("%s\n" % ("\n".join([x for x in tmp_input if x is not None])))
            overall_out_fh.write("%s\n" % ("\n".join([x for x in tmp_input if x is not None])))
            out.close()

            command = "%s %s %s %s --keep_dxsm --compress_output 1> %s_pipeline.out 2> %s_pipeline.err" % (config['python_interpreter'], config['pipeline_script'], pipeline_input_file, config['pipeline_config_file'], namespace, namespace)
            sys.stderr.write("%s\n" % command)
            sys.stderr.flush()
            sp.call(command, shell=True)

    elif namespace in ['bg', 'fg', 'pg']:

        genome_id_fasta_path = {}  # {genome id}->path to genome fasta
        for line in open(config["%s_genome_paths" % namespace], 'r'):
            data = line.strip().split("\t")
            genome_id_fasta_path[int(data[1])] = data[0]

        sim_org_info = {}  # {taxid}->{'reporting_id':, 'reporting_name':, 'taxid':, 'nucl_type':, 'genome_id_info': {genome id}->[accession, genome path]}
        for line in open(config["%s_minhash_meta" % namespace]):  # <genome index> <accession> <taxid> <name> <length>
            data = line.strip().split("\t")
            taxid = int(data[2])
            for tx in ncbi_tax.get_path(taxid):
                reporting_info = taxid_reporting_name_info.get("%d" % tx, None)
                if reporting_info is not None and reporting_info["reporting_id"] in org_info:
                    oi = org_info[reporting_info["reporting_id"]]
                    sim_info = sim_org_info.setdefault(taxid, {'reporting_id': reporting_info["reporting_id"], 'reporting_name': reporting_info["reporting_name"], 'taxid': taxid, 'nucl_type': oi["nucl_type"], 'genome_id_info': {}})
                    genome_id = int(data[0])
                    sim_info['genome_id_info'][genome_id] = [data[1], genome_id_fasta_path[genome_id]]
        # flatten sequence information into lists for multiprocessing simulations
        simulation_info = []
        for taxid, sim_info in sim_org_info.items():
            info_dict = {}
            info_dict['genome_id_info'] = sim_info["genome_id_info"]
            info_dict['taxid'] = taxid
            info_dict['output_dir'] = output_dir
            info_dict['read_length'] = config["read_length"]
            info_dict['namespace'] = namespace
            info_dict['depth'] = config["%s_sim_depth" % namespace]
            info_dict['nucl_type'] = sim_info["nucl_type"]
            simulation_info.append(info_dict)

        pipeline_input_lines = []
        p = mp.Pool(config['ncpus_sim'])
        results = p.imap(simulate_genomes, simulation_info)
        for res in results:
            pipeline_input_lines.extend(res)
        p.close()
        p.join()

        # run pipeline, process 100 at a time
        for tmp_input in itertools.zip_longest(*[iter(pipeline_input_lines)]*100):
            pipeline_input_file = "%s/%s.input.txt" % (output_dir, namespace)
            out = open(pipeline_input_file, 'w')
            out.write("%s\n" % ("\n".join([x for x in tmp_input if x is not None])))
            overall_out_fh.write("%s\n" % ("\n".join([x for x in tmp_input if x is not None])))
            out.close()

            command = "%s %s %s %s --keep_dxsm --compress_output 1> %s_pipeline.out 2> %s_pipeline.err" % (config['python_interpreter'], config['pipeline_script'], pipeline_input_file, config['pipeline_config_file'], namespace, namespace)
            sys.stderr.write("%s\n" % command)
            sys.stderr.flush()
            sp.call(command, shell=True)

    else:
        sys.stderr.write("FATAL: %s not in ['16S', '18S', 'viral', 'bg', 'fg', 'pg']\n" % namespace)
        sys.stderr.flush()
        sys.exit(1)


inputs = {}  # {namespace}->{reporting id}->{'reporting_name':, 'taxid':, 'nucl_type':}
for line in open(args.input_file, 'r'):
    data = line.strip().split("\t")
    namespace = data[0]
    if namespace not in ['16S', '18S', 'viral', 'bg', 'fg', 'pg']:
        namespace = data[0].upper()
        if namespace not in ['16S', '18S', 'viral', 'bg', 'fg', 'pg']:
            sys.stderr.write("FATAL: %s not in ['16S', '18S', 'viral', 'bg', 'fg', 'pg']\n" % namespace)
            sys.stderr.flush()
            sys.exit(1)
    repid_dict = inputs.setdefault(namespace, {})
    repid_dict[data[1]] = {'reporting_name': data[2], 'taxid': int(data[3]), 'nucl_type': data[4].lower()}

for namespace, org_info in inputs.items():
    namespace_output_dir = "%s/%s" % (args.output_dir.rstrip("/"), namespace)
    try:
        os.mkdir(namespace_output_dir)
    except FileExistsError:
        pass
    sys.stderr.write("processing namespace %s\n" % namespace)
    sys.stderr.flush()
    overall_out_fh = open("%s/%s.all_input.txt" % (args.output_dir.rstrip("/"), namespace), 'w')
    start = time.time()
    drive_sims(namespace, org_info, args.output_dir.rstrip("/"), config, overall_out_fh)
    sys.stderr.write("TIME to process %s\n" % namespace)
    sys.stderr.flush()
    overall_out_fh.close()