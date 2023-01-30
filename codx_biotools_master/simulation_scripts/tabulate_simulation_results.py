import sys
import glob
import argparse
import json
import multiprocessing as mp
import copy

parser = argparse.ArgumentParser(description="Tabulate bg simulation results")
parser.add_argument("simulation_output_dir", type=str, help="Simulation output directory, complete path to top level")
parser.add_argument("tgs_file", type=str, help="tgs file")
parser.add_argument("ncbi_module_path", type=str, help="path to ncbi module code")
parser.add_argument("ncbi_path", type=str, help="folder to merged.dmp, nodes.dmp and names.dmp files")
parser.add_argument("--reporting_info", type=str, help="complete path to reporting info json", default=None)
parser.add_argument("--genome_file_map", type=str, help="map from simulation id to reference used to simulate", default=None)  # only use with bg namespace
parser.add_argument("--assembly_summary", type=str, help="ncbi assembly summary file of references used to simulate", default=None)  # only use with bg namespace
parser.add_argument("--no_group", action='store_true', help="set to not group by geneid when computing coverage", default=False)
parser.add_argument("--n_processes", type=int, help="number of processes used to parse dxsm output", default=5)
args = parser.parse_args()


reporting_obj = None
if args.reporting_info is not None:
    reporting_obj = json.load(open(args.reporting_info, 'r'))


def process_dxsm_output(simulation_info):  # [sim_folder, sim_taxid, dxsm, taxid_geneid_coverage, ncbi_tax, seqid_taxids, seqid_geneid, seqid_group]
    sim_folder = simulation_info[0]
    sim_taxid = simulation_info[1]
    dxsm = simulation_info[2]
    taxid_geneid_coverage = simulation_info[3]
    ncbi_tax = simulation_info[4]
    seqid_taxids = simulation_info[5]
    seqid_geneid = simulation_info[6]
    seqid_group = simulation_info[7]

    for taxid, geneid_info in taxid_geneid_coverage.items():
        for gid, coverage_info in geneid_info.items():
            coverage_info['coverage'] = 0.0
            coverage_info['n_covered_bases'] = 0

    if ncbi_tax.get_species_taxid_if_exists(sim_taxid) not in taxid_geneid_coverage:
        sys.stdout.write("WARNING: No database sequences, will skip for %d\t%s\n" % (sim_taxid, ncbi_tax.get_name(sim_taxid)))
        return sim_folder, sim_taxid, 0, {}

    sim_group_info = set()
    obj = json.load(open(dxsm, 'r'))
    for seqinfo in obj["sequence_coverage_information"]:
        taxid = seqid_taxids[seqinfo['sequence_name']]
        geneid = seqid_geneid[seqinfo['sequence_name']]
        n_covered_bases = seqinfo['n_covered_bases']
        sequence_length = seqinfo['sequence_length']
        coverage = seqinfo['coverage']
        geneid_info = taxid_geneid_coverage[taxid]
        coverage_info = geneid_info[geneid]

        if taxid == sim_taxid:  # add to sim_group_info
            sim_group_info.update(seqid_group.get(seqinfo["sequence_name"], []))

        if coverage > coverage_info['coverage']:
            coverage_info['coverage'] = coverage
            coverage_info['n_covered_bases'] = n_covered_bases
            coverage_info['sequence_length'] = sequence_length
    taxid_coverage = {}  # {taxid}->[n_covered_bases, total_sequence_length]
    for taxid, geneid_info in taxid_geneid_coverage.items():
        coverage_info = taxid_coverage.setdefault(taxid, [0, 0])
        for gid, g_coverage_info in geneid_info.items():
            coverage_info[0] += g_coverage_info['n_covered_bases']
            coverage_info[1] += g_coverage_info['sequence_length']

    highest_coverage = 0
    nonspecific_coverages = {}  # {taxid}->coverage
    for taxid, coverage_info in sorted(taxid_coverage.items(), key=lambda x: (x[1][0] / x[1][1]), reverse=True):
        if coverage_info[0] < 1:
            break
        coverage = coverage_info[0] / coverage_info[1]
        if taxid == sim_taxid or taxid in sim_group_info:
            if coverage > highest_coverage:
                highest_coverage = coverage
        else:
            nc = nonspecific_coverages.get(taxid, 0.0)
            if coverage > nc:  # keep highest coverage for each nonspecific taxid
                nonspecific_coverages[taxid] = coverage
    return sim_folder, sim_taxid, highest_coverage, nonspecific_coverages


def process_dxsm_output_no_group(simulation_info):  # [sim_folder, sim_taxid, dxsm, taxid_geneid_coverage, ncbi_tax, seqid_taxids, seqid_geneid]
    sim_folder = simulation_info[0]
    sim_taxid = simulation_info[1]
    dxsm = simulation_info[2]
    taxid_geneid_coverage = simulation_info[3]
    ncbi_tax = simulation_info[4]
    seqid_taxids = simulation_info[5]
    seqid_geneid = simulation_info[6]

    taxid_coverage = {}  # taxid->[n covered bases, sequence length, coverage]

    if ncbi_tax.get_species_taxid_if_exists(sim_taxid) not in taxid_geneid_coverage:
        sys.stdout.write("WARNING: No database sequences, will skip for %d\t%s\n" % (sim_taxid, ncbi_tax.get_name(sim_taxid)))
        return sim_folder, sim_taxid, 0, {}

    obj = json.load(open(dxsm, 'r'))
    for seqinfo in obj["sequence_coverage_information"]:
        taxid = seqid_taxids[seqinfo['sequence_name']]
        geneid = seqid_geneid[seqinfo['sequence_name']]
        n_covered_bases = seqinfo['n_covered_bases']
        sequence_length = seqinfo['sequence_length']
        coverage = seqinfo['coverage']
        prev_coverage = taxid_coverage.get(taxid, [0, 0, 0])
        if coverage > prev_coverage[2]:
            taxid_coverage[taxid] = [n_covered_bases, sequence_length, coverage]

    highest_coverage = 0
    nonspecific_coverages = {}  # {taxid}->coverage
    for taxid, coverage_info in sorted(taxid_coverage.items(), key=lambda x: (x[1][0] / x[1][1]), reverse=True):
        if coverage_info[0] < 1:
            break
        coverage = coverage_info[2]
        if taxid == sim_taxid:
            if coverage > highest_coverage:
                highest_coverage = coverage
        else:
            nc = nonspecific_coverages.get(taxid, 0.0)
            if coverage > nc:  # keep highest coverage for each nonspecific taxid
                nonspecific_coverages[taxid] = coverage
    return sim_folder, sim_taxid, highest_coverage, nonspecific_coverages

sys.path.append(args.ncbi_module_path)
from ncbi_taxonomy_utils import *
ncbi_tax = ncbi_taxonomy(args.ncbi_path.rstrip("/") + "/merged.dmp", args.ncbi_path.rstrip("/") + "/nodes.dmp", args.ncbi_path.rstrip("/") + "/names.dmp")

# "group": {"id": "0c6eba84d676a4f13f34ee9dc7ded66d", "taxid": 351195, "name": "Salimicrobium", "taxonomy": "351195:Salimicrobium;186817:Bacillaceae;1385:Bacillales;91061:Bacilli;1239:Firmicutes;1783272:Terrabacteria group;2:Bacteria;131567:cellular organisms;1:root", "members": [{"name": "Salimicrobium album", "taxid": 50717, "taxono
seqid_taxids = {}
seqid_geneid = {}
taxid_geneid_coverage = {}  # {taxid}->{geneid}->{coverage, n_covered_bases, sequence_length}
seqid_group = {}  # {seqid}->set(group member taxids)
for line in open(args.tgs_file):
    data = line.strip().split("\t")
    taxid = ncbi_tax.get_species_taxid_if_exists(int(data[0]))
    seqid_taxids[int(data[2])] = taxid
    seqid_geneid[int(data[2])] = int(data[1])
    obj = json.loads(data[3])
    geneid_info = taxid_geneid_coverage.setdefault(taxid, {})
    geneid_info.setdefault(int(data[1]), {'coverage': 0.0, 'n_covered_bases': 0, 'sequence_length': obj['seq_length']})
    if "group" in obj:
        for member_obj in obj.get("members", []):
            seqid_group.setdefault(int(data[2]), set([])).add(member_obj["taxid"])

accession_taxid = {}
if args.assembly_summary is not None:
    for line in open(args.assembly_summary):
        if line[0] == "#":  # skip comment lines
            continue
        data = line.strip().split("\t")
        accession = data[0]
        taxid = int(data[6])
        accession_taxid[accession] = taxid

sim_folder_taxid = {}
if args.genome_file_map is not None:
    for line in open(args.genome_file_map):
        data = line.strip().split("\t")
        accession = "_".join(data[1].split("/")[-1].split("_")[:2])
        sim_folder_taxid[data[0]] = accession_taxid[accession]

accumulated_results = {}  # {expected taxid}->{'correct_coverages':[], 'nonspecific_coverages':[]}
dxsm_files = glob.glob(args.simulation_output_dir.rstrip("/") + "/**/*.dxsm.out")
simulation_jobs = []  # [[sim_folder, sim_taxid, dxsm, taxid_geneid_coverage, ncbi_tax, seqid_taxids, seqid_geneid]]
for dxsm in dxsm_files:
    sim_folder = dxsm.split("/")[-2]
    sim_taxid = 0
    if args.genome_file_map is not None:  # if the folder doesn't map using the genome_file_map, then it is the seqid
        sim_taxid = sim_folder_taxid[sim_folder]
    else:
        sim_taxid = seqid_taxids[int(sim_folder)]
    job_info = [sim_folder,
                sim_taxid,
                dxsm,
                taxid_geneid_coverage,
                ncbi_tax,
                seqid_taxids,
                seqid_geneid,
                seqid_group
            ]
    simulation_jobs.append(job_info)

coverages = {}  # {taxid}->{'positive_coverage':[], 'nonspecific_coverage':[]}
p = mp.Pool(args.n_processes)
results = []
if args.no_group is False:
    results = p.imap(process_dxsm_output, simulation_jobs, 300)
else:
    results = p.imap(process_dxsm_output_no_group, simulation_jobs, 300)

for res in results:
    sim_folder = res[0]
    sim_taxid = res[1]
    coverage = res[2]
    nonspecific_coverages = res[3]
    coverages.setdefault(sim_taxid, {'positive_coverage': [], 'nonspecific_coverage': []})['positive_coverage'].append(coverage)
    for taxid, nspc in nonspecific_coverages.items():
        coverages.setdefault(taxid, {'positive_coverage': [], 'nonspecific_coverage': []})['nonspecific_coverage'].append(nspc)
    sys.stderr.write("Done with %s\n" % sim_folder)
    sys.stderr.flush()

open("all_coverages.json", 'w').write("%s" % (json.dumps(coverages)))


# gather cutoffs and compute true positives / false positives
out = open("sim_results_by_taxid", 'w')
for taxid, coverage_info in coverages.items():
    nonspecific_coverages = sorted(coverage_info['nonspecific_coverage'], reverse=True)
    tmp_org_rep_info = reporting_obj.get("%d" % taxid, {})
    cutoff = tmp_org_rep_info.get("RNA_stringent", 0.0)
    if args.genome_file_map is not None:
        cutoff = tmp_org_rep_info.get("DNA_stringent", 0.0)
    if reporting_obj is None and len(nonspecific_coverages) > 0:
        cutoff = nonspecific_coverages[0]  # choose highest nonspecific coverage as cutoff
    TP = 0
    FN = 0
    for cov in coverage_info['positive_coverage']:
        if cov > cutoff:
            TP += 1
        else:
            FN += 1
    total = TP + FN
    if total < 1:
        total = 1
    out.write("%d\t%f\t%d\t%d\t%f\n" % (taxid, TP / total, TP, FN, cutoff))
out.close()

