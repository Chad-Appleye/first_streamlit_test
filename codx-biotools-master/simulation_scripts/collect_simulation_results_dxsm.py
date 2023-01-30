import sys
import glob
import multiprocessing as mp
import json
import argparse
import os

parser = argparse.ArgumentParser(description="Collect simulation results from results of dxsm output")
parser.add_argument("simulation_output_dir", type=str, help="Simulation output directory, complete path to top level")
parser.add_argument("config_file", type=str, help="config file name")
parser.add_argument("namespace", type=str, help="specify namespace.  can be one of [16S, 18S, viral, bg]")
args = parser.parse_args()

config_args = json.load(open(args.config_file, 'r'))

os.chdir(config_args["top_level_dir"])

"""
load ncbi taxonomy
"""
sys.path.append(config_args["ncbi_tax_utils"])
from ncbi_taxonomy_utils import *
ncbi_taxonomy_path = config_args["ncbi_tax_path"].rstrip("/")
ncbi_tax = ncbi_taxonomy("%s/merged.dmp" % ncbi_taxonomy_path, "%s/nodes.dmp" % ncbi_taxonomy_path, "%s/names.dmp" % ncbi_taxonomy_path)


def get_reporting_info(reporting_json, taxid):
    reporting_name_info = {}
    for tx in ncbi_tax.get_path(taxid):
        if "%d" % tx in reporting_json:
            reporting_name_info = reporting_json["%d" % tx]
            break
    return reporting_name_info


def get_group_member_taxids(tgs_obj):
    group_info = tgs_obj.get("group", None)
    member_taxids = []
    if group_info is not None:
        members_info = group_info.get("members", [])
        for mi in members_info:
            if mi.get("taxid", -1) > 2:
                member_taxids.append(mi["taxid"])
    return member_taxids


taxid_reporting_name_info = json.load(open(config_args["taxid_reporting_names"], 'r'))

tgs_file = None
if args.namespace == "bg":
    tgs_file = config_args["bact_gen_tgs"]
elif args.namespace == "16S":
    tgs_file = config_args["bact_16S_tgs"]
elif args.namespace == "18S":
    tgs_file = config_args["fp_18S_tgs"]
elif args.namespace == "viral":
    tgs_file = config_args["viral_tgs"]
else:
    sys.stderr.write("FATAL: no class type recognized for namespace %s\n" % args.namespace)
    sys.exit(1)

seqid_taxid_obj = {}  # {seqid}->[taxid, obj]
for line in open(tgs_file, 'r'):
    data = line.strip().split("\t")
    info = json.loads(data[3])
    seqid_taxid_obj[int(data[2])] = [int(data[0]), info]

dxsm_files = glob.glob("%s/**/*.dxsm.out" % (args.simulation_output_dir.rstrip("/")), recursive=True)
sys.stderr.write("Number of dxsm output files to process: %d\n" % len(dxsm_files))
# {"sequence_coverage_information":[{"sequence_name":50328, "n_reads":45, "seque

for dxsm in dxsm_files:
    seqid = int(dxsm.split("/")[-2])
    taxid_info = seqid_taxid_obj[seqid]
    reporting_info = get_reporting_info(taxid_reporting_name_info, taxid_info[0])
    obj = json.load(open(dxsm, 'r'))
    max_seqid = -1
    max_taxid_info = [-1, {}]
    max_coverage = 0.0
    max_nonspecific_coverage = 0.0  # only include if not part of a group or previous reporting id
    max_group_taxids = {}
    max_nonspecific_seqid = -1
    max_nonspecific_taxid_info = [0, None]
    max_nonspecific_reporting_info = {}
    for seq_info in obj["sequence_coverage_information"]:
        txi = seqid_taxid_obj[seq_info["sequence_name"]]
        ri = get_reporting_info(taxid_reporting_name_info, txi[0])
        if txi[0] == taxid_info[0] or ri.get("reporting_id") == reporting_info.get("reporting_id"):
            if seq_info["coverage"] > max_coverage:
                max_taxid_info = txi
                max_seqid = seq_info["sequence_name"]
                max_coverage = seq_info["coverage"]
                max_group_taxids = get_group_member_taxids(txi[1])
        else:
            if seq_info["coverage"] > max_nonspecific_coverage:
                tmp_seqid = seq_info["sequence_name"]
                tmp_txi = seqid_taxid_obj[tmp_seqid]
                if tmp_txi[0] not in max_group_taxids:
                    max_nonspecific_seqid = tmp_seqid
                    max_nonspecific_coverage = seq_info["coverage"]
                    max_nonspecific_taxid_info = tmp_txi
                    max_nonspecific_reporting_info = get_reporting_info(taxid_reporting_name_info, tmp_txi[0])
    """
    simulated seqid
    max correct seqid
    max correct taxid
    max correct reporting name
    max correct coverage
    max nonspecific seqid
    max nonspecific taxid
    max nonspecific reporting name
    max nonspecific coverage
    """
    sys.stdout.write("%d\t%d\t%d\t%s\t%f\t%d\t%d\t%s\t%f\n" % (
        seqid,
        max_seqid,
        max_taxid_info[0],
        reporting_info.get("reporting_name", ncbi_tax.get_name(max_taxid_info[0])),
        max_coverage,
        max_nonspecific_seqid,
        max_nonspecific_taxid_info[0],
        max_nonspecific_reporting_info.get("reporting_name", ncbi_tax.get_name(max_nonspecific_taxid_info[0])),
        max_nonspecific_coverage
        )
    )

