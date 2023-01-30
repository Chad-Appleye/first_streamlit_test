import sys
import glob
import argparse
import os
import json
import gzip
import multiprocessing as mp

"""
Calculates performance metrics of the simulations. 

Note:  True negatives are not included from other namespaces other than the intended one.  Reporting ids are excluded 
that were part of the expected taxid group
"""



parser = argparse.ArgumentParser(description="Drive simulations of organisms in input file")
parser.add_argument("input_file", type=str, help="input file where each line is <expected taxid><tab><summary output (with standard file naming)>")
parser.add_argument("pipeline_json", type=str, help="json file the pipeline used to produce simulation results")
parser.add_argument("profile_json", type=str, help="json file with profile organisms")
parser.add_argument("output_file", type=str, help="output file name, give complete path")
parser.add_argument("seq_type", type=str, choices=['RNA', 'DNA'], help="Has to be in ['RNA', 'DNA'] to be valid.  Evaluation only proceeds for either RNA or DNA, can't do both at the same time.  Input file can contain both, but will filter out disagreeing input\n")
parser.add_argument("--ncpus", type=int, help="number of summary files to process in parallel. default=1\n", default=1)
args = parser.parse_args()

config_args = json.load(open(args.pipeline_json, 'r'))

# change working directory to top level that contains all the database paths
os.chdir(config_args["top_level_dir"])

# taxid_reporting_name_info = json.load(open(config_args["taxid_reporting_names"], 'r'))

"""
load profile information
"""
profile_org_repids = {}  #{namespace}->{reporting ids}
profile_taxid_info = {}  # {taxid}->reporting info
profile_reporting_id_info = {}  # {repid}->reporting info
medically_relevant_ids = set()  # repid
for line in open(args.profile_json):
    obj = json.loads(line.strip())
    if obj["mr_status"] is True:
        medically_relevant_ids.add(obj["compound_id"])
    for taxid in obj["taxids"]:
        profile_taxid_info[taxid] = obj
    profile_reporting_id_info[obj["compound_id"]] = obj
    dna_info = obj.get("DNA", "NA").upper()
    if dna_info != "NA":
        if dna_info == 'BG':
            profile_org_repids.setdefault("bg", set()).add(obj["compound_id"])
        elif dna_info == 'FG':
            profile_org_repids.setdefault("fg", set()).add(obj["compound_id"])
        elif dna_info == 'PG':
            profile_org_repids.setdefault("pg", set()).add(obj["compound_id"])
        elif dna_info == 'VIRAL' or obj.get("class_type", None) == "viral":
            profile_org_repids.setdefault("viral", set()).add(obj["compound_id"])
        else:
            sys.stderr.write("WARNING:\tNo dna namespace found for %s\n" % (line.strip()))
            sys.stderr.flush()
    rna_info = obj.get("RNA", "NA").upper()
    if rna_info != "NA":
        if rna_info == 'VIRAL' or obj.get("class_type", None) == "viral":
            profile_org_repids.setdefault("viral", set()).add(obj["compound_id"])
        elif rna_info == '16S':
            profile_org_repids.setdefault("16S", set()).add(obj["compound_id"])
        elif rna_info == '18S':
            profile_org_repids.setdefault("18S", set()).add(obj["compound_id"])
        else:
            sys.stderr.write("WARNING:\tNo rna namespace found for %s\n" % (line.strip()))
            sys.stderr.flush()
    if obj.get("class_type", None) == "viral":
        profile_org_repids.setdefault("viral", set()).add(obj["compound_id"])


def get_reporting_info(taxid):
    reporting_info = None
    for tx in ncbi_tax.get_path(taxid):
        reporting_info = profile_taxid_info.get(tx, None)
        if reporting_info is not None:
            break
    return reporting_info

def get_reporting_info_by_repid(repid):
    return profile_reporting_id_info.get(repid, None)


def taxid_in_path(tx1, tx2):
    """
    Check to see if taxid 2 is in path of taxid 1
    :param tx1:  Taxid
    :param tx2:  Taxid
    :return: True or False
    """
    for tx in ncbi_tax.get_path(tx1):
        if tx == tx2:
            return True
    return False

"""
load ncbi taxonomy
"""
sys.path.append(config_args["ncbi_tax_utils"])
from ncbi_taxonomy_utils import *
ncbi_taxonomy_path = config_args["ncbi_tax_path"].rstrip("/")
ncbi_tax = ncbi_taxonomy("%s/merged.dmp" % ncbi_taxonomy_path, "%s/nodes.dmp" % ncbi_taxonomy_path, "%s/names.dmp" % ncbi_taxonomy_path)


# summary keys:
# 'total_genes', 'name', 'version', 'class_type', 'ncbi_name', 'taxonomy', 'taxid', 'group_coverage', 'gene_info',
# 'compound_id', 'core_coverage', 'dbs', 'reporting_id', 'total_bases', 'read_count', 'mean_consensus_pid', 'overall_covered_bases',
#  'coverage', 'subclass', 'mr_status'
# group keys:


def parse_evaluate_summary_out(info_dict, expected_reporting_info, expected_sim=True):
    # open file and set variables relating to gzipped or not
    tp = {}
    fp = {}
    fn = {}
    tn = {}
    disregarded_tn_repids = {}
    gz = False
    if os.path.isfile(info_dict["summary_out"]) is False:
        sys.stderr.write("%s doesn't exist, will skip\n" % info_dict["summary_out"])
        sys.stderr.flush()
        if expected_sim is True:
            fn[expected_reporting_info["compound_id"]] = {"coverage": 0.0}
        return tp, fp, fn, tn, disregarded_tn_repids
    fh = open(info_dict["summary_out"])
    if ".gz" in info_dict["summary_out"]:
        fh = gzip.open(info_dict["summary_out"])
        gz = True
    disregard_group_repid = set()
    expected_coverage = 0.0
    # first pass through summary file and find repids to exclude from performance calculations.  Also find expected coverage
    for line in fh:
        tmp_line = line
        if gz is True:
            tmp_line = line.decode()
        obj = json.loads(tmp_line.strip())
        expected_org_group = False
        if obj["reporting_id"] == expected_reporting_info["compound_id"]:  # if organism is the one that is expected
            expected_org_group = True
            expected_coverage = obj["coverage"]
        if expected_org_group is False:  # go through group to see is expected organism is a group member
            for tmp_goi in obj.get("group", {}).get("members", []):
                if tmp_goi["reporting_id"] == expected_reporting_info["compound_id"]:
                    expected_org_group = True
                    break
        if expected_org_group is True:  # add all reporting ids to be excluded from TN calculations
            disregard_group_repid.add(obj["reporting_id"])
            for tmp_goi in obj.get("group", {}).get("members", []):
                disregard_group_repid.add(tmp_goi["reporting_id"])
    fh.close()

    # second pass is to collect performance information
    found_expected = False
    fh = open(info_dict["summary_out"])
    if ".gz" in info_dict["summary_out"]:
        fh = gzip.open(info_dict["summary_out"])
        gz = True
    for line in fh:
        tmp_line = line
        if gz is True:
            tmp_line = line.decode()
        obj = json.loads(tmp_line.strip())
        if obj["reporting_id"] not in profile_reporting_id_info:  # do not count anything that is not part of the profile
            continue
        if obj["reporting_id"] not in medically_relevant_ids:  # skip any ids that are not medically relevant
            continue
        if obj["reporting_id"] == expected_reporting_info["compound_id"]:  # if organism is the one that is expected
            found_expected = True
            cutoff = expected_reporting_info["%s_stringent" % info_dict["seq_type"]]
            if obj["coverage"] >= cutoff:  # true positive
                tp[obj["reporting_id"]] = {"coverage": obj["coverage"]}
            else:  # false negative
                fn[obj["reporting_id"]] = {"coverage": obj["coverage"]}
        else:  # possible false positive or true negative
            if obj["reporting_id"] in disregard_group_repid:  # skip any members in group of expected organism
                continue
            rep_info = get_reporting_info(obj["taxid"])
            if rep_info is None:
                sys.stderr.write("WARNING: no reporting info for taxid %d, will skip\n"%(obj["taxid"]))
                sys.stderr.flush()
                continue
            ck = "%s_stringent" % info_dict["seq_type"]
            cutoff = rep_info.get(ck, None)
            if cutoff is None:
                sys.stderr.write("FATAL: key %s doesn't exist in %s\n" % (ck, rep_info))
                sys.stderr.flush()
                sys.exit(1)
            if obj["coverage"] >= cutoff:  # false positive
                fp[obj["reporting_id"]] = {"coverage": obj["coverage"],
                                           "namespace": info_dict["namespace"],
                                           "expected_class_type": expected_reporting_info.get("class_type", "NA"),
                                           "expected_nucleic_acid": expected_reporting_info.get("nucleic_acid", "NA"),
                                           "expected_specific_threshold": expected_reporting_info["%s_stringent" % expected_reporting_info["seq_type"]],
                                           "expected_reporting_id": expected_reporting_info["compound_id"],
                                           "expected_taxid": expected_reporting_info["taxids"][0],
                                           "expected_coverage": expected_coverage,
                                           "expected_name": expected_reporting_info["reporting_name"],
                                           "taxid": obj["taxid"],
                                           "name": obj["name"]}
            else:  # true negative
                tn[obj["reporting_id"]] = {"coverage": obj["coverage"], "name": obj["name"], "taxid": obj["taxid"]}
    fh.close()
    if found_expected is False:
        fn[expected_reporting_info["compound_id"]] = {"coverage": expected_coverage}

    for repid in profile_org_repids.get(info_dict["namespace"], []):
        if repid not in medically_relevant_ids:  # skip any non medically relevant organisms
            continue
        if repid in tp:
            pass
        elif repid in fp:
            pass
        elif repid in tn:
            pass
        elif repid in disregard_group_repid:
            rep_info = get_reporting_info_by_repid(repid)
            taxid = -1  # default taxid value, need to check for missing taxids
            if len(rep_info["taxids"]) > 0:
                taxid = rep_info["taxids"][0]
            disregarded_tn_repids[repid] = {"name": rep_info["reporting_name"], "taxid": taxid}
        else:
            rep_info = get_reporting_info_by_repid(repid)
            taxid = -1  # default taxid value, need to check for missing taxids
            if len(rep_info["taxids"]) > 0:
                taxid = rep_info["taxids"][0]
            tn[repid] = {"coverage": 0, "name": rep_info["reporting_name"], "taxid": taxid}

    return tp, fp, fn, tn, disregarded_tn_repids


def evaluate_summary_output(info_dict):
    """
    :param info_dict:  {'summary_out': summary_file,
                        'accession': file accession,
                        'expected_taxid': taxid,
                        'tgs_file': tgs_file,
                        'namespace': namespace,
                        'seq_type': seq_type,
                        'expected_rep_info': expected_reporting_info}
    :return: info_dict, true_positives, false_positives, false_negatives, disregarded_taxa
    """

    #true_positive = {}  # {reporting_id}->{coverage:, taxid:, name:}
    #true_negatives = {}  # {reporting_id}->{coverage:, taxid:, name:}
    #false_positives = {}  # {reporting_id}->{coverage:, taxid:, name:}
    #false_negative = {}  # {reporting_id}->{coverage:, taxid:, name:}

    true_positives, false_positives, false_negative, true_negatives, disregarded_repids = parse_evaluate_summary_out(info_dict, info_dict['expected_rep_info'])
    # check for cross namespace false positives
    namespace_postfix = {
        "RNA":{"16S": ["rna.bacterial.dxsm.out.summary", config_args["bact_16S_tgs"]],
               "18S": ["rna.fungal_parasite.dxsm.out.summary", config_args["fp_18S_tgs"]],
               "viral": ["rna.viral.dxsm.out.summary", config_args["viral_tgs"]]
               },
        "DNA":{"bg": ["dna.bacterial.dxsm.out.summary", config_args["bact_gen_tgs"]],
               "fg": ["dna.fungal_parasite.fungal.dxsm.out.summary", "None"],
               "pg": ["dna.fungal_parasite.parasite.dxsm.out.summary", "None"],
               "viral": ["dna.viral.dxsm.out.summary", config_args["viral_tgs"]]
               }
    }
    for seq_type, file_endings in namespace_postfix.items():  # find cross namespace FP
        if info_dict["seq_type"] == seq_type:
            for tmp_namespace, file_ending in file_endings.items():
                if tmp_namespace != info_dict["namespace"]:
                    tmp_info_dict = {
                        'summary_out': info_dict['accession'] + "." + file_ending[0],
                        'accession': info_dict['accession'],
                        'expected_taxid': info_dict["expected_taxid"],
                        'tgs_file': file_ending[1],
                        'namespace': tmp_namespace,
                        'seq_type': seq_type
                    }
                    tp, fp, fn, tn, dr = parse_evaluate_summary_out(tmp_info_dict, info_dict['expected_rep_info'], expected_sim=False)
                    for repid, info in fp.items():
                        false_positives[repid] = info

    return info_dict, true_positives, false_positives, false_negative, true_negatives, disregarded_repids


summary_info_input = []
for line in open(args.input_file, 'r'):
    data = line.strip().split("\t")
    summary_file = data[1]
    taxid = int(data[0])
    namespace = None
    tgs_file = None
    accession = None
    composition_file = None
    seq_type = None

    expected_reporting_info = get_reporting_info(taxid)
    if expected_reporting_info is None:
        sys.stderr.write("WARNING: no reporting info found for taxid %d, will skip\n" % (taxid))
        sys.stderr.flush()
        continue

    if expected_reporting_info["mr_status"] is False:  # skip anything that is not medically relevant
        continue

    if "rna.bacterial.dxsm.out.summary" in summary_file:
        namespace = "16S"
        tgs_file = config_args["bact_16S_tgs"]
        accession = ".".join(summary_file.split(".")[:-5])
        seq_type = "RNA"
    elif "rna.fungal_parasite.dxsm.out.summary" in summary_file:
        namespace = "18S"
        tgs_file = config_args["fp_18S_tgs"]
        accession = ".".join(summary_file.split(".")[:-5])
        seq_type = "RNA"
    elif "rna.viral.dxsm.out.summary" in summary_file:
        namespace = "viral"
        tgs_file = config_args["viral_tgs"]
        accession = ".".join(summary_file.split(".")[:-5])
        seq_type = expected_reporting_info.get("nucleic_acid", "NA")
    elif "dna.viral.dxsm.out.summary" in summary_file:
        namespace = "viral"
        tgs_file = config_args["viral_tgs"]
        accession = ".".join(summary_file.split(".")[:-5])
        seq_type = expected_reporting_info.get("nucleic_acid", "NA")
    elif "dna.bacterial.dxsm.out.summary" in summary_file:
        namespace = "bg"
        tgs_file = config_args["bact_gen_tgs"]
        accession = ".".join(summary_file.split(".")[:-5])
        seq_type = "DNA"
    elif "dna.fungal_parasite.parasite.dxsm.out.summary" in summary_file:
        tgs_file = "None"
        namespace = "pg"
        accession = ".".join(summary_file.split(".")[:-6])
        seq_type = "DNA"
    elif "dna.fungal_parasite.fungal.dxsm.out.summary" in summary_file:
        tgs_file = "None"
        namespace = "fg"
        accession = ".".join(summary_file.split(".")[:-6])
        seq_type = "DNA"

    if namespace is None:  # skip if no namespace is found
        sys.stderr.write("WARNING: no namespace detected for %s, will skip\n" % (summary_file))
        sys.stderr.flush()
        continue
    elif seq_type != args.seq_type:  # skip if the seq type is incorrect
        sys.stderr.write("INCORRECT seq_type for %s, will skip\n" % summary_file)
        sys.stderr.flush()
        continue
    else:
        expected_reporting_info["seq_type"] = seq_type
        info_obj = {
            'summary_out': summary_file,
            'accession': accession,
            'expected_taxid': taxid,
            'tgs_file': tgs_file,
            'namespace': namespace,
            'seq_type': seq_type,
            'expected_rep_info': expected_reporting_info
        }
        summary_info_input.append(info_obj)


results_json = {}  # {reporting_id}->{namespace}->{name:, taxid:, TP:, FP:, FN:, TN:, DR:}
                   # {reporting_id}->{meta}->{namespace:, class_type:, nucleic_acid:, specific_threshold:, sensitive_threshold:}
p = mp.Pool(args.ncpus)
results = p.imap(evaluate_summary_output, summary_info_input, 1)
for res in results:
    info_dict = res[0]
    expected_reporting_info = info_dict["expected_rep_info"]
    tp = res[1]
    fp = res[2]
    fn = res[3]
    tn = res[4]
    disregarded_repids = res[5]
    expected_coverage = 0.0

    # add meta info
    namespace_res = results_json.setdefault(expected_reporting_info["compound_id"], {})
    namespace_res.setdefault("meta", {
        "namespace": info_dict["namespace"],
        "class_type": expected_reporting_info.get("class_type", "NA"),
        "nucleic_acid": expected_reporting_info.get("nucleic_acid", "NA"),
        "specific_threshold": expected_reporting_info["%s_stringent" % expected_reporting_info["seq_type"]],
        "sensitive_threshold": expected_reporting_info["%s_gray" % expected_reporting_info["seq_type"]]
    })

    # True positives
    for repid, info in tp.items():
        namespace_res = results_json.setdefault(repid, {})
        namespace_res.setdefault(info_dict["namespace"], {"name": expected_reporting_info["reporting_name"], "taxid": expected_reporting_info["taxids"][0], "TP": [], "FP":[], "FN":[], "TN": 0, "DR": 0})["TP"].append(info["coverage"])
        expected_coverage = info["coverage"]
    # False negatives
    for repid, info in fn.items():
        namespace_res = results_json.setdefault(repid, {})
        namespace_res.setdefault(info_dict["namespace"], {"name": expected_reporting_info["reporting_name"], "taxid": expected_reporting_info["taxids"][0], "TP": [], "FP":[], "FN":[], "TN": 0, "DR": 0})["FN"].append(info["coverage"])
        expected_coverage = info["coverage"]
    # False positives
    for repid, info in fp.items():
        namespace_res = results_json.setdefault(repid, {})
        res_dict = namespace_res.setdefault(info["namespace"], {"name": info["name"], "taxid": info["taxid"], "TP": [], "FP":[], "FN": [], "TN": 0, "DR": 0})
        res_dict["FP"].append({
            "coverage": info["coverage"],
            "expected_coverage": expected_coverage,
            "expected_class_type": info["expected_class_type"],
            "expected_nucleic_acid": info["expected_nucleic_acid"],
            "expected_specific_threshold": info["expected_specific_threshold"],
            "expected_reporting_id": info["expected_reporting_id"],
            "expected_namespace": info_dict["namespace"],
            "expected_taxid": info["expected_taxid"],
            "expected_name": info["expected_name"],
            "summary_out": info_dict["summary_out"]
        })
    # True negatives
    for repid, info in tn.items():
        namespace_res = results_json.setdefault(repid, {})
        namespace_res.setdefault(info_dict["namespace"], {"name": info["name"], "taxid": info["taxid"], "TP": [], "FP":[], "FN":[], "TN": 0, "DR": 0})["TN"] += 1
    # disregarded repids
    for repid, info in disregarded_repids.items():
        namespace_res = results_json.setdefault(repid, {})
        namespace_res.setdefault(info_dict["namespace"], {"name": info["name"], "taxid": info["taxid"], "TP": [], "FP":[], "FN":[], "TN": 0, "DR": 0})["DR"] += 1

json.dump(results_json, open(args.output_file, 'w'))