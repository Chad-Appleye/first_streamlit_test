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

parser = argparse.ArgumentParser(description="Minhash search results from binned reads")
parser.add_argument("input_binner", type=str, help="input binner file from txtll output, with reads named <seq id>.<pos>.<read id>")
parser.add_argument("minhash_config", type=str, help="minhash database")
parser.add_argument("n_processes", type=int, help="number of processes to run concurrently")
parser.add_argument("output_directory", type=str, help="output directory for simulations, will create seqid specific outputs under this")
parser.add_argument("parent_taxid", type=int, help="parent taxid under which to collect reads")
args = parser.parse_args()

"""
json config file

minhash_binary:         path to minhash searching executable
minhash_db:             path to minhash database
minhash_metainfo:       path to minhash database meta information
minhash_nkmers:         number of kmers to use from query sample
ncbi_tax_utils:         path to folder with ncbi_taxonomy_utils
ncbi_tax_path:          path to folder containing merged.dmp, nodes.dmp, names.dmp
"""

config_args = json.load(open(args.minhash_config, 'r'))

"""
load ncbi taxonomy
"""
sys.path.append(config_args["ncbi_tax_utils"])
from ncbi_taxonomy_utils import *
ncbi_taxonomy_path = config_args["ncbi_tax_path"].rstrip("/")
ncbi_tax = ncbi_taxonomy("%s/merged.dmp" % ncbi_taxonomy_path, "%s/nodes.dmp" % ncbi_taxonomy_path, "%s/names.dmp" % ncbi_taxonomy_path)

cached_taxid_paths = {}  # {taxid}->[path]
def get_taxid_path(taxid):
    path = cached_taxid_paths.get(taxid, None)
    if path is None:
        path = ncbi_tax.get_path(taxid)
        cached_taxid_paths[taxid] = path
    return path


def process_binner_block(block_info):
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
        for tx in get_taxid_path(taxid):
            taxid_propagated_count[tx] += 1
            if args.parent_taxid == tx:
                out_fa.write(">%s\t%d\n%s\n" % (data[0], taxid, data[6]))
    f.close()
    out_fa.close()



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

