import sys
import argparse
import pysam
import random
from itertools import cycle
import numpy as np


def simulate_reads(sequence_name, sequence, depth, read_length):
    fasta_read_blocks = []
    n_reads = int((len(sequence) * depth) / read_length)
    if n_reads < depth:  # simulate at least depth number of reads
        n_reads = depth
    sim_end_position = (len(sequence) - read_length) + 1
    if sim_end_position < 0:  # check if sequence is less than read length, if so, then use 0 as the position
        sim_end_position = 0
    for i in range(n_reads):
        random_start = random.randint(0, sim_end_position)
        fasta_read_blocks.append(">%s.%d.%d\n%s\n" % (sequence_name, random_start, i, sequence[random_start:random_start+read_length]))
    return fasta_read_blocks


def simulate_reads_tiled(sequence_name, sequence, stride, read_length, cover_end=False):
    fasta_read_blocks = []
    i = 0
    read_count = 0
    intended_coverage = int(read_length / stride)
    if len(sequence) < (read_length + (stride * intended_coverage)):  # if sequence is too short to get coverage
        return simulate_reads(sequence_name, sequence, intended_coverage, read_length)
    while (i + read_length - 1) < len(sequence):
        fasta_read_blocks.append(">%s.%d.%d\n%s\n" % (sequence_name, i, read_count, sequence[i:i + read_length]))
        read_count += 1
        i += stride
    # add read to cover end of sequence if needed, and if cover_end is True
    if cover_end:
        j = len(sequence) - read_length
        if j % stride != 0:
            fasta_read_blocks.append(">%s.%d.%d\n%s\n" % (sequence_name, j, read_count, sequence[j:]))
    return fasta_read_blocks


def simulate_reads_tiled_fh(file_handle, sequence_name, sequence, stride, read_length, cover_end=False):
    i = 0
    read_count = 0
    intended_coverage = int(read_length / stride)
    if len(sequence) < (read_length + (stride * intended_coverage)):  # if sequence is too short to get coverage
        for read_block in simulate_reads(sequence_name, sequence, intended_coverage, read_length):
            file_handle.write(read_block)
    while (i + read_length - 1) < len(sequence):
        file_handle.write(">%s.%d.%d\n%s\n" % (sequence_name, i, read_count, sequence[i:i + read_length]))
        read_count += 1
        i += stride
    # add read to cover end of sequence if needed, and if cover_end is True
    if cover_end:
        j = len(sequence) - read_length
        if j % stride != 0:
            file_handle.write(">%s.%d.%d\n%s\n" % (sequence_name, j, read_count, sequence[j:]))


def simulate_reads_tiled_by_depth(sequence_name, sequence, depth, read_length, file_handle=None):
    """
    Create simulated reads evenly spaced across the entire sequence.
    When file handle is set to an existing writable file handle, the return value of fasta_read_blocks
    will be empty since the reads will be written to a file instead.  Useful for large sequences (genomes)
    """
    fasta_read_blocks = []
    num_reads = round(len(sequence) * depth / read_length)
    last_read_coord = len(sequence) - read_length
    stride = last_read_coord / num_reads
    if stride < 1: stride = 1
    coords = [int(x) for x in np.around(np.arange(0, last_read_coord, stride))] + [last_read_coord]
    if len(coords) < num_reads:
        start_coords = []
        for coord in cycle(coords):
            start_coords.append(coord)
            if len(start_coords) >= num_reads + 2: break
    else:
        start_coords = coords
    read_count = 0
    for sc in start_coords:
        if file_handle is None:
            fasta_read_blocks.append(">%s.%d.%d\n%s\n" % (sequence_name, sc, read_count, sequence[sc:sc + read_length]))
        else:
            file_handle.write(">%s.%d.%d\n%s\n" % (sequence_name, sc, read_count, sequence[sc:sc + read_length]))
        read_count += 1
    return fasta_read_blocks


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulate single-ended, perfect reads from at specified depth from each sequence in the input file.  Output is ><seq name>.<pos>.<read id>")
    parser.add_argument("input_file", type=str, help="Fasta file, can be zipped")
    parser.add_argument("--depth", type=int, help="Depth at which to simulate.  Default=10", default=10)
    parser.add_argument("--read_length", type=int, help="Length of simulated reads", default=150)
    args = parser.parse_args()

    for seq in pysam.FastxFile(args.input_file):
        sys.stdout.write("%s" % "".join(simulate_reads(seq.name, seq.sequence, args.depth, args.read_length)))
        