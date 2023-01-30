import time
import argparse
import re
import os
import pickle

# Removing these imports as they are causing problems. 
# import Levenshtein
# from oligotools import *

import pandas as pd
import numpy as np
from pandas.core.algorithms import value_counts
import plotly.graph_objects as go

from Bio.SeqIO.FastaIO import SimpleFastaParser

from oligotools import reverse_complement as rc
from oligotools import expand_ambiguous_dna as expand, coprimers_to_dataframe
from tqdm import tqdm


def _get_fasta_stats(handle):
    """parse the fasta file
    :param handle: opened fasta file

    :output: number of records, record length
    """
    num_records = 0
    for _, seq in SimpleFastaParser(handle):
        num_records += 1
    record_len = len(seq)
    return num_records, record_len


def parse_alignment(fasta_path, return_array=False):
    """parse fasta alignment file
    :param fasta_path: path to fasta alignemnt file
    :param return_array: return both numpy array and list of lists

    :output: alignment sequences as list of lists, returns alignment as
    numpy array if return_array=True
    """
    if return_array:
        with open(fasta_path) as f:
            num_records, record_len = _get_fasta_stats(f)
        seq_array = np.empty((num_records, record_len), dtype=str)

    with open(fasta_path) as f:
        name_ls = []
        
        seq_ls = []
        idx = 0
        for name, seq in SimpleFastaParser(f):
            name_ls.append(name)
            if return_array:
                seq_array[idx, :] = list(seq)
            seq_ls.append(seq)
            idx += 1
    if return_array:
        return seq_ls, seq_array
    else:
        return seq_ls

    
def flatten(iterable):
    """ return a 1-D list for any 2-D iterable
    """
    return [j for i in iterable for j in i]


def find_insertion_deletion(sequence):
    """ Using Regex, find the '-' caracters in the sequence that correspond to insertion or deletion by identifying if 
    the '-' has sequenced bases before and after it
    
    Params:
    input: sequence (str)
    out: sequence with 'Z' in place of '-' which are insertions or deletions (str)
    """
    lstripped = sequence.lstrip('-')
    start = len(sequence) - len(lstripped)
    stripped = lstripped.rstrip('-')
    end = len(lstripped) - len(stripped)
    return '-' * start + stripped.replace('-', 'Z') + '-' * end


def mismatch_args(coprimer_path, alignment_path, sheet_name=0, pseudo_genome=None):
    """ 
    read excel or csv file to get the CoPrimer information to pass to MismatchVisualization class
    Params:
        coprimer_path: path to the csv or xlsx file containing the CoPrimer information from the CoPrimer design software. Should include columns: OligoName, Sequence, Gap
        alignment_path: path to the alignment (fasta.aln) file used in the design on the CoPrimers
        sheet_name: (Optional) if passing xlsx file, can specify sheet name. See pandas.read_excel() documentation for more. 
        pseudo_genome: (str, path, SeqIO.Record) (optional) -- if pseudo_genome is already determined from Entropy or other means, it can be passed here 
            which will speed up calculation time significantly
    Out:
        dictionary with kwargs that can be passed to the MismatchVisualization class using ** method of expansion
    """

    if coprimer_path.endswith('.csv') or coprimer_path.endswith('.txt'):
        print('reading csv file...')
        df = pd.read_csv(coprimer_path)
    elif coprimer_path.endswith('.xlsx'):
        print('reading excel...')
        df = coprimers_to_dataframe(coprimer_path, sheet_name=sheet_name)
    else:
        raise Exception('invalid file type passed')
    
    if pseudo_genome != None:
        #check if it's a file
        if os.path.isfile(pseudo_genome):
            if os.path.splitext(pseudo_genome)[1] == '.fasta':
                pseudo_genome = get_seq_from_fasta(pseudo_genome)
            else:
                with open(pseudo_genome, 'r') as f:
                    pseudo_genome = f.readline()
        elif isinstance(pseudo_genome, str):
            pass
        elif isinstance(pseudo_genome, Bio.SeqRecord.SeqRecord):
            pseudo_genome = str(pseudo_genome.seq)
        else:
            raise Exception(f'unrecognized type given for pseudo_genome. {type(pseudo_genome)} given. Must be path to file, string sequence, or Bio.SeqRecord')

    #split the forwards and reverses
    forwards, reverses = split_forwards_reverses(df)

    #build dict
    kwargs = {
        'alignment_path': alignment_path,
        'forward_primer_sequences': forwards.Sequence.tolist(),
        'reverse_primer_sequences': reverses.Sequence.tolist(), 
        'forward_primer_names': forwards.OligoName.tolist(), 
        'reverse_primer_names': reverses.OligoName.tolist(),
        'forward_gap_lengths': forwards.Gap.tolist(),
        'reverse_gap_lengths': reverses.Gap.tolist(), 
        'pseudo_genome': pseudo_genome
    }
    print('')
    return kwargs


class Primers:
    """A class to represent primers and an alignment
    Attributes
    ----------
    alignment: ls
        The input alignment
    primers: ls
        The input primer sequences
    expanded_primers: ls
        Every primer that is represented by the given ambigous sequence
    primer_regions: tuple
        (start, end) position of the primer relative to the alignment
    expanded_primer_regions: list of tuples
        same as primer regions, but as list with same length as 
        expanded_primers
    """
    def __init__(self, alignment, primers, orientation='forward',
            iscoprimer=True, coprimer_gaps=None):
        if not isinstance(primers, list):
            self.primers = [primers]
        else:
            self.primers = primers
        if iscoprimer:
            if coprimer_gaps is None:
                raise ValueError(
                    'Gap lengths must be supplied with CoPrimer sequences.')
            if not isinstance(coprimer_gaps, list):
                coprimer_gaps = [coprimer_gaps]
            self.primers = [self._linearize_coprimer(coprimer, gap) 
                            for coprimer, gap in 
                            zip(self.primers, coprimer_gaps)]
        if orientation == 'reverse_complement':
            self.primers = [rc(primer) for primer in self.primers]
        expanded_primers = []
        primer_regions = []
        expanded_primer_regions = []
        for primer in self.primers:
            current_primer_region = self._get_primer_region(alignment, primer)
            current_expanded_primers = expand(primer, ignore_n=True)
            primer_regions.append(current_primer_region)
            expanded_primers.extend(current_expanded_primers)
            expanded_primer_regions.extend([current_primer_region] 
                                            * len(current_expanded_primers))
        self.expanded_primers = expanded_primers
        self.primer_regions = primer_regions
        self.expanded_primer_regions = expanded_primer_regions

    def _get_primer_region(self, alignment, primer):
        start = -1
        primer_region = None
        for seq in alignment:
            for expanded_primer in expand(primer, ignore_n=True):
                pattern = re.compile(expanded_primer.replace('N', '.'))
                match = pattern.search(seq)
                if match:
                    start, end = match.span()
                    primer_region = (start, end)
                    break
                else:
                    continue
        if primer_region is None:
            raise Exception(f"primer {primer} not found in alignment")
        return primer_region

    def _linearize_coprimer(self, coprimer_seq, gap_len):
        seq_ls = [seq.split(']') for seq in coprimer_seq.split('[')]
        seq_ls_flat = [item.upper() for sublist in seq_ls for item in sublist]
        accepted_bases = {'A', 'C', 'G', 'T', 'R', 'Y', 'S', 'W'
                        'K', 'M', 'B', 'D', 'H', 'V', 'N'}
        primer_coprimer = [item for item in seq_ls_flat  
                        if set(item).issubset(accepted_bases) 
                        and not item == len(item) * 'A']
        linearized_coprimer = (primer_coprimer[1] + 'N' * gap_len 
            + primer_coprimer[0])
        return linearized_coprimer
    

class PrimerCoverage:
    """A class to represent primers and an alignment

    Attributes
    ----------
    alignment: ls
        The input alignment
    forward_primers: ls
        The input forward primer sequences
    reverse_primers: ls
        The input reverse primer sequences
    forward_names: ls
        the names of the forward primers. eg. 2631.F1
    reverse_names: ls
        the names of the reverse primers. eg. 2631.R1
    probes: ls
        The input probe sequences
    coverage: float
        The fraction of sequences for which all input primers/probes
        have an exact match
    matching_seqs: numpy array
        The sequences with a perfect match for all input primers/probes
        have an exact match
    nonmatching_seqs: numpy array
        The sequences not matching all input primers/probes
    nonmatching_forward_distances: numpy array
        The number of mismatches in the forward primers
    nonmatching_reverse_distances: numpy array
        The number of mismatches in the reverse primers
    nonmatching_probe_distances: numpy array
        The number of mismatches in the probes
    nonmatching_forward_mismatch_positions: numpy array
        The positions of the mismatches for the forward primers
    nonmatching_reverse_mismatch_positions: numpy array
        The positions of the mismatches for the reverse primers
    nonmatching_probe_mismatch_positions: numpy array
        The positions of the mismatches for the probes

    Methods
    -------
    display_report:
        displays an html report of the mismatch distribution
    """
    def __init__(self, alignment_path, forward_primers, reverse_primers, forward_names, reverse_names, probes=None,
                 probe_orientation='forward', iscoprimer=True, 
                 forward_gaps=None, reverse_gaps=None):
        self.forward_names = forward_names
        self.reverse_names = reverse_names
        self.alignment = parse_alignment(alignment_path)
        self.forward_primers = Primers(self.alignment, forward_primers, 
                                       iscoprimer=iscoprimer, 
                                       coprimer_gaps=forward_gaps)
        self.reverse_primers = Primers(self.alignment, reverse_primers,
                                       orientation='reverse_complement',
                                       iscoprimer=iscoprimer,
                                       coprimer_gaps=reverse_gaps)
        if probes is not None:
            if probe_orientation == 'reverse_complement':
                self.probes = Primers(self.alignment, probes, 
                                      orientation='reverse_complement',
                                      iscoprimer=False)
            else:
                self.probes = Primers(self.alignment, probes,
                                      iscoprimer=False)
        else:
            self.probes = None
        (self.coverage, self.matching_seqs, self.nonmatching_seqs, 
         self.nonmatching_forward_distances, self.nonmatching_reverse_distances, 
         self.nonmatching_probe_distances, self.in_alignment,
         self.nonmatching_forward_mismatch_positions, 
         self.nonmatching_reverse_mismatch_positions,
         self.nonmatching_probe_mismatch_positions) = self._get_coverage()

    def _get_alignment_seq_matches(self, primer_object):
        in_alignment = []
        distances = []
        mismatch_positions = []
        for seq in self.alignment:
            expanded_primer_in_current_seq = []
            expanded_primer_distances = []
            expanded_primer_mismatch_pos = []
            for expanded_primer, region in zip(primer_object.expanded_primers,
                    primer_object.expanded_primer_regions):
                expanded_primer_in_current_seq.append(
                    all([segment in seq[region[0]:region[1]]
                        for segment in expanded_primer.split('N')])
                )
                expanded_primer_distances.append(
                    Levenshtein.distance(
                        expanded_primer, seq[region[0]:region[1]]
                        ) - expanded_primer.count('N')  # don't count Ns
                )
                expanded_primer_mismatch_pos.append(
                    [region[0] + i for i in np.arange(len(expanded_primer)) 
                     if expanded_primer[i] != seq[region[0]:region[1]][i]
                     and expanded_primer[i] != 'N'])
            in_alignment.append(any(expanded_primer_in_current_seq))
            distances.append(min(expanded_primer_distances))
            mismatch_positions.append(
                expanded_primer_mismatch_pos[
                    np.argmin(expanded_primer_distances)])
        return in_alignment, distances, mismatch_positions




    def _get_coverage(self):
        if self.probes is None:
            (in_alignment_forward, forward_distances,
                forward_mismatch_positions) = self._get_alignment_seq_matches(
                    self.forward_primers)
            (in_alignment_reverse, reverse_distances,
                reverse_mismatch_positions) = self._get_alignment_seq_matches(
                    self.reverse_primers)
            in_alignment_both = []
            for i in zip(*[in_alignment_forward, in_alignment_reverse]):
                in_alignment_both.append(all(i))
            matching_seqs = np.array(self.alignment)[
                np.array(in_alignment_both)]
            nonmatching_seqs = np.array(self.alignment)[
                ~np.array(in_alignment_both)]
            nonmatching_forward_distances = np.array(forward_distances)[
                ~np.array(in_alignment_forward)]
            nonmatching_reverse_distances = np.array(reverse_distances)[
                ~np.array(in_alignment_reverse)]
            nonmatching_probe_distances = np.empty((1,1))
            nonmatching_forward_mismatch_positions = np.array(
                forward_mismatch_positions)[~np.array(in_alignment_forward)]
            nonmatching_reverse_mismatch_positions = np.array(
                reverse_mismatch_positions)[~np.array(in_alignment_reverse)]
            nonmatching_probe_mismatch_positions = np.empty((1,1))
            coverage = sum(in_alignment_both) / len(in_alignment_both)
        else:
            (in_alignment_forward, forward_distances,
                forward_mismatch_positions) = self._get_alignment_seq_matches(
                    self.forward_primers)
            (in_alignment_reverse, reverse_distances,
                reverse_mismatch_positions) = self._get_alignment_seq_matches(
                    self.reverse_primers)
            (in_alignment_probes, probe_distances,
                probe_mismatch_positions) = self._get_alignment_seq_matches(
                    self.probes)
            in_alignment_both = []
            for i in zip(*[in_alignment_forward, in_alignment_reverse, 
                    in_alignment_probes]):
                in_alignment_both.append(all(i))
            matching_seqs = np.array(self.alignment)[
                np.array(in_alignment_both)]
            nonmatching_seqs = np.array(self.alignment)[
                ~np.array(in_alignment_both)]
            nonmatching_forward_distances = np.array(forward_distances)[
                ~np.array(in_alignment_forward)]
            nonmatching_reverse_distances = np.array(reverse_distances)[
                ~np.array(in_alignment_reverse)]
            nonmatching_probe_distances = np.array(probe_distances)[
                ~np.array(in_alignment_probes)]
            nonmatching_forward_mismatch_positions = np.array(
                forward_mismatch_positions)[~np.array(in_alignment_forward)]
            nonmatching_reverse_mismatch_positions = np.array(
                reverse_mismatch_positions)[~np.array(in_alignment_reverse)]
            nonmatching_probe_mismatch_positions = np.array(
                probe_mismatch_positions)[~np.array(in_alignment_probes)]
            coverage = sum(in_alignment_both) / len(in_alignment_both)
        return (coverage, matching_seqs, nonmatching_seqs, 
                nonmatching_forward_distances, nonmatching_reverse_distances,
                nonmatching_probe_distances, in_alignment_both,
                nonmatching_forward_mismatch_positions,
                nonmatching_reverse_mismatch_positions,
                nonmatching_probe_mismatch_positions)

    def _generate_display_figs(self):
        pass

    def display_report(self):
        """Display report as HTML figures and text in console.
        """
        print(f"Coverage\n{100 * self.coverage:0.1f}%")
        print(f"\nMean nonmatching distances:")
        fig = go.Figure()
        fig.add_trace(go.Histogram(x=self.nonmatching_forward_distances,
            name='Forward Primers'))
        fig.add_trace(go.Histogram(x=self.nonmatching_reverse_distances,
            name='Reverse Primers'))
        fig.update_xaxes(title='Number of Mismatches')
        fig.update_yaxes(title='Number of Sequences')
        if self.probes is not None:
            fig.add_trace(go.Histogram(x=self.nonmatching_probe_distances,
                name='Probe'))
            print(
                f"\nProbes\n"
                f"{np.mean(self.nonmatching_probe_distances):0.1f}")
        print(f"\nForward primers\n"
            f"{np.mean(self.nonmatching_forward_distances):0.1f}")
        print(f"\nReverse primers\n"
            f"{np.mean(self.nonmatching_reverse_distances):0.1f}")
        fig.write_html('number_of_mismatches.HTML')

        fig2 = go.Figure()
        xvals_forward = [item for sublist in 
            self.nonmatching_forward_mismatch_positions for item in sublist]
        xvals_reverse = [item for sublist in
            self.nonmatching_reverse_mismatch_positions for item in sublist]
        fig2.add_trace(go.Histogram(x=xvals_forward,
            name='Forward mismatch positions', xbins=dict(size=1)))
        fig2.add_trace(go.Histogram(x=xvals_reverse,
            name='Reverse mismatch positions', xbins=dict(size=1)))
        fig2.update_xaxes(title='Mismatch Position',
            ticktext=list(
                self.alignment[0][min(xvals_forward):max(xvals_reverse) 
                                        + 1]), 
                tickvals=np.arange(min(xvals_forward), max(xvals_reverse) + 1),
                tickmode='array', tickangle=0, tickfont=dict(size=10))
        fig2.update_yaxes(title='Number of Sequences')
        if self.probes is not None:
            xvals_probe = [item for sublist in 
                self.nonmatching_probe_mismatch_positions
                for item in sublist]
            fig2.add_trace(go.Histogram(x=xvals_probe,
                                        name='Probe mismatch positions', 
                                        xbins=dict(size=1)))
        fig2.write_html('position.HTML')


class MismatchVisualization:
    """ A class for creating data for javascript visualization tool of CoPrimer binding regions and mismatches
    
    Attributes
    -----------
    alignment_path: str
        path to the alignment file given
    forward_primer_names: ls
        names of the forward CoPrimers
    forward_gap_lengths: ls
        list of int values of gap lengths between primer and capture on the CoPrimer
    forward_primer_sequences: ls
        list of the sequences of the linearized forward CoPrimers
    forward_primer_regions: ls 
        list of tuple values for the genome position of the binding regions for the CoPrimers
    reverse_primer_names: ls
        names of the reverse CoPrimers
    reverse_gap_lengths: ls 
        same as forward_gap_lengths, but for the reverse CoPrimers
    reverse_primer_sequences: ls
        list of the sequences of the linearized reverse CoPrimers
    reverse_primer_regions: ls 
        same as forward_primer_regions, but for the reverse CoPrimers
    num_records: int
        the number of sequences analyzed 
    record_len: int
        the length of the sequences analyzed
    pseudo_genome: str
        sequence with the most commonly-occurring, non-ambiguous base at each position
    mismatch_df: DataFrame
        pandas DataFrame with the mismatch values at each position in the genome
    
    Methods
    ----------
    
    get_data: None
        writes 3 CSV files in a data folder that are used in the javascript visualization tool
    amplicon_sequence: str
        returns the sequence of the amplicon needed to cover all the CoPrimers in the class. Can be used to get the positive control sequence """

    def __init__(self, alignment_path, forward_primer_sequences, reverse_primer_sequences,
        forward_primer_names, reverse_primer_names, forward_gap_lengths, reverse_gap_lengths, pseudo_genome=None):

        # determine if Primers are CoPrimers or traditional by checking gap lengths
        if sum(forward_gap_lengths) == 0:
            primer_type = 'traditional'
        elif sum(forward_gap_lengths) > 0:
            primer_type = 'coprimer'
        else:
            raise Exception('Cannot have negative gap lengths')
        self.primer_type = primer_type

    
        #define class variables
        self.alignment_path = alignment_path
        # self.forward_primer_sequences = forward_primer_sequences
        self.forward_primer_names = forward_primer_names
        self.forward_gap_lengths = forward_gap_lengths
        # self.reverse_primer_sequences = reverse_primer_sequences
        self.reverse_primer_names = reverse_primer_names
        self.reverse_gap_lengths = reverse_gap_lengths

        # attributes calculated in class
        self.alignment = parse_alignment(self.alignment_path)
        
        if pseudo_genome != None:
            self.pseudo_genome = pseudo_genome
            print('building sequence matrix...')
            self.seq_array = seq_matrix(self.alignment_path)
        else:
            self.pseudo_genome = self._pseudo_genome()
        

        # get linearized primer sequences
        if primer_type == 'traditional':
            forward_primers = forward_primer_sequences
        else:
            forward_primers = []
            for (primer, gap) in zip(forward_primer_sequences, forward_gap_lengths):
                forward_primers.append(linearize_coprimer(primer, False, gap))
        self.forward_primer_sequences = forward_primers

        # repeat for reverses
        if primer_type == 'traditional':
            reverse_primers = [rc(primer) for primer in reverse_primer_sequences]
        else:
            reverse_primers = []
            for (primer, gap) in zip(reverse_primer_sequences, reverse_gap_lengths):
                reverse_primers.append(linearize_coprimer(primer, True, gap))
        self.reverse_primer_sequences = reverse_primers


        # foward primer regions
        self.forward_primer_regions = self._get_primer_region(self.forward_primer_sequences, self.forward_primer_names, reverse=False)
        
        # reverse primer regions
        self.reverse_primer_regions = self._get_primer_region(self.reverse_primer_sequences, self.forward_primer_names, reverse=True)

        #set up a primers dataframe using pd.DataFrame.from_dict method with a transposition to allow differential array lengths 
        primers_dict = {
            'forward_names': self.forward_primer_names,
            'forward_sequences': self.forward_primer_sequences,
            'forward_gaps': self.forward_gap_lengths, 
            'forward_regions': self.forward_primer_regions,
            'reverse_names': self.reverse_primer_names,
            'reverse_sequences': self.reverse_primer_sequences, 
            'reverse_gaps': self.reverse_gap_lengths, 
            'reverse_regions': self.reverse_primer_regions,
            'reverse_compliment_sequences': [rc(primer) for primer in self.reverse_primer_sequences]
        }

        with open(self.alignment_path) as f:
            self.num_records, self.record_len = self._get_fasta_stats(f)

        self.primers_df = pd.DataFrame.from_dict(primers_dict, orient='index').T

        self.mismatch_df = self._calculate_mismatches()


    def _get_fasta_stats(self, handle):
        """Parses fasta file and returns the number of records and the length of records

        :return num_records: number of alignments in the file
        :return record_len: length of the alignments 
        """
        num_records = 0
        for _, seq in SimpleFastaParser(handle):
            num_records += 1
        record_len = len(seq)

        return num_records, record_len


    def _count_mismatches(self, data):
        """returns the number of bases not equal to the most common base
        """
        ambiguous = ['N', 'M', 'R', 'W', 'S', 'Y', 'K', 'V', 'H', 'D', 'B']
        ser = pd.Series([i for i in data if i not in ambiguous])
        counts = ser.value_counts()
        counts.drop('-', inplace=True)

        common_base = counts.index[0]
        num_seqs = sum(counts)
        try:
            common_mismatch = counts.index[1]
            common_mismatch_freq = counts[1]/num_seqs * 100
        except IndexError:
            common_mismatch = 'None'
            common_mismatch_freq = 'N/A'
        common_freq = counts[0]
        try:
            mismatch_count = sum(counts[1:])
        except IndexError:
            mismatch_count = 0
        return common_base, common_freq, mismatch_count, common_mismatch, common_mismatch_freq, num_seqs
    

    def _common_base(self, data):
        """get the most common base for the position passed
        """
        try:
            series = pd.Series(data).value_counts()
            return str(series.index[0])
        except IndexError:
            return '-'


    def _pseudo_genome(self):
        # convert alignment to numpy array
        with open(self.alignment_path) as f:
            self.num_records, self.record_len = self._get_fasta_stats(f)
            print(f'Calculating mismatches for {self.num_records} records of length {self.record_len}')


        with open(self.alignment_path) as f:
            print('\nanalyzing sequences...')
            name_ls = []
            self.seq_array = np.empty((self.num_records, self.record_len), dtype=str)
            idx = 0
            for name, seq in tqdm(SimpleFastaParser(f), total=self.num_records):
                name_ls.append(name)
                self.seq_array[idx, :] = list(find_insertion_deletion(seq))
                idx += 1

        pseudo_genome = []
        print('\nanalyzing genome..') 
        for i in tqdm(np.arange(self.record_len)):
            pseudo_genome.append(self._common_base([el for el in self.seq_array[:, i] if el in ['A', 'T', 'G', 'C']]))

        pseudo_genome = ''.join(pseudo_genome)

        return pseudo_genome

    
    def _calculate_mismatches(self):

        print('\ncalculating mismatches..')

        # trim the pseudo genome to the range where primers bind +- 20
        start = min(flatten(self.forward_primer_regions)) - 20
        end = max(flatten(self.reverse_primer_regions)) + 20

        # make sure start and end are within the range of the sequence index
        if start < 0:
            start = 0
        if end > self.record_len:
            end = self.record_len
        
        # mismatch_genome = self.pseudo_genome[start: end]
        positions = []
        bases = []
        num_seqs = []
        mismatch_counts = []
        common_mismatches = []
        common_mismatch_freqs = []

        for i in tqdm(np.arange(end-start)):
            pos = start + i + 1
            positions.append(pos) 
            common_base, common_freq, mismatch_count, common_mismatch, common_mismatch_freq, num_seq = self._count_mismatches((self.seq_array[:, pos - 1]))
            bases.append(common_base)
            mismatch_counts.append(mismatch_count)
            common_mismatches.append(common_mismatch)
            common_mismatch_freqs.append(common_mismatch_freq)
            num_seqs.append(num_seq)

        df = pd.DataFrame({
            'positions': positions,
            'bases': bases,
            'mismatches': mismatch_counts, 
            'number_of_sequences': num_seqs, 
            'common_mismatch': common_mismatches,
            'common_mismatch_freq': common_mismatch_freqs, 
            'max_sequences': self.num_records
        })

        return df


    def _get_primer_region(self, primers, primer_names, reverse):
        print('getting primer regions')
        if reverse == True:
            if self.primer_type == 'traditional':
                primers = [rc(primer) for primer in primers]
            else:
                pass
        
        
        regions = []
        for primer, name in zip(primers, primer_names):

            start = -1
            primer_region = None
            try:
                reg = binding_location(primer, self.pseudo_genome, mismatches=0, method='span')
            except:
                try:
                    reg = binding_location(primer, self.pseudo_genome, mismatches=3, method='span') #if there's a mismatch, allow up to 3 mismatches
                except:
                    reg = None
            if reg != None:
                start, end = reg
                primer_region = (start, end)
                regions.append(primer_region)
            else:
                print(f'did not find {name}')
                continue
        return regions


    def get_data(self, folder=None):
        """write the data needed to generate the javascript-based visualization. Will write 3 CSV files to specified folder
        Params:
            folder: The path to the folder where files will be written. If none specified, a folder named "data" will be created in current working directory
        """
        if folder is None:
            folder = './data'
        if not os.path.isdir(folder):
            os.makedirs(folder, exist_ok=True)
        # create list of dictionaries 
        forward_primers = []

        for i in range(len(self.forward_primer_sequences)):
            sequence = self.forward_primer_sequences[i]
            if "N" in sequence:
                primer, capture = re.split('N+', sequence)
            else:
                primer = sequence
                capture = ''
            forward_dict = {
                'name': self.forward_primer_names[i], 
                'sequence': sequence,
                'start_position': self.forward_primer_regions[i][0],
                'stop_position': self.forward_primer_regions[i][1],
                'primer_length': len(primer), 
                'capture_length': len(capture)
            }
            forward_primers.append(forward_dict)


        reverse_primers = []

        for i in range(len(self.reverse_primer_sequences)):
            sequence = self.reverse_primer_sequences[i]
            if "N" in sequence:
                primer, capture = re.split('N+', sequence)
            else:
                primer = sequence
                capture = ''
            reverse_dict = {
                'name': self.reverse_primer_names[i], 
                'sequence': sequence,
                'start_position': self.reverse_primer_regions[i][0],
                'stop_position': self.reverse_primer_regions[i][1],
                'primer_length': len(primer), 
                'capture_length': len(capture)
            }
            reverse_primers.append(reverse_dict)

        # create dataframes 

        mismatches_df = self.mismatch_df
        forward_df = pd.DataFrame(forward_primers)
        reverse_df = pd.DataFrame(reverse_primers)
        
        # write all to a csv
        mismatches_df.to_csv(os.path.join(folder, 'mismatches.csv'), index=False)
        forward_df.to_csv(os.path.join(folder, 'forwards.csv'), index=False)
        reverse_df. to_csv(os.path.join(folder,'reverses.csv'), index=False)

        print('Files written')
    

    def amplicon_sequence(self, overlap_length=1, primers=None):
        """method to get the sequence for the amplicon. will give the sequence from the start of the forward with the lowest position
        to the end of the reverse with the highest position plus the overlap lengths on each side 
        Params: 
            overlap_length -- The number of bases of to return beyond the amplicon region, default=1
            primers (lst, optional) -- the names of the forward and reverse primers for which to return the sequence region. (Default None. returns all primers in class)
        Out:
            str -- sequence for the amplicon of the regions
        """
        if primers is not None:
            if type(primers) != list:
                raise Exception('Primers must be a list')
            
            #build dataframe for forward and reverse primers
            forwards = {
                'name': self.forward_primer_names,
                'position': self.forward_primer_regions
            }

            reverses = {
                'name': self.reverse_primer_names,
                'position': self.reverse_primer_regions
            }
            forward_primer_df = pd.DataFrame(forwards)
            reverse_primer_df = pd.DataFrame(reverses)

            #keep only primers specified in the method call
            selected_forwards = forward_primer_df[forward_primer_df['name'].isin(primers)]
            selected_reverses = reverse_primer_df[reverse_primer_df['name'].isin(primers)]

            forward_regions = selected_forwards['position']
            reverse_regions = selected_reverses['position']
        

        else:
            forward_regions = self.forward_primer_regions
            reverse_regions = self.reverse_primer_regions 

        start_position = min(forward_regions)[0] - overlap_length #(overlap_length) bases before the start of the min forward primer
        end_position = max(reverse_regions)[1] + overlap_length 

        #make sure the positions are in index range
        if start_position < 0:
            print('The start position given was less than 0, the sequence given begins at the start')
            start_position = 0
        if end_position >= len(self.pseudo_genome):
            end_position = len(self.pseudo_genome) - 1 #since 0 index has to be 1 less than the length of the sequence
        
        return self.pseudo_genome[start_position: end_position]
        

    def save_class(self, filename=None):
        """ Save class instance as binary pickle object. The class can be reloaded with pickle.load()
        Parameters:
            filename (str/PathLike, Optional) -- file name to save the class to (default None)
        Return:
            None - class will be written as binary pickle object """
        if filename is None:
            filename = 'mismatch_viz_class.pkl'
    
        with open(filename, 'wb') as file:
            pickle.dump(self, file, -1)
        print(f'class saved as {filename}')


    def check_amplicon(self, sequence):
        """ 
        Check if the primers in the class will bind to the specified sequence and return a list of the primers that will bind
        Arguments:
            sequence (str) -- sequence to check if primers bind to it
        Return:
            list -- returns list of all primer names that bind to the sequence    
        """
        df = self.primers_df
        matching_forwards = []
        matching_reverses = []
        non_binding_forwards = []
        non_binding_reverses = []
    
        for i, row in df.iterrows():
            # check if forward matches
            pattern = re.compile(row['forward_sequences'].replace('N', '.'))
            match = pattern.search(sequence)
            if match:
                matching_forwards.append(row['forward_names'])
            else:
                non_binding_forwards.append(row['forward_names'])
            
            # check if reverse matches
            pattern = re.compile(row['reverse_compliment_sequences'].replace('N', '.'))
            match = pattern.search(sequence)
            if match:
                matching_reverses.append(row['reverse_names'])
            else:
                non_binding_reverses.append(row['reverse_names'])
        
        if len(non_binding_forwards) == 0 and len(non_binding_reverses) == 0:
            print('All primers bind to the sequence')
        else:
            print(f'non-binding forwards: {non_binding_forwards} \nnon-binding_reverses: {non_binding_reverses}')
        
        total_binding = matching_forwards + matching_reverses

        if len(total_binding) == 0:
            print('None of the primers bind to that sequence')
        
        return total_binding


## Needs work. Create input file that can be parsed. Make a separate function.
# if __name__ == "__main__":
#     parser = argparse.ArgumentParser(
#         description=('Calculate mismatch histograms for primers against a '
#             'given alignment'))
#     parser.add_argument('-i', '--input_alignment',
#         help='path to input .aln alignment file', required=True)
#     parser.add_argument('-p', '--primer_file',
#         help=('path to tab-delimited primer file with format '
#             '<forward|reverse|probe><\t><coprimer|primer><\t>'
#             '<gap length if a coprimer '
#             'or probe orientation if a probe (forward or reverse_complement)>'
#             '<\t><sequence><\n>'))
#     parser.add_argument('-o', '--output',
#         help='path to output html file', required=True)
#     args = parser.parse_args()
#     primer_info = pd.read_csv(args.primer_file, header=None, sep='\t')
#     coprimers = primer_info[primer_info[1] == 'coprimer']
#     primers = primer_info[primer_info[1] == 'primer']
#     for primer_df in [coprimers, primers]:
#         if not primer_df.empty:
#             forward_primers = primer_df.loc[
#                 primer_df[0] == 'forward', 3].tolist()

#             reverse_primers = primer_df.loc[
#                 primer_df[0] == 'reverse', 3].tolist()
#             probes = primer_df.loc[
#                 primer_df[0] == 'probe', 3].tolist()
#     coverage = PrimerCoverage()



# test out the class

