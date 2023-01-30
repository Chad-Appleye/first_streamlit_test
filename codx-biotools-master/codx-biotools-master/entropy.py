from copy import Error
from os.path import basename
from unicodedata import normalize
import pandas as pd
from pandas.core.indexes.base import Index
import scipy.stats
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from tqdm import tqdm  # progress bar   
import pickle
import os
import random

from pathlib import Path
from collections import Counter
from Bio.SeqIO.FastaIO import SimpleFastaParser
from oligotools import calculate_gc as gc, calculate_tm as tm, dotplot, find_missing_x, find_missing_y, largest_gap, reverse_complement, trim_to_tm, binding_location
from math import log10, floor


def find_insertion_deletion(sequence):
    """ Using Regex, find the '-' caracters in the sequence that correspond to insertion or deletion by identifying if 
    the '-' has sequenced bases before and after it
    
    Params:
    input: sequence (str)
    out: sequence with 'J' in place of '-' which are insertions or deletions (str)
    """
    lstripped = sequence.lstrip('-')
    start = len(sequence) - len(lstripped)
    stripped = lstripped.rstrip('-')
    end = len(lstripped) - len(stripped)
    return '-' * start + stripped.replace('-', 'J') + '-' * end


def mismatch_count(align1, align2, score, begin, end):
    gap = 0
    mismatch = 0
    for a, b in zip(align1[begin:end], align2[begin:end]):
        if a == b:
            continue
        elif a == '-' or b == '-':
            gap += 1
        else:
            mismatch += 1
    return gap, mismatch


def write_target_dict(target_dict, file_name): 
    """ Write the target dictionary to a text file in an easy-to-read format. 
    Parameters:
        target_dict (dict) -- Dictionary of the targets from frind_targets() method of Entropy class
        file_name (str/Path) -- file name of the text file to be written
        
    Returns:
        None- text file will be written
    """
    #make the directories to the file, if necessary
    path = Path(file_name)
    parent = path.parent
    os.makedirs(parent, exist_ok=True)

    # iterate through directory and write it
    with open(file_name, 'w+') as f:
        for key, val in target_dict.items():
            if type(val) == dict:
                f.write(str(key) + '\n')
                for key2, val2 in val.items():
                    f.write('\t' + str(key2) + ': ' + str(val2) + '\n')
            else:
                f.write(str(key) + ': ' + str(val) + '\n')
    
    print('target file written at: {}'.format(file_name))


def depth_counts(seq_array, normalize=True):
    num_seqs = seq_array.shape[0]
    seq_array_bool = seq_array != '-' # convert the array to a boolean array 1 if not = '-', 0 otherwise
    counts = np.sum(seq_array_bool, axis=0) # get 1-D count of the sum for each position
    if normalize is True:
        return counts / num_seqs # normalize the counts
    else:
        return counts


def graph_depth(seq_array, normalize=True):
    """ 
    Graph the quality of sequences (proportion of sequences that were sequenced successfully) along the span of the sequence. 

    Arguments:
        seq_array (ndarray) -- a numpy array of aligned sequences
    Returns:
        Plotly figure 

    """
    num_seqs = seq_array.shape[0]
    counts = depth_counts(seq_array, normalize=normalize)

    fig = px.line(counts)
    fig.update_layout(title=f'sequence depth chart. Sequences: {num_seqs}', xaxis_title='Genome Position', yaxis_title='Seqeunced Proportion', showlegend=False)
    
    return fig


def all_mismatches(target_dict, entropy_class):
    """ 
    Describe the mismatches for every sequence in a dictionary. Generally used to evaluate targets from Entropy.find_targets()

    Arguments:
        target_dict (dict) -- dictionary containing target sequences to be evaluated
        entropy_class (obj) -- instantiated entropy class for the targets being passed
    Return: 
        (dict) target names and mismatch count Series -- mismatch calculations will be printed to the console
    """

    ent = entropy_class

    if type(target_dict) != dict:
        raise TypeError('target_dict must be a dictionary')
    
    out_dict = {}

    # iteratively evaluate each sequence
    for key, val in target_dict.items():
        if type(val) == dict:
            for key2, val2 in val.items():
                if 'seq' in key2.lower():
                    mismatch = ent.target_mismatch_counts(val)
                    print(f'{key}:\n{mismatch}')
                    out_dict.update({key2: mismatch})
                else:
                    continue
        elif 'seq' in key.lower():
            mismatch = ent.target_mismatch_counts(target_dict)
            print(f'{target_dict}:\n{mismatch}')
            out_dict.update({key: mismatch})
        else:
            continue
    return out_dict



class Entropy: 

    def __init__(self, fasta_path, target_name='Target'):
        """Instantiating an Entropy class will automatically calculated the entropy for 
        each position on the provided alignment files. 

        :param target_name: name of organism and/or gene of interest
        :param fasta_path: path to the .fasta.aln file where the alignments are stored

        :return: None, class will be instantiated and output can be specified using methods
        """
        self.target_name = target_name
        self.fasta_path = fasta_path
        self.df, self.pseudo_sequence, self.mismatches = self._generate_df() #dataframe of entropy values
        

    def _get_fasta_stats(self, handle):
        """Parses fasta file and returns the number of records and the length of records

        :return num_records: number of alignments in the file
        :return record_len: length of the alignments 
        """
        num_records = 0
        sequences = []
        for _, seq in SimpleFastaParser(handle):
            num_records += 1
            sequences.append(seq)
        record_len = len(sequences[0])

        sequences = [find_insertion_deletion(i) for i in sequences]

        return num_records, record_len, sequences


    def _entropy(self, data):
        """Returns an entropy calculation from input data
        :param data: Pandas series

        :return: entropy value
        """
        data = pd.Series(data)
        counts = data.value_counts()
        return scipy.stats.entropy(counts)
    
    
    def show(self, window=30, range=None, title=None):
        """show graph of entropy values
        Arguments:
            window (int, Optional) -- The number of positions to include in the rolling average (Default 30)
            range (lst, optional) -- The x-axis range if the entire genome is not desired (Default None)
            title (str, optional) -- If title different from target_name in class instance, pass title to update graph title (Default None)
        Return:
            None. Graph will be shown in console
        """
        fig = self._generate_fig(window, range, title)
        fig.show()
    

    def write_file(self, file_name=None, window=30, range=None, title=None):
        """Write html file with the graph of entropy values
        Arguments:
            file_name (str, Optional) -- file name to write html graph (Default ./{target_name}_entropy.html)
            window (int, Optional) -- The number of positions to include in the rolling average (Default 30)
            range (lst, optional) -- The x-axis range if the entire genome is not desired (Default None)
            title (str, optional) -- If title different from target_name in class instance, pass title to update graph title (Default None)
        Return:
            None. HTML graph file will be written 
        """
        fig = self._generate_fig(window, range, title)
        if file_name == None:
            fig.write_html(f'{self.target_name}_entropy.html')
        else: 
            if not file_name.endswith('.html'):
                file_name = file_name.split('.')[0] + '.html'
            parent_dir = Path(file_name).parent
            os.makedirs(parent_dir, exist_ok=True)
            basename = os.path.basename(file_name)
            filename = os.path.join(parent_dir, basename)
            try:
                fig.write_html(filename)
            except:
                print("That is not a valid file name")


    def _generate_df(self):
        """Creates the dataframe of the entropy values and a pseudo sequence  of most common bases
        Output:
            df: dataframe of the entropy values
            pseudo-sequence: sequence of the most common bases at each position of the alignment
            mismatches: list of the number of mismatches at each position
        """
        # convert alignment to numpy array
        with open(self.fasta_path) as f:
            
            self.num_records, self.record_len, all_sequences = self._get_fasta_stats(f)
            print(f'\n\nlooking at {self.num_records} sequences of length {self.record_len}...')

        with open(self.fasta_path) as f:
            name_ls = []
            seq_array = np.empty((self.num_records, self.record_len), dtype=str)
            idx = 0
            print('\nProcessing alignments...\nThis may take a while')
            for name, seq in tqdm(SimpleFastaParser(f), total=self.num_records, unit='seq'):  
                name_ls.append(name)
                seq_array[idx, :] = list(find_insertion_deletion(seq)) #inserted function here to change '-' to 'J' when it represents an insertion or deletion
                idx += 1
            self.seq_array = seq_array #save as class attribute so it can be accessed later
        
        #calculate entropy for each position
        positions = []
        entropy_values = []
        pseudo_sequence = []
        mismatches = []
        allowable_bases = ['A', 'T', 'G', 'C', 'J']
        print('\ncomputing entropy...')
        for i in tqdm(np.arange(self.record_len), unit='pos', desc='Genome'):
            positions.append(i + 1)
            entropy_values.append(self._entropy([el for el in seq_array[:, i] if el in allowable_bases]))
            pseudo_sequence.append(self._common_base([el for el in seq_array[:, i] if el in allowable_bases]))
            mismatches.append(self._count_mismatches([el for el in seq_array[:, i]]))
        #create dataframe of entropy values and rolling average values
        df = pd.DataFrame({'position': positions, 'entropy': entropy_values})
        return df, ''.join(pseudo_sequence), mismatches


    def _generate_fig(self, window, range=None, title=None):
        """ generate the plotly figure object with genome location on the x axis and entropy values on the y axis
        Parameters:
            window (int) -- The number of positions to include in the rolling average
            range (lst, optional) -- The x-axis range if the entire genome is not desired (Default None)
            title (str, optional) -- If title different from target_name in class instance, pass title to update graph title (Default None)
        Return:
            plotly graph object
        """
        df = self.df
        df['window'] = df.entropy.rolling(window).mean()
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=df.position, y=df.entropy, mode='lines', name='Entropy Values'))
        fig.add_trace(go.Scatter(x=df.position, y=df.window, line=dict(color='red', width=3), name=f"Rolling Average = {window}", showlegend=True))
        fig.update_layout(title=self.target_name,
                        xaxis_title='genome position',
                        yaxis_title='Entropy value')
        if range is not None:
            if type(range) != list:
                raise Exception('range type must be list with 2 values')
            fig.update_xaxes(range=range)
        if title is not None:
            fig.update_layout(title=title)
        return fig


    def graph_entropy_and_depth(self, window=30, range=None, title=None, depth_normalize=False):
        """ Graph the entropy and depth in a single graph 
        Arguments:
            window (int) -- rolling average window value (Default 30)
            range (list) -- x-axis range of the graph (Default None)
            title (str) -- title of the graph, if None, a title will be generated (Default None)
            depth_normalize (bool) -- if True, the depth values will be normalized between 0 and 1
        Return 
            Plotly figure object
        """

        try:
            int(window)
        except:
            raise TypeError('window must be an integer')

        df = self.df
        df['window'] = df.entropy.rolling(window).mean()
        df['depth'] = depth_counts(self.seq_array, normalize=depth_normalize)
        
        fig = go.Figure()

        fig.add_trace(go.Scatter(x=df['position'], y=df['depth'], name='depth', yaxis='y1'))
        fig.add_trace(go.Scatter(x=df['position'], y=df['window'], name=f'entropy rolling average = {window}', yaxis='y2', opacity=.8))
        fig.add_trace(go.Scatter(x=df['position'], y=df['entropy'], name='entropy values', yaxis='y2', opacity=.7))
    

        fig.update_layout(
        yaxis=dict(
            title='depth values', 
            side='left'
        ), 
        yaxis2=dict(
            title='entropy values', 
            side='right',
            overlaying='y'
        ))
        if range is not None:
            range_error = 'range must be a list of 2 values'
            if type(range) != list or len(range) != 2:
                raise Exception(range_error)
            fig.update_xaxes(range=range)
            
        if title is not None:
            fig.update_layout(title=title)
        else:
            title = f'Entropy and depth values for {self.num_records} sequences of {self.target_name}'
            fig.update_layout(title=title)

        return fig

    
    def _common_base(self, data):
        """get the most common base for the position passed
        """
        try:
            return Counter(data).most_common(1)[0][0]
        except IndexError:
            return '-'

    
    def _count_mismatches(self, data):
        """returns the number of bases not equal to the most common base
        """
        most_common = Counter(data).most_common(1)[0][0]
        return sum(1 for i in data if i not in [most_common, 'N', 'W', '-', 'Y', 'R'])


    def get_sequence(self, start, end):
        """returns a pseudo sequence of the input alignment where each base is the most common 
        base at that position in the alignment

        Args:
            start: The starting position in the alignment for the sequence
            end: the ending position in the alignment for the sequence
        
        Output: 
            pseudo-sequence: sequence at the specified positions of the most common bases in the alignment
        """
        return self.pseudo_sequence[start:end+1]

    
    def mismatch(self, start, end):
        """analyze the mismatches 
        Args:
            start: The starting position in the alignment for the sequence
            end: the ending position in the alignment for the sequence
        Output:
            plotly histogram of the number of mismatches at the specified location
        """
        data = self.mismatches[start:end+1]
        position = np.arange(start, end+1)
        sequence = list(self.pseudo_sequence[start:end+1])
        df = pd.DataFrame(data=data, columns=['mismatches'])
        df['position'] = position
        df['base'] = sequence
        lower_range = min(df['mismatches']) - min(df['mismatches']) * .2
        upper_range = max(df['mismatches']) + max(df['mismatches']) * .2
        title = f'{self.target_name} {start} to {end} for {self.num_records} sequences'
        fig = px.bar(df, x='position', y='mismatches')
        fig.update_layout(title=title, 
            xaxis_title='Sequence', 
            yaxis_title='Number of Mismatches',
            xaxis=dict(
                nticks=len(data),
                tickvals = df['position'],
                ticktext = df['base']
            ),
            yaxis=dict(range=[lower_range, upper_range]

            ))
        fig.write_html(f'{self.target_name}_mismatch.html')


    def save_class(self, filename=None):
        if filename == None:
            filename = f'{self.target_name}_entropy_class.pkl'
        else: 
            filename = filename
        
        with open(filename, 'wb') as output:
            pickle.dump(self, output, -1)
        print(f'class saved as {filename}')  


    def find_targets(self, target_length=50, gc_cutoff=40, percent=25, max_amplicon=120, max_entropy=None, start=None, stop=None, method='max'):
        """
        Method to find target regions with lowest peak entropy values for CoPrimer design 

        Params:
            target_length (int) -- the length of the sequence to be used in design of the CoPrimers (Default 50)
            gc_cutoff (int) -- the minimum acceptable GC percentage of the target sequence (Default 40)
            percent (int) -- the percentage of the data to use. Use a higher percentage if no results are given (Default 25)
            max_amplicon (int) -- the maximum allowed amplicon length from start of forward position to end of reverse position (Default 120)
            max_entropy (int) If specified, only targets with a maximum entropy value below the specified value will be returned. Values range from 0 to 1 (Default None)
            start (int) -- Starting position allowed for target -useful if designing on specific gene target (Default None)
            stop (int) -- Ending position allowed for target -useful if designing on specific gene target (Default None)
            method (str) -- method to quantify entropy in target window - one of ["max", "median", "mean"] - Default max

        Return: 
            Dictionary containing the target regions, sequences, and GC%

        """
        # Initial error handling
        # make sure the target_length * 2 isn't larger than the max_amplicon argument as this will result in no targets being found
        if target_length * 2 > max_amplicon:
            raise Exception('target_length of {} and max_amplicon of {} not allowed. The target_length X 2 must be <= max_amplicon.'.format(target_length, max_amplicon))

        if percent <= 0 or percent > 100:
            raise Exception('percent = {} not allowed, percent must be a value between 0 and 100'.format(percent))

        if target_length < 50:
            print('A target length < 50 is not likely to return a valid CoPrimer design.')

        methods = ["max", "median", "mean"]
        method_cleaned = str(method).strip().lower()
        if method_cleaned not in methods:
            raise Exception(f'method: {method} not recognized')
        
        else:
            method = method_cleaned
        
        new_col = f'rolling_{method}'

        sequence = self.pseudo_sequence 

        df = self.df.copy()
        
        calculations = {
            'mean': df.entropy.rolling(window=target_length).mean(),
            'max': df.entropy.rolling(window=target_length).max(), 
            'median': df.entropy.rolling(window=target_length).median()
        }

        for key in calculations.keys():
            df[f'rolling_{key}'] = calculations[key]

        df = df.dropna()
        df = df.sort_values(by=new_col, ascending=True)

        # filter dataframe based on start and stop arguments
        if start is not None:
            if type(start) != int:
                raise Exception('start postion must be integer')
            df = df[df['position']>=(start + target_length)]
        
        if stop is not None: 
            if type(stop) != int:
                raise Exception('end position must be integer')
            df = df[df['position']<=stop]

        # use only the lowest specified entropy values
        rows_to_keep = round(len(df) * percent/100)

        df = df.iloc[:rows_to_keep]

        # check max_entropy, if specified
        if max_entropy is not None: 
            df = df[df[new_col]<=max_entropy]

        # find the first target region that meets criteria
        print('\nAnalyzing possible targets..')
        print('May finish before progress bar is complete')
        for i in tqdm(range(len(df))):
            row_1 = df.iloc[i]
            end_position_1 = int(row_1.position) #position is 1 ahead of index, so this will be inclusive 
            start_position_1 = end_position_1 - (target_length - 1) #the number at the index is included in the max calculation, so only the previous (target_length -1) bases will also be included
            target_sequence_1 = sequence[start_position_1 -1: end_position_1]
            gc_percent_1 = gc(target_sequence_1)

            # it has to meet the GC percent cutoff specified
            if gc_percent_1 < gc_cutoff:
                continue

            #don't allow '-' or 'J' to be in the target sequence - this will filter out frequently unread regions at the beginning or end of reads and common insertions and deletions
            if '-' in target_sequence_1 or 'J' in target_sequence_1:
                continue

            try: 
                target_1 = {
                    'start_position': start_position_1,
                    'end_position': end_position_1, 
                    'sequence': target_sequence_1, 
                    'GC_percent': gc_percent_1, 
                    'maximum_entropy': round(row_1.rolling_max, 2 - int(floor(log10(abs(row_1.rolling_max))))-1) #round to 2 sig figs
                }
            except ValueError:
                target_1 = {
                    'start_position': start_position_1,
                    'end_position': end_position_1, 
                    'sequence': target_sequence_1, 
                    'GC_percent': gc_percent_1, 
                    'maximum_entropy': round(row_1.rolling_max, 1)
            }
            # find a second target within the specified range from the amplicon_length
            for j in range(i + 1, len(df)):
                row_2 = df.iloc[j]
                end_position_2 = int(row_2.position)
                if abs(end_position_2 - end_position_1) < target_length: #make sure the second target is at least 'target_length' away from the first target
                    continue 
                start_position_2 = end_position_2 - (target_length - 1)
                

                # calculate the amplicon length
                if start_position_1 < start_position_2:
                    amplicon_length = end_position_2 - start_position_1
                else: 
                    amplicon_length = end_position_1 - start_position_2
                
                # skip if the amplicon length exceeds the specified max
                if amplicon_length > max_amplicon:
                    continue

                # make sure the G/C content meets criteria
                target_sequence_2 = sequence[start_position_2 - 1: end_position_2]
                gc_percent_2 = gc(target_sequence_2)
                if gc_percent_2 < gc_cutoff:
                    continue

                # don't allow '-' or 'J' in the target sequence to filter out common insertions and deletions as well as commonly unread regions before and after the sequence
                if '-' in target_sequence_2 or 'J' in target_sequence_2:
                    continue

                try: 
                    target_2 = {
                        'start_position': start_position_2,
                        'end_position': end_position_2, 
                        'sequence': target_sequence_2, 
                        'GC_percent': gc_percent_2, 
                        'maximum_entropy': round(row_2.rolling_max, 2 - int(floor(log10(abs(row_2.rolling_max))))-1) #round to 2 sig figs
                    }
                except ValueError:
                    target_2 = {
                        'start_position': start_position_2,
                        'end_position': end_position_2, 
                        'sequence': target_sequence_2, 
                        'GC_percent': gc_percent_2, 
                        'maximum_entropy': round(row_2.rolling_max, 1)
                }
                
                if start_position_1 < start_position_2:
                    targets_dict = {
                        'forward': target_1,
                        'reverse': target_2
                    }
                else: 
                    targets_dict = {
                        'forward': target_2,
                        'reverse': target_1
                    }

                # add the largest possible amplicon length to the dictionary
                targets_dict['max_amplicon'] = targets_dict['reverse']['end_position'] - targets_dict['forward']['start_position']

                print('\nFinished')
                return targets_dict

        # if this is reached, it means 2 targets couldn't be found with the given criteria
        raise Exception('Could not find targets matching those criteria. Try widening search criteria')


    def target_mismatch_counts(self, target_sequence):
        """ 
        Counts the number of mismatches that each sequence in the class alignment would have with the given target and returns a pandas Series count. Mismatches and gaps are both treated as an error
        Params:
            target_sequence: The target sequence to analyze, generally a target given from the find_targets() method. Can be str or dict
        Output:
            Pandas Series value_counts of the percentage of sequences with that number of 'errors' for the given target region 
        
        """

        regular_bases = ['A', 'T', 'G', 'C']

        #use all sequences in the alignment
        sequences = self.seq_array
        genome = self.pseudo_sequence

        if type(target_sequence) == dict:
            if 'sequence' in target_sequence.keys():
                target_sequence = target_sequence['sequence']
            else: 
                raise Exception("'sequence' is not one of the dictionary keys")

        elif type(target_sequence) != str:
            raise Exception('Unrecognized data type for target_sequence') 

        target_sequence = target_sequence.strip().upper()


        range, is_reverse = binding_location(target_sequence, genome, mismatches=3, method='span', return_if_reverse=True)

        if is_reverse == True:
            target_sequence = reverse_complement(target_sequence)

        errors = []

        for seq in tqdm(sequences):
            seq = ''.join(seq[range[0]: range[-1]])
            if seq == '-' * len(seq):
                continue
            if seq == target_sequence: #if the sequence is found, there are no errors
                errors.append(0)
                continue
            else:
                num_errors = 0
                for a, b in zip(target_sequence, seq):
                    if a not in regular_bases:
                        continue
                    elif a != b:
                        if a == 'I':
                            continue
                        if b in (regular_bases + ['J']):
                            num_errors +=1 #treat gaps and mismatches both as errors
                        else:
                            continue
                errors.append(num_errors)

        
        errors_series = pd.Series(errors)

        return errors_series.value_counts(normalize=True, dropna=False).sort_index() * 100


    def all_mismatch_counts(self, target_dict):
        """ 
        Describe the mismatches for every sequence in a dictionary. Generally used to evaluate targets from Entropy.find_targets()

        Arguments:
            target_dict (dict) -- dictionary containing target sequences to be evaluated
        Return: 
            None -- mismatch calculations will be printed to the console
        """

        all_mismatches(target_dict, self)


    def write_targets(self, file_name, target_length=50, gc_cutoff=40, percent=25, max_amplicon=120, max_entropy=None, start=None, stop=None):
        """ Method to find targets and write them to a text file. Combines find_targets method and write_target_dict function
        Arguments: 
            file_name (str/Path) -- location where the target will be written
            target_length (int) -- the length of the sequence to be used in design of the CoPrimers (Default 50)
            gc_cutoff (int) -- the minimum acceptable GC percentage of the target sequence (Default 40)
            percent (int) -- the percentage of the data to use. Use a higher percentage if no results are given (Default 25)
            max_amplicon (int) -- the maximum allowed amplicon length from start of forward position to end of reverse position (Default 120)
            max_entropy (int) If specified, only targets with a maximum entropy value below the specified value will be returned. Values range from 0 to 1 (Default None)
            start (int) -- Starting position allowed for target -useful if designing on specific gene target (Default None)
            stop (int) -- Ending position allowed for target -useful if designing on specific gene target (Default None)

        Return: 
            targets dict. text files will be written
         """

        if not file_name.endswith('.txt'):
            file_name = file_name.split('.')[0] + '.txt'

        targets = self.find_targets(target_length=target_length, gc_cutoff=gc_cutoff, percent=percent, max_amplicon=max_amplicon, max_entropy=max_entropy, start=start, stop=stop)
        if len(targets) == 0:
            raise Exception('No targets to write, try widening parameters')
        write_target_dict(target_dict=targets, file_name=file_name)

        return targets


    def depth(self, normalize=False):
        """ 
        Graph the quality of sequences (proportion of sequences that were sequenced successfully) along the span of the sequence. 

        Arguments:
            seq_array (ndarray) -- a numpy array of aligned sequences
            normalize (bool) -- normalize the data as proportion (Default False)
        Returns:
            Plotly figure 
        """  

        return graph_depth(self.seq_array, normalize=normalize)


    def primer_design(self, primer_tm, probe_tm, amplicon_length, percent=25, max_probe_length=45, start=None, stop=None, quiet=False):
        """ 
        design regular primers with taqman-style probe on low entropy region. 
        Arguments:
            primer_tm (int) -- predicted melting temperature of the primers
            probe_tm (int) -- predicted melting temperature of the probe 
            amplicon_length (int) -- the length in bp of the amplicon for the design
            percent (int) -- the percentage of the data to use. Use a higher percentage if no results are given (Default 25)
            max_probe_length (int) -- maximum allowable probe length (Default 45)
            start (int) -- starting position of potential designs
            stop (int) -- ending position of potential designs
            quiet (bool) -- if calling method iteratively, set quiet to True to silence the output on each iteration
        Return: 
            dictionary with the specifications of the primers and probe
        """

        try:
            primer_tm = int(primer_tm)
            probe_tm = int(probe_tm)
            max_probe_length = int(max_probe_length)
            amplicon_length = int(amplicon_length -1)
            if start != None:
                start = int(start)
            if stop != None:
                stop = int(stop)
        except:
            raise 
        
        df = self.df.copy()

        df['rolling_max'] = df['entropy'].rolling(window=amplicon_length).max()
        df.dropna(inplace=True)
        df.sort_values(by='rolling_max', ascending=True, inplace=True)

        # filter dataframe based on start and stop arguments
        if start is not None:

            df = df[df['position']>=(start + amplicon_length)]
        
        if stop is not None: 
            df = df[df['position']<=stop]

        # use only the lowest specified entropy values
        rows_to_keep = round(len(df) * percent/100)

        df = df.iloc[:rows_to_keep]

        # try design iteratively 
        for i, row in df.iterrows():
            end_position = int(row.position)
            start_position = int(end_position - amplicon_length)
            sequence = self.pseudo_sequence[start_position-1: end_position] # since position is ahead of index by 1, this will be inclusive of the end position
            forward_primer = trim_to_tm(sequence, primer_tm, method='right', return_tm=True, quiet=quiet)
            reverse_primer = trim_to_tm(sequence, primer_tm, method='left', return_tm=True, quiet=quiet)
            forward_seq, forward_tm = forward_primer.values()
            reverse_seq, reverse_tm = reverse_primer.values()

            # define the probe region to be the middle section of the amplicon separated by 2 bp from the forward and reverse to account for exonuclease activity
            probe_region = sequence[(len(forward_seq) + 2): -(len(reverse_seq)+2)]
            
            if tm(probe_region) < probe_tm: # if the probe TM is already too small, this can't be a viable design. Likely, the sequence is too short or the GC content is too small
                continue
            else:
                probe = trim_to_tm(probe_region, probe_tm, method='right', return_tm=True, quiet=quiet)
                probe_seq, probe_tm_actual = probe.values()
            if len(probe_seq) > max_probe_length: # lower the allowable TM by 0.5 to solve too long probe sequence problem if it is very close
                probe_seq = trim_to_tm(probe_seq, probe_tm - 0.5, method='right') 
                if len(probe_seq) > max_probe_length:
                    continue

            # prepare dictionaries of the design    
            forward_primer.update({
                'start_position': start_position,
                'end_position': start_position + len(forward_seq),
                'GC_content': gc(forward_seq), 
                'length': len(forward_seq)
            })
        
            reverse_primer.update({
                'reverse_complement': reverse_complement(reverse_primer['sequence']),
                'start_position': end_position - len(reverse_seq),
                'end_position': end_position,
                'GC_content': gc(reverse_seq), 
                'length': len(reverse_seq)
            })

            probe.update({
                'sequence': probe_seq,
                'tm': tm(probe_seq),
                'start_position': start_position + len(forward_seq) + 2,
                'end_position': start_position + len(forward_seq) +2 + len(probe_seq),
                'GC_content': gc(probe_seq),
                'length': len(probe_seq)
            })

            bases = ['A', 'T', 'G', 'C']

            actual_amplicon_seq = self.pseudo_sequence[start_position -1 : end_position] # this is the actual amplicon, not added adjacent bases or extended for gblock synthesis
            # amplicon_seq here describes the gblock sequence, which must be at least 125 for synthesis
            try:
                
                amplicon_seq = self.pseudo_sequence[start_position-2: end_position+1] # adding on adjacent bases
            
            except IndexError:
                try:
                    amplicon_seq = 'C' + self.pseudo_sequence[start_position-1: end_position+1]
                except IndexError:
                    amplicon_seq = self.pseudo_sequence[start_position-1: end_position] + 'C'
            
            # now extend the gblock sequence until it reaches at least 125 length
            extend = 1
            while len(amplicon_seq) < 125:
                try:
                    amplicon_seq = self.pseudo_sequence[start_position - 2 - extend: end_position + 1 + extend]
                except IndexError:
                    random_bases = []
                    for i in range(extend):
                        random_bases.append(random.choice(bases))
                    amplicon_seq = amplicon_seq + ''.join(random_bases)
                extend += 1
            
            amplicon = {
                'sequence': actual_amplicon_seq,
                'gblock_seq': amplicon_seq,
                'actual_length' : len(actual_amplicon_seq),
                'gblock_length': len(amplicon_seq)
            }
            design = {
                'forward': forward_primer,
                'probe': probe,
                'reverse': reverse_primer,      
                'amplicon': amplicon
            }
            return design
        if quiet == True:
            return None
        else:
            print('No design could be found with the given parameters')
            return None 
    

    def minimal_amplicon_design(self, primer_tm, probe_tm, min_amplicon,  max_amplicon, percent=25, max_probe_length=45, start=None, stop=None):
        """ 
        Iteratively calls primer_design until a successful design is returned with the smallest possible amplicon that meets the parameters
        Arguments:
            primer_tm (int) -- predicted melting temperature of the primers
            probe_tm (int) -- predicted melting temperature of the probe 
            min_amplicon (int) -- minimum allowable amplicon size
            max_amplicon (int) -- maximum allowable amplicon size
            percent (int) -- the percentage of the data to use. Use a higher percentage if no results are given (Default 25)
            max_probe_length (int) -- maximum allowable probe length (Default 45)
            start (int) -- starting position of potential designs
            stop (int) -- ending position of potential designs
        Return: 
            dictionary with the specifications of the primers and probe
        """
        try:
            int(min_amplicon)
            abs(min_amplicon)
        except:
            raise Exception('min_amplicon must be an integer value')

        try:
            int(max_amplicon)
            abs(max_amplicon)
        except ValueError:
            raise Exception('max_amplicon must be an integer value')

        # pick a starting amplicon length to begin with. 
        amplicon_length = min_amplicon

        for i in range(amplicon_length, max_amplicon + 1):
            design = self.primer_design(primer_tm, probe_tm, i, percent, max_probe_length, start, stop, quiet=True)
            if design is None:
                continue
            else:
                return design
        
        print('No primers could be designed given those parameters')
        return None
        

    def find_copter_target(self, target_length=150, gc_cutoff=45, start=None, stop=None, method='mean'):
        """ 
        Find a target for CopTer design. 
        Will search the entire genome or given region and return the window with the lowest entropy calculated by the given method

        Arguments: 
            target_length (int) -- length of the target to return for CopTer design -Default 150
            gc_cutoff (int) -- lowest acceptable GC percent for the region - Default 45
            start (int) -- position of the genome to start the search - Default None will search entire genome
            stop (int) -- position of the genome to end the search - Default None will search entire genome
            meethod (str) -- quantification method for entropy, one of ["max", "median", "mean"] - Default mean
        
        Return: Dict with the target sequence, region, GC percent, and position
        """
        window = target_length
        method = method.lower().strip()
        accepted_methods = ['mean', 'median', 'max']

        if method not in accepted_methods:
            raise Exception(f'method {method} not possible, must be one of {accepted_methods}')

        sequence = self.pseudo_sequence 

        df = self.df.copy()
        calculations = {
            'mean': df.entropy.rolling(window=window).mean(),
            'max': df.entropy.rolling(window=window).max(), 
            'median': df.entropy.rolling(window=window).median()
        }

        new_col = f'rolling_{method}'

        for key in calculations.keys():
            df[f'rolling_{key}'] = calculations[key]

        # df[new_col] = calculations[method]
        df = df.dropna()
        df = df.sort_values(by=new_col, ascending=True)

        # filter dataframe based on start and stop arguments
        if start is not None:
            if not isinstance(start, int) and not isinstance(start, np.int64):
                raise Exception('start postion must be integer')
            df = df[df['position']>=(start + target_length)]
        
        if stop is not None: 
            try:
                stop = int(stop)
            except:
                raise Exception('end position must be integer')
            df = df[df['position']<=stop]

        # find the target region that meets criteria
        print('\nAnalyzing possible targets..')
        print('May finish before progress bar is complete')
        for i in tqdm(range(len(df))):
            row = df.iloc[i] # row of interest in the dataframe
            quantification = round(float(row[new_col]), 2)
            end_position = int(row.position) #position is 1 ahead of index, so this will be inclusive 
            start_position = end_position - (target_length - 1) #the number at the index is included in the max calculation, so only the previous (target_length -1) bases will also be included
            target_sequence = sequence[start_position -1: end_position]
            gc_percent = gc(target_sequence)
            med_ent = round(float(row['rolling_median']), 2)
            mean_ent = round(float(row['rolling_mean']), 2)
            max_ent = round(float(row['rolling_max']), 2)

            # it has to meet the GC percent cutoff specified
            if gc_percent < gc_cutoff:
                continue

            #don't allow '-' or 'J' or 'Z' to be in the target sequence - this will filter out frequently unread regions at the beginning or end of reads and common insertions and deletions
            if any(i in target_sequence for i in ['-', 'J', 'Z']):
                continue

            try:
                target = {
                    'start_position': start_position, 
                    'end_position': end_position, 
                    'target_length': end_position - start_position, 
                    'sequence': target_sequence, 
                    'GC_percent': gc_percent, 
                    'median_entropy': med_ent, 
                    'mean_entropy': mean_ent, 
                    'max_entropy': max_ent, 
                    'method_used:': f'lowest {method} entropy'
                }
                return target
            except:
                continue
         

    def dotplot(self, sequence, word_size=15, pseudo_range=None, range=None):
        """ 
        Add dotplot onto the entropy graph. The dotplot is used to compare regions of similarity between two sequences

        Arguments: 
            sequence (SeqRecord, str) -- sequence to compare to the pseudo_sequence. This sequence will make up the y-axis of the dotplot
            word_size (int) -- word size to use in matching sub-sequences. Larger word_size will result in fewer chance matches (Default=15)
            pseudo_range (list, tuple) -- range of the pseudo sequence to use. If None, the whole sequence will be used. 
            range (list, tuple) -- range of the graph to display
        return: 
            Plotly figure
        """

        if pseudo_range != None:
            sequence1 = self.pseudo_sequence[pseudo_range[0]: pseudo_range[-1]]
        else: 
            sequence1 = self.pseudo_sequence

        fig1, dp_df = dotplot(sequence1, sequence, word_size=word_size)
        self.dotplot_df = dp_df
        
        window = 30
        df = self.df
        df['window'] = df.entropy.rolling(window).mean()
        
        fig = go.Figure()

        fig.add_trace(go.Scatter(x=df['position'], y=df['window'], name=f'entropy rolling average = {window}', yaxis='y2', opacity=.8))
        fig.add_trace(go.Scatter(x=df['position'], y=df['entropy'], name='entropy values', yaxis='y2', opacity=.7))
        fig.add_trace(go.Scatter(x=dp_df['sequence1'], y=dp_df['sequence2'], name='dotplot', yaxis='y1', mode='markers', marker=dict(size=4)))
    

        fig.update_layout(
        yaxis=dict(
            title='comparative sequence index', 
            side='left'
        ), 
        yaxis2=dict(
            title='entropy values', 
            side='right',
            overlaying='y'
        ))
        if range is not None:
            range_error = 'range must be a list of 2 values'
            if type(range) not in [list, tuple] or len(range) != 2:
                raise Exception(range_error)
            fig.update_xaxes(range=range)
            


        title = f'Entropy and dotplot for {self.num_records} sequences of {self.target_name}. Word size = {word_size}'
        fig.update_layout(title=title)

        return fig


    def largest_dotplot_gaps(self):
        """ 
        Find the largest gap in the x and y dimensions of the dotplot. These will be the largest areas of heterogeneity 
        Return x_gap(list), y_gap(list) 
        """

        try:
            df = self.dotplot_df
        except:
            raise Exception('You must call show_dotplot first to run the calculations')
        y_missing = find_missing_y(df)
        x_missing = find_missing_x(df)
        x_gap = largest_gap(x_missing)
        y_gap = largest_gap(y_missing)

        return x_gap, y_gap


    def homology_analysis(self, primer, include_mismatches=False):
        """ 
        Peform mismatch analysis and return dataframe with the homology percentage as is needed for FDA submission

        Arguments: 
            primer (str): sequence of the primer to analyze
            include_mismatches (bool) -- if True, the returned dataframe with have a mismatches column as well as homology (Default False)
        Return: 
            Pandas DataFrame with percent homology and percent of sequences associated with each homology level
        """
        mismatches = self.target_mismatch_counts(primer)
        primer_len = len(primer)
        
        mismatch_numbers = []
        percents = []
        for i in mismatches.index: 
            mismatch_numbers.append(i)
        for j in mismatches:
            percents.append(round(j, 2))

        # create df 
        df = pd.DataFrame({
            'Mismatches': mismatch_numbers, 
            'Percent of Sequences': percents, 
        })

        df['Percent Homology'] = df['Mismatches'].apply(lambda x: round((primer_len - x)/primer_len * 100, 2))
        if include_mismatches:
            return df
        else:
            return df[['Percent Homology', 'Percent of Sequences']]





class ProteinEntropy:


    def __init__(self, fasta_path, target_name='Target'):
        """Instantiating an Entropy class will automatically calculated the entropy for 
        each position on the provided alignment files. 

        :param target_name: name of organism and/or gene of interest
        :param fasta_path: path to the .fasta.aln file where the alignments are stored

        :return: None, class will be instantiated and output can be specified using methods
        """
        self.target_name = target_name
        self.fasta_path = fasta_path
        self.df, self.pseudo_sequence, self.mismatches = self._generate_df() #dataframe of entropy values
        

    def _get_fasta_stats(self, handle):
        """Parses fasta file and returns the number of records and the length of records

        :return num_records: number of alignments in the file
        :return record_len: length of the alignments 
        """
        num_records = 0
        sequences = []
        for _, seq in SimpleFastaParser(handle):
            num_records += 1
            sequences.append(seq)
        record_len = len(sequences[0])

        sequences = [find_insertion_deletion(i) for i in sequences]

        return num_records, record_len, sequences


    def _entropy(self, data):
        """Returns an entropy calculation from input data
        :param data: Pandas series

        :return: entropy value
        """
        data = pd.Series(data)
        counts = data.value_counts()
        return scipy.stats.entropy(counts)
    
    
    def show(self, window=30, range=None, title=None):
        """show graph of entropy values
        Arguments:
            window (int, Optional) -- The number of positions to include in the rolling average (Default 30)
            range (lst, optional) -- The x-axis range if the entire genome is not desired (Default None)
            title (str, optional) -- If title different from target_name in class instance, pass title to update graph title (Default None)
        Return:
            None. Graph will be shown in console
        """
        fig = self._generate_fig(window, range, title)
        fig.show()
    

    def write_file(self, file_name=None, window=30, range=None, title=None):
        """Write html file with the graph of entropy values
        Arguments:
            file_name (str, Optional) -- file name to write html graph (Default ./{target_name}_entropy.html)
            window (int, Optional) -- The number of positions to include in the rolling average (Default 30)
            range (lst, optional) -- The x-axis range if the entire genome is not desired (Default None)
            title (str, optional) -- If title different from target_name in class instance, pass title to update graph title (Default None)
        Return:
            None. HTML graph file will be written 
        """
        fig = self._generate_fig(window, range, title)
        if file_name == None:
            fig.write_html(f'{self.target_name}_entropy.html')
        else: 
            if not file_name.endswith('.html'):
                file_name = file_name.split('.')[0] + '.html'
            parent_dir = Path(file_name).parent
            os.makedirs(parent_dir, exist_ok=True)
            basename = os.path.basename(file_name)
            filename = os.path.join(parent_dir, basename)
            try:
                fig.write_html(filename)
            except:
                print("That is not a valid file name")


    def _generate_df(self):
        """Creates the dataframe of the entropy values and a pseudo sequence  of most common bases
        Output:
            df: dataframe of the entropy values
            pseudo-sequence: sequence of the most common bases at each position of the alignment
            mismatches: list of the number of mismatches at each position
        """
        # convert alignment to numpy array
        with open(self.fasta_path) as f:
            
            self.num_records, self.record_len, all_sequences = self._get_fasta_stats(f)
            print(f'\n\nlooking at {self.num_records} sequences of length {self.record_len}...')

        with open(self.fasta_path) as f:
            name_ls = []
            seq_array = np.empty((self.num_records, self.record_len), dtype=str)
            idx = 0
            print('\nProcessing alignments...\nThis may take a while')
            for name, seq in tqdm(SimpleFastaParser(f), total=self.num_records, unit='seq'):  
                name_ls.append(name)
                seq_array[idx, :] = list(find_insertion_deletion(seq)) #inserted function here to change '-' to 'J' when it represents an insertion or deletion
                idx += 1
            self.seq_array = seq_array #save as class attribute so it can be accessed later
        
        #calculate entropy for each position
        positions = []
        entropy_values = []
        pseudo_sequence = []
        mismatches = []
        print('\ncomputing entropy...')
        for i in tqdm(np.arange(self.record_len), unit='pos', desc='Genome'):
            positions.append(i + 1)
            entropy_values.append(self._entropy([el for el in seq_array[:, i] if el != '-']))
            pseudo_sequence.append(self._common_base([el for el in seq_array[:, i] if el not in ['O', 'U']]))
            mismatches.append(self._count_mismatches([el for el in seq_array[:, i]]))
        #create dataframe of entropy values and rolling average values
        df = pd.DataFrame({'position': positions, 'entropy': entropy_values})
        return df, ''.join(pseudo_sequence), mismatches


    def _generate_fig(self, window, range=None, title=None):
        """ generate the plotly figure object with genome location on the x axis and entropy values on the y axis
        Parameters:
            window (int) -- The number of positions to include in the rolling average
            range (lst, optional) -- The x-axis range if the entire genome is not desired (Default None)
            title (str, optional) -- If title different from target_name in class instance, pass title to update graph title (Default None)
        Return:
            plotly graph object
        """
        df = self.df
        df['window'] = df.entropy.rolling(window).mean()
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=df.position, y=df.entropy, mode='lines', name='Entropy Values'))
        fig.add_trace(go.Scatter(x=df.position, y=df.window, line=dict(color='red', width=3), name=f"Rolling Average = {window}", showlegend=True))
        fig.update_layout(title=self.target_name,
                        xaxis_title='genome position',
                        yaxis_title='Entropy value')
        if range is not None:
            if type(range) != list:
                raise Exception('range type must be list with 2 values')
            fig.update_xaxes(range=range)
        if title is not None:
            fig.update_layout(title=title)
        return fig

    
    def _common_base(self, data):
        """get the most common base for the position passed
        """
        try:
            return Counter(data).most_common(1)[0][0]
        except IndexError:
            return '-'

    
    def _count_mismatches(self, data):
        """returns the number of bases not equal to the most common base
        """
        most_common = Counter(data).most_common(1)[0][0]
        return sum(1 for i in data if i not in [most_common, '-'])


    def get_sequence(self, start, end):
        """returns a pseudo sequence of the input alignment where each base is the most common 
        base at that position in the alignment

        Args:
            start: The starting position in the alignment for the sequence
            end: the ending position in the alignment for the sequence
        
        Output: 
            pseudo-sequence: sequence at the specified positions of the most common bases in the alignment
        """
        return self.pseudo_sequence[start:end+1]

    
    def mismatch(self, start, end):
        """analyze the mismatches 
        Args:
            start: The starting position in the alignment for the sequence
            end: the ending position in the alignment for the sequence
        Output:
            plotly histogram of the number of mismatches at the specified location
        """
        data = self.mismatches[start:end+1]
        position = np.arange(start, end+1)
        sequence = list(self.pseudo_sequence[start:end+1])
        df = pd.DataFrame(data=data, columns=['mismatches'])
        df['position'] = position
        df['base'] = sequence
        lower_range = min(df['mismatches']) - min(df['mismatches']) * .2
        upper_range = max(df['mismatches']) + max(df['mismatches']) * .2
        title = f'{self.target_name} {start} to {end} for {self.num_records} sequences'
        fig = px.bar(df, x='position', y='mismatches')
        fig.update_layout(title=title, 
            xaxis_title='Sequence', 
            yaxis_title='Number of Mismatches',
            xaxis=dict(
                nticks=len(data),
                tickvals = df['position'],
                ticktext = df['base']
            ),
            yaxis=dict(range=[lower_range, upper_range]

            ))
        fig.write_html(f'{self.target_name}_mismatch.html')


    def save_class(self, filename=None):
        if filename == None:
            filename = f'{self.target_name}_entropy_class.pkl'
        else: 
            filename = filename
        
        with open(filename, 'wb') as output:
            pickle.dump(self, output, -1)
        print(f'class saved as {filename}')  


    def find_target(self, target_length=100, percent=25, max_entropy=None, start=None, stop=None, method='max'):
        """
        Method to find target region with the lowest peak entropy to use as a reference sequence to analyze nucleic acid entropy

        Params:
            target_length (int) -- the length of the sequence to be used to get a reference nucleic acid sequence(Default 100)
            percent (int) -- the percentage of the data to use. Use a higher percentage if no results are given (Default 25)
            max_amplicon (int) -- the maximum allowed amplicon length from start of forward position to end of reverse position (Default 120)
            start (int) -- Starting position allowed for target -useful if designing on specific gene target (Default None)
            stop (int) -- Ending position allowed for target -useful if designing on specific gene target (Default None)
            method (str) -- 'max', 'mean', or 'median' method used to determine the best region. Max will give the region with the lowest maximum entropy (Default max)

        Return: 
            dict of sequence and position

        """
        # Initial error handling

        if percent <= 0 or percent > 100:
            raise Exception('percent = {} not allowed, percent must be a value between 0 and 100'.format(percent))

        if target_length < 100:
            print('A nucleic acid sequence < 300 may be less likely to produce a good region for CoPrimer design')
        
        if target_length > len(self.pseudo_sequence):
            print('the specified target length is larger than the sequence used in the analysis')
            target_length = len(self.pseudo_sequence)
        

        sequence = self.pseudo_sequence 

        df = self.df

        method = method.strip().upper()

        if method == 'MAX':
            df['rolling'] = df.entropy.rolling(window=target_length).max()
        elif method == 'MEAN':
            df['rolling'] = df.entropy.rolling(window=target_length).mean()
        elif method == 'MEDIAN':
            df['rolling'] = df.entropy.rolling(window=target_length).median()
        else:
            raise Exception('method got an unexpected value. Must be one of ["mean", "max", "median"]')
        


        df = df.dropna()
        df = df.sort_values(by='rolling', ascending=True)

        # filter dataframe based on start and stop arguments
        if start is not None:
            if type(start) != int:
                raise Exception('start postion must be integer')
            df = df[df['position']>=start]
        
        if stop is not None: 
            if type(stop) != int:
                raise Exception('end position must be integer')
            df = df[df['position']<=stop]

        # use only the lowest specified entropy values
        rows_to_keep = round(len(df) * percent/100)

        df = df.iloc[:rows_to_keep]

        # check max_entropy, if specified
        if max_entropy is not None: 
            df = df[df['rolling']<=max_entropy]

        # find the first target region that meets criteria
        print('\nAnalyzing possible targets..')
        for i in range(len(df)):
            row = df.iloc[i]
            end_position = int(row.position) #position is 1 ahead of index, so this will be inclusive 
            start_position = end_position - (target_length - 1) #the number at the index is included in the max calculation, so only the previous (target_length -1) bases will also be included
            target_sequence = sequence[start_position -1: end_position]


            #don't allow '-' or 'J' to be in the target sequence - this will filter out frequently unread regions at the beginning or end of reads and common insertions and deletions
            if '-' in target_sequence or 'J' in target_sequence:
                continue
            entropy_metric = row['rolling']
            print(f'{method.lower()} entropy: {entropy_metric}\nStarting position: {start_position}')

            target_dict = {
                'sequence': target_sequence, 
                'method_used': method,
                'entropy': entropy_metric,
                'start_position': start_position,
                'end_position': end_position
            }

            return target_dict

        else:
            # if this is reached, it means a target couldn't be found with the given criteria
            raise Exception('Could not find a target matching those criteria. Try widening search criteria')


    def target_mismatch_counts(self, target_sequence):
        """ 
        Counts the number of mismatches that each sequence in the class alignment would have with the given target and returns a pandas Series count. Mismatches and gaps are both treated as an error
        Params:
            target_sequence: The target sequence to analyze, generally a target given from the find_targets() method. Can be str or dict
        Output:
            Pandas Series value_counts of the percentage of sequences with that number of 'errors' for the given target region 
        
        """

        #use all sequences in the alignment
        sequences = self.seq_array
        genome = self.pseudo_sequence

        if type(target_sequence) == dict:
            if 'sequence' in target_sequence.keys():
                target_sequence = target_sequence['sequence']
            else: 
                raise Exception("'sequence' is not one of the dictionary keys")


        elif type(target_sequence) != str:
            raise Exception('Unrecognized data type for target_sequence, must be string') 

        target_sequence = target_sequence.strip().upper()


        start = genome.find(target_sequence)

        if start == -1:
            raise Exception('that target sequence is not in the sequence of most common amino acids')
        else:
            range = (start, start + len(target_sequence))

        errors = []

        for seq in tqdm(sequences):
            seq = ''.join(seq[range[0]: range[-1]])
            if seq in genome: #if the sequence is found, there are no errors
                errors.append(0)
                continue

            num_errors = 0
            for a, b in zip(seq, target_sequence):
                if a == '-':
                    continue
                elif a != b:
                    num_errors +=1 #treat gaps and mismatches both as errors
            errors.append(num_errors)
        
        errors_series = pd.Series(errors)

        return errors_series.value_counts(normalize=True, dropna=False).sort_index() * 100


    def write_target(self, file_name, target_length=100, percent=25, max_entropy=None, start=None, stop=None, method='max'):
        """ Method to find targets and write them to a text file. Combines find_targets method and write_target_dict function
        Arguments: 
            file_name (str/Path) -- location where the target will be written
            target_length (int) -- the length of the sequence to be used in design of the CoPrimers (Default 50)
            percent (int) -- the percentage of the data to use. Use a higher percentage if no results are given (Default 25)
            max_entropy (int) If specified, only targets with a maximum entropy value below the specified value will be returned. Values range from 0 to 1 (Default None)
            start (int) -- Starting position allowed for target -useful if designing on specific gene target (Default None)
            stop (int) -- Ending position allowed for target -useful if designing on specific gene target (Default None)
            method (str) -- 'max', 'mean', or 'median' method used to determine the best region. Max will give the region with the lowest maximum entropy (Default max)

        Return: 
            targets dict. text files will be written
         """

        if not file_name.endswith('.txt'):
            file_name = file_name.split('.')[0] + '.txt'

        target = self.find_target(target_length=target_length, percent=percent, max_entropy=max_entropy, start=start, stop=stop, method=method)
        if len(target) == 0:
            raise Exception('No targets to write, try widening parameters')
        write_target_dict(target_dict=target, file_name=file_name)

        return target


    def depth(self, normalize=True):
        """ 
        Graph the quality of sequences (proportion of sequences that were sequenced successfully) along the span of the sequence. 

        Arguments:
            seq_array (ndarray) -- a numpy array of aligned sequences
            normalize (bool, optional) -- normalize the depth counts (default True)
        Returns:
            Plotly figure 
        """  

        return graph_depth(self.seq_array, normalize=normalize)
