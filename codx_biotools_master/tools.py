import numpy as np
import json
import os
import glob
import openpyxl
import subprocess as sp
import pandas as pd
import pickle

from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
from numpy.random import default_rng

def replace_all(text, dic):
    """Replaces characters in a string using a dictionary"""
    for i, j in dic.items():
        text = text.replace(i, j)
    return text


def running_average_cumsum(seq, window=100):
    s = np.insert(np.cumsum(seq), 0, [0])
    return (s[window :] - s[:-window]) * (1. / window)


def slidingWindow(sequence, winSize, step=1):
    """Returns a generator that will iterate through
    the defined chunks of input sequence.  Input sequence
    must be iterable."""
    winSize = int(winSize)
    # Verify the inputs
    try:
        it = iter(sequence)
    except TypeError:
        raise Exception("**ERROR** sequence must be iterable.")
    if not ((type(winSize) == type(0)) and (type(step) == type(0))):
        raise Exception("**ERROR** type(winSize) and type(step) must be int.")
    if step > winSize:
        raise Exception("**ERROR** step must not be larger than winSize.")
    if winSize > len(sequence):
        raise Exception(
            "**ERROR** winSize must not be larger than sequence length.")

    # Pre-compute number of chunks to emit
    numOfChunks = int(((len(sequence) - winSize) / step) + 1)

    # Do the work
    for i in range(0, numOfChunks * step, step):
        yield sequence[i:i + winSize]

            
def remove_dir(dir_path):
    """ Will check if directory exists, then remove the directory and all its contents 
    Arguments
        dir_path (str) -- path to the directory
    Return
        None

    """
    if os.path.isdir(dir_path):
        files = glob.glob(dir_path + '/**/*', recursive=True)
        sub_dirs = []
        for file in files:
            if os.path.isfile(file):
                os.remove(file)
            else:
                try:
                    os.rmdir(file)
                except:
                    sub_dirs.append(file)

        if len(sub_dirs) > 0:
            for dir in sub_dirs:
                os.rmdir(dir)
        os.rmdir(dir_path)

        
def get_sheet_names(path):
    """ 
    Get the sheet names of an Excel workbook

    Arguments
        path (str, obj) -- path to the Excel workbook
    Reutrn:
        list of sheet names    
    """
    wb = openpyxl.load_workbook(path)
    return wb.sheetnames


def value_counts(iter, return_dict=False):
    """ 
    Count the number of occurences for each unique value in an interable object

    Arguments:
        iter (iterable) -- Any iterable object such as list or tuple
        return_dict (Bool) -- if True, a dictionary of unique values and the number of occurances will be returned. If False, nothing is returned (Default False)
    
    Return:
        if return_dict is True - Return a dictionary, else no return
    """
    unique_vals = []
    counts = {}
    for i in iter:
        if i not in unique_vals:
            unique_vals.append(i)
    for j in unique_vals:
        counts.update({j: iter.count(j)})
    
    counts = dict(sorted(counts.items(), key=lambda x: x[1], reverse=True))
    for key, val in counts.items():
        print(f'{key}: {val}')
    if return_dict:
        return counts


def to_clipboard(text):
    """ 
    Copy text to computer clipboard *Note: written specifically for MacOS. Other operating systems may behave differently
    Arguments:
        text (str) -- any string to copy to the clipboard
    Return:
        No return - text will be copied to the machine clipboard and can be pasted in any application
    """
    text = str(text)
    sp.run('pbcopy', universal_newlines=True, input=text)


def print_dict(dict):
    """ 
    Use json to print an easy-to-read layout of a dictionary

    Arguments:
        dict (dict) -- any dictionary 
    Return: 
        None
    """

    print(json.dumps(dict, sort_keys=True, indent=3))


def rand_subset(df, length):
    """ 
    Get a random subset of a pandas DataFrame of a specified length. This is useful to reduce the amount of data to make it more manageable but not introduce 
    selection bias from temporal or geographic trends in the data. 

    The data returned will only be reduced in axis=0, no column data will be eliminated

    Arguments:
        df (pandas DataFrame) -- Dataframe of any dimensionality
        length (int) -- length of the resulting DataFrame - if larger than input data, no change will be made

    Return:
        (DataFrame) -- Data with reduced rows -- a random subset of the input data
    """
    rng = default_rng()
    if isinstance(df, pd.DataFrame): 
        pass
    else:
        raise TypeError('df must be of type pandas.DataFrame')

    if length >= len(df):
        print('The length passed is larger than the input data, original DataFrame returned')

    ints = rng.choice(len(df), size=length, replace=False)
    return df.iloc[ints].copy(deep=True)


def load_class(path):
    with open(path, 'rb') as f:
        return pd.read_pickle(f)


def parent_dir(path):
    """ 
    Takes path to a file and returns the path to the directory that file is in

    arguments: 
        path (str, path-like): path to file. If path is a directory, the parent of the directory will be given

    returns: 
        path to the parent directory 
    """
        
    
    return '/'.join(path.split('/')[:-1])


def flatten(obj):
    """ 
    Flatten a 2-dimensional list or tuple and return a 1-dimensional list of the values
    """
    return [j for i in obj for j in i]

def read_target_dict(target_path):
    """ 
    Read the target file from write_target_dict and return the target dictionary
    """

    def format_dict(line):
        """ 
        inner function to format dictionaries read from text file
        """
        key, val = line.split(':')
        key, val = key.strip(), val.strip()
        if '.' in val:
            dtype = float
        else:
            dtype = str
        if dtype == float:
            return {key: float(val)}
        else:
            try:
                return {key: int(val)}
            except:
                return {key: val}

    forward = reverse = probe = False
    target_dict = {}
    with open(target_path, 'r') as f:

        for i, line in enumerate(f):
            if 'forward' in line.lower():
                forward = i
                forward_dict = {}
                continue
            elif 'reverse' in line.lower():
                reverse = i
                reverse_dict = {}
                continue
            elif 'probe' in line.lower():
                probe = i
                probe_dict = {}
                continue

    with open(target_path, 'r') as f:
        for i, line in enumerate(f):
            if i > forward and i < reverse and line.startswith('\t'):
                forward_dict.update(format_dict(line))
            if (i > reverse and probe == False and line.startswith('\t')) or (i > reverse and i < probe and line.startswith('\t')):
                reverse_dict.update(format_dict(line))
            if probe != False and i > probe and line.startswith('\t'):
                probe_dict.update(format_dict(line))
            if 'max_amplicon' in line:
                amplicon = format_dict(line)
    dict_names = ['forward', 'reverse']
    dicts = [forward_dict, reverse_dict]
    if probe != False:
        dicts.append(probe_dict)
        dict_names.append('probe')                
    for dictionary, name in zip(dicts, dict_names): 
        target_dict.update({name: dictionary})
    if amplicon:
        target_dict.update(amplicon)
    return target_dict
