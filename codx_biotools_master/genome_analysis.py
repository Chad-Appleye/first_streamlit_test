import os
import glob
import subprocess as sp
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import pickle

from entropy import Entropy, depth_counts
from Bio.SeqIO.FastaIO import SimpleFastaParser


def scan_large_genome(reference, sequences, email, window=100000, out_dir=None):
    """ 
    Calculate entropy in peice-wise method along a large genome to determine entropy values and most sequenced regions - ViralMSA.py must be in global path

    Uses ViralMSA to align each segment of the reference sequence. The entropy class for each segment will be saved. Once, a candidate region has been identified, load the saved class to 
    find targets in that region. 

    Arguments:
        reference (str) -- path to the reference FASTA of the whole genome
        sequnces (str) -- path to the multi-FASTA of all available sequences
        email (str) -- email to use with ViralMSA
        window (int, optional) -- the window size to use along the genome (default 100000)
        out_dir (str, optional) -- the output directory where all the files will be written. None will write a directory in the current working directory (default None)

    Return:
        pandas DataFrame object with depth and entropy for the full sequence
    
    """
    try:
        int(window)
    except:
        raise Exception('window must be an integer')

    if out_dir is None:
        out_dir = 'split_genome'
    

    # create subdir for references
    ref_dir = os.path.join(out_dir, 'references')
    os.makedirs(ref_dir, exist_ok=True)

    # load the reference sequence
    with open(reference, 'r+') as f:
        for name, seq in SimpleFastaParser(f):
            refseq = seq
            refname = name
    
    # split the reference sequence
    counter = 0
    while len(refseq) >= window:
        slice = refseq[:window]
        refseq = refseq[window:]
        out_path = os.path.join(ref_dir, f'reference_{counter}.fasta')
        name = f'> {refname}  reference_{counter}'
        with open(out_path, 'w') as r:
            r.write(f'{name} \n{slice}')
        counter += 1
    else:
        counter +=1
        name = f'> {refname}  reference_{counter}'
        out_path = os.path.join(ref_dir, f'reference_{counter}.fasta')
        with open(out_path, 'w') as r:
            r.write(f'{name} \n{refseq}')
    
    # run alignment

    alns = os.path.join(out_dir, 'alignments') # alignments dir
    references = glob.glob(f'{ref_dir}/*.fasta')
    os.makedirs(alns, exist_ok=True)

    for ref in references:
        basename = os.path.basename(ref).split('.')[0]
        out = os.path.join(alns, basename)
        cmd = f'ViralMSA.py -r {ref} -s {sequences} -e {email} -o {out}'
        sp.check_call(cmd, shell=True)
    
    # calculate the entropy for each alignment

    alignments = glob.glob(f'{alns}/**/*.aln', recursive=True)
    ent_dir = os.path.join(out_dir, 'entropies')
    os.makedirs(ent_dir, exist_ok=True)


    data = [] # store the dataframes

    for aln in alignments:
        basename = aln.split('/')[-2]
        ref_number = int(basename.split('_')[-1])
        start_position = ref_number * window + 1
        positions = np.arange(start_position, start_position + window, 1)
        ent = Entropy(aln, target_name=basename)
        folder = os.path.join(ent_dir, basename)
        os.makedirs(folder, exist_ok=True)
        positions = positions[:len(ent.df)]
        ent.df['position'] = positions # redefine the position in the entropy class
        index_adjust = window * ref_number
        ent.pseudo_sequence = '-' * index_adjust + ent.pseudo_sequence # adjust the index of the sequence
        ent.save_class(os.path.join(folder, basename + '.pkl'))
        ent.write_file(os.path.join('folder', basename + '.html'))
        data_dict = {
            'position': positions, 
            'depth': depth_counts(ent.seq_array, normalize=False),
            'entropy': ent.df['entropy']
        }
        tmp_df = pd.DataFrame(data_dict)
        data.append(tmp_df)
    df = pd.concat(data)
    df = df.sort_values(by='position')
    return df

    
def graph_entropy_depth(df, window=50):
    """ 
    Takes dataframe output from scan_large_genome and calculates the rolling entropy values, returns a plotly figure
    object showing the depth and entropy values graphed on seperate y-axes

    Arguments:
        df (pandas DataFrame) -- DataFrame object output from scan_large_genome
        window (int, optional) -- rolling average window for entropy values (default 50)
    Return:
        plotly figure object
    """

    df['rolling_entropy'] = df['entropy'].rolling(window).mean()

    fig = go.Figure()

    fig.add_trace(go.Scatter(x=df['position'], y=df['depth'], name='depth', yaxis='y1'))
    fig.add_trace(go.Scatter(x=df['position'], y=df['rolling_entropy'], name=f'entropy window = {window}', yaxis='y2', opacity=.8))

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
    
    return fig


def test_region(region, entropy_class, target_length=50, max_amplicon=200, gc_cutoff=40):
    """ 
    get statistics about test designs on a target
    Arguments:
        region (iterable of length 2) -- start and stop position of the region to test
        entropy_class (Entropy object) -- entropy class object or path to saved object with .pkl extension of corresponding entropy class
        target_length (int, optional) -- length of target to return (default 50)
        max_amplicon (int, optional) -- max_amplicon for the Entropy.find_targets method (defualt 200)
        gc_cutoff (int 0 to 100, optional) -- lowest allowed gc percent for a target (default 40)

    return:
        dictionary of targets. Statistics will be printed to the console
    """

    if type(entropy_class) == str:
        if entropy_class.endswith('.pkl'):
            with open(entropy_class, 'rb') as f:
                entropy_class = pickle.load(f)
    elif type(entropy_class) == 'entropy.Entropy':
        pass
    else:
        raise ValueError('Not a valid data type for entropy_class')
    
    if len(region) != 2:
        raise ValueError('region must be an iterable of 2 values')
    
    start = region[0]
    stop = region[-1]

    targets = entropy_class.find_targets(target_length=target_length, max_amplicon=max_amplicon, gc_cutoff=gc_cutoff, start=start, stop=stop)

    forward = targets['forward']
    reverse = targets['reverse']

    print(f'forward mismatches: \n{entropy_class.target_mismatch_counts(forward)}')
    print(f'reverse mismatches: \n{entropy_class.target_mismatch_counts(reverse)}')

    return targets

    

    

