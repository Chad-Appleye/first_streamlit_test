#!/Users/dallondurfey/opt/anaconda3/envs/codx/bin/python

import argparse
import os
from os.path import split
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pandas.core.indexes.base import Index
import plotly.express as px
import glob
import entropy


from sklearn.decomposition import PCA 
from sklearn.cluster import KMeans


def build_seq_matrix(aln_file):
    """ Create a sequence matrix where rows are individual sequences and columns are positions along the genome
    Params: 
        aln_file: path to a fasta.aln file -- generally output from ViralMSA 
        
    Output: 
        numpy array sequence matrix nxr where n is the number of sequences given and r is the length of the genome, 
        list of formatted sequences, 
        list of data about each corresponding sequence"""

    #put sequences in to a list
    if aln_file.endswith('.aln'):
        with open(aln_file, 'r') as f:
            sequences = f.read().splitlines()
    else:
        raise Exception('The file passed was not an alignment file')
    
    meta_data = []
    formatted_sequences = []
    matrix = []
    for seq in sequences:
        if seq[0] == '>':
            meta_data.append(seq)
            continue
        matrix.append(list(seq))
        formatted_sequences.append(seq)

    return np.array(matrix), formatted_sequences, meta_data


def perform_pca(seq_matrix, meta_data, n_components, plot, explained_variance, out_dir):
    """ Perform a principal componenet analysis (PCA) on a sequence matrix using one-hot-encoding for values at each position in the genome
    Params: 
        seq_matrix (ndarray) -- output from build_seq_matrix
        meta_data (list) -- data about each sequence in the matrix
        n_components (int) -- Number of components to use in the PCA
        plot (bool) -- True will produce a 2-D plot of the first 2 principal components
        explained_variance (bool) -- True will calculated the explained variance for the first 5 principal components and produce a line plot
        
    Output: 
        returns a pandas DataFrame object with each sequence and the principal components specified"""
        

    figures_folder = os.path.join(out_dir, 'figures')
    os.makedirs(figures_folder, exist_ok=True)
    
    #create pandas df
    seq_df = pd.DataFrame(seq_matrix)
    bool_df = pd.get_dummies(seq_df) #use one-hot-encoding to get a boolean matrix 

    #define the pca model
    pca = PCA(n_components=n_components, random_state=42)

    df_pca = pca.fit_transform(bool_df)

    columns = []
    for i in range(n_components):
        columns.append(f'PC{i+1}')

    PCA_df = pd.DataFrame(df_pca, columns=columns)   
    PCA_df['sequence'] = meta_data

    if explained_variance == True:
        #fit the model to a few more components to geta better idea of the explained variance
        pca2 = PCA(n_components=5, random_state=42)
        pca2.fit(bool_df)
        ev = pca2.explained_variance_ratio_
        print('Explained Variance: \n{}'.format(ev))
        df_ev = pd.DataFrame(ev, index=range(1, 6), columns=['Explained Variance'])
        ev_fig = px.line(df_ev, x=df_ev.index, y='Explained Variance')
        ev_fig.write_html(os.path.join(figures_folder, 'expalined_variance.html'))
        ev_fig.update_layout(
            xaxis = dict(
                tickmode = 'array', 
                tickvals = list(range(1, 6))
            )
        )
        ev_fig.show()



    if plot == True:
        if n_components != 2:
            print('Plotting only the first 2 principal componenets...')
        PCA_df.plot.scatter(x='PC1', y='PC2')
        plt.title('Principal Component Analysis')
        
    
    return PCA_df


def cluster(df, n_clusters, cluster_plot, out_dir):
    """ Use k-means clustering to group sequences into their subgroups
    Params: 
        df (DataFrame) -- pandas DataFrame object output from perform_pca(). 
        n_clusters (int) -- number of clusters to use in the kmeans model
        cluster_plot (bool) -- Create a cluster plot, same as the 2-D PCA plot, but the markers are colored by their kmeans group
        
    Output: returns DataFrame with sequences and their associated clusters """
    if type(df) != pd.core.frame.DataFrame:
        raise Exception('dataframe for cluster must be a pandas dataframe')

    kmeans = KMeans(n_clusters=n_clusters, random_state=42, verbose=1)

    #drop non-numberical data columns
    cluster_data = df.select_dtypes(['number'])

    kmeans.fit_transform(cluster_data)

    labels = kmeans.labels_

    #add labels to the df
    df['cluster'] = labels
    #convert to string dtype so that colorscale is discrete in plot
    df['cluster'] = df['cluster'].astype(str)


    if cluster_plot == True:
        figure_folder = os.path.join(out_dir, 'figures')
        os.makedirs(figure_folder, exist_ok=True)
        fig = px.scatter(df, x='PC1', y='PC2', color='cluster', hover_data=['sequence'], 
        title='Principal component analysis with K-Means clustering')
        fig.write_html(os.path.join(figure_folder, 'cluster_fig.html'))
        fig.show()
    return df


def split_aln(cluster_df, sequences, meta_data, out_dir):
    """ write a new fasta.aln file for each of the groups output from the cluster analysis
    Params:
        df: Pandas DataFrame object with the sequences, their data, and the cluster value assigned to it from K-means output from cluster()
        sequences: list of sequences from the original fasta.aln file
        meta_data: list of the sequence data from the fasta.aln file
    Ouput:
        returns path to the split_alignments folder for use further in the workflow. A file will be written for each cluster
         """
    alignments_folder = os.path.join(out_dir, 'split_alignments')
    if not os.path.isdir(alignments_folder):
        os.mkdir(alignments_folder)

    dfs = [x for _, x in cluster_df.groupby('cluster')]

    for df in dfs:
        cluster = df['cluster'].iloc[0]
        file_name = os.path.join(alignments_folder, 'cluster' + cluster + '.fasta.aln')
        with open(file_name, 'w+') as f:
            for seq, data in zip(sequences, meta_data):
                if data in df['sequence'].tolist():
                    f.write(data + '\n' + seq + '\n')
    
    print('alignments split')
    return alignments_folder


def cluster_entropies(alignments_path, organism, out_dir):
    """ 
    Calculate the entropy profile for each alignment file
    Params:
        alignments_path: list of paths or path to directory containing alignment files
        organism (str) -- the target organism that will appear in the heading of graphs created

    Output: 
        No return. html plotly graphs with entropy profiles will be written, entropy classes will be saved as binary and writen as .pkl files
    """
    if organism is None:
        organism = 'Target'

    if type(alignments_path) == list:
        clusters = alignments_path.sort()
    else:
        clusters = sorted(glob.glob(alignments_path + '/*.fasta.aln'))
    
    if len(clusters) == 0:
        raise Exception('no alignments found, must be .aln file')

    print('calculating entropy for the following files: {}'.format(clusters))

    entropies = os.path.join(out_dir, 'entropies')
    entropy_graphs = os.path.join(entropies, 'graphs')

    os.makedirs(entropy_graphs, exist_ok=True)

    for i in range(len(clusters)):
        cluster = clusters[i]
        cluster_ent = entropy.Entropy(cluster, f'{organism}_cluster_{i}')
        cluster_ent.save_class(os.path.join(entropies, f'cluster_{i}.pkl'))
        cluster_ent.write_file(os.path.join(entropy_graphs, f'cluster_{i}.html'))


def main(aln_file, n_components=2, pca_plot=False, explained_variance=True, n_clusters=3, cluster_plot=True, organism=None, out_dir=None):

    """ 
    Use unsupervized machine learning and principal component analysis to group sequences into mathematically similar groups

    Arguments:
        aln_file: path to a fasta.aln file -- generally output from ViralMSA 
        n_components (int) -- Number of components to use in the PCA (Default 2)
        pca_plot (bool) -- True will produce a 2-D plot of the first 2 principal components (Default False)
        explained_variance (bool) -- True will calculated the explained variance for the first 5 principal components and produce a line plot (Default True)
        n_clusters (int) -- number of clusters to use in the kmeans model (Default 3)
        cluster_plot (bool) -- Create a cluster plot, same as the 2-D PCA plot, but the markers are colored by their kmeans group (Defualt True)
        organism (str) -- the target organism that will appear in the heading of graphs created (Default None)
        out_dir (str) -- path to the output directory where all files will be written. If None, files and subdirectories will be created in working directory (Default None)

    """

    #make out_dir if necessary
    if out_dir is None:
        out_dir = './'
    os.makedirs(out_dir, exist_ok=True)

    #first, build the sequence matrix and get the lists of sequences and data
    print('building sequence matrix...')
    seq_matrix, sequences, meta_data = build_seq_matrix(aln_file)

    #perform pca on the sequences matrix
    print('calculating PCA...')
    PCA_df = perform_pca(seq_matrix, meta_data, n_components, pca_plot, explained_variance, out_dir)

    #perform k-means cluster analysis
    print('Performing clustering...')
    cluster_df = cluster(PCA_df, n_clusters, cluster_plot, out_dir)

    #split the alignment file into n_clusters number of files
    print('splitting alignment files...')
    alignments_folder = split_aln(cluster_df, sequences, meta_data, out_dir)

    #calculate entropy for each cluster
    cluster_entropies(alignments_folder, organism, out_dir)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='use PCA and clustering to find genetic subgroups from an alignment file. Output those as pickled entropy classes and entropy graphs, as well as the seperated .aln files for each group')
    parser.add_argument('aln', type=str, help='path to the alignment (.aln) file, generally output from ViralMSA')
    parser.add_argument('-comp', '--n_components', type=int, default=2, help='(Optional) n-components for the PCA model, default=2')
    parser.add_argument('--pca_plot', type=bool, default=False, help='(Optional) output a plot of the PCA in 2 dimensions without cluster coloring, default=False')
    parser.add_argument('-ev', '--explained_variance', type=bool, default=True, help='(Optional) Output a plot of the explained variance from the PCA model for the first 5 components')
    parser.add_argument('-clust', '--n_clusters', type=int, default=3, help='(Optional) number of centroids or clusters in the K-Means cluster analysis, default=3')
    parser.add_argument('--cluster_plot', type=bool, default=True, help='(Optional) Output a plot of the PCA with cluster colors, default=True')
    parser.add_argument('-O', '--organism', type=str, default=None, help='The name of the organism -- or gene marker -- targeted that will be in the title of the graphs and the name of the entropy classes')
    parser.add_argument('-o' '--out_dir', type=str, default=None, help='path to output directory')
    args = parser.parse_args()
    main(args.aln, args.n_components, args.pca_plot, args.explained_variance, args.n_clusters, args.cluster_plot, args.organism, args.out_dir)


