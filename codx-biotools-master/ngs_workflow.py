import glob
from operator import pos
import pandas as pd 
import os
import matplotlib.pyplot as plt 
import subprocess
import numpy as np
import plotly.express as px
import allel
from tqdm import tqdm

from plotly.subplots import SubplotRef, make_subplots
from Bio.SeqIO.FastaIO import SimpleFastaParser


def make_fasta(csv_dir, output_file=None):

    """
    Generate a FASTA file from multiple csv files of sequence data. 

    Params: 
        csv_dir (str/PathLike) -- path to the directory containing the .csv files
        output_file (str/PathLike) -- path and name of the output file. Default is output.fasta

    Return: 
        None. FASTA file of the data contained in all of the csv files will be written 
    """

    all_files = glob.glob(os.path.join(csv_dir, "/**/*.csv"), recursive=True)

    df = pd.concat(pd.read_csv(f, header=None) for f in all_files)

    for i in range(len(df)):
        if i % 2:
            df.iloc[i] = df.iloc[i][0].replace('<', '').replace('>', '')

    if output_file:
        df.to_csv(output_file, sep='\n', index=False, header=False)
    else: 
        df.to_csv('output.fasta', sep='\n', index=False, header=False)


def get_references(fasta_file, names=None, outuput_dir=None, multi_fasta=True):
    """ 
    Makes reference fastas from the input fasta by replacing the IUPAC code and brackets at the SNP location with unambiguous bases
    Params: 
        fasta_file (str/PathLike) -- path to the fasta file that will be used to make the reference
        names (lst, optional) -- the names of the SNPs or targets to keep in the reference file, if None all will be kept (default None)
        output_file (str, optional) -- the path and/or name of the output fasta reference file (default None)
        multi_fasta (bool, optional) -- write a multi-fasta file with all references if True, write a fasta reference file for each reference if False (default True)

    Return: 
        pandas DataFrame with name and sequence. fasta reference file will be written
    """
    # IUPAC codes to replace with unambiguous bases
    replacements = {
        '[R]': 'A',
        '[Y]': 'C', 
        '[S]': 'C',
        '[W]': 'A',
        '[K]': 'G',
        '[M]': 'A',
        '[B]': 'C',
        '[D]': 'A',
        '[H]': 'A', 
        '[V]': 'A'
    }
    seqs = [] # list of dictionaries to build dataframe

    with open(fasta_file) as f:
        for name, seq in SimpleFastaParser(f):
            for key, value in replacements.items():
                seq = seq.replace(key, value)
            seqs.append({'name': name, 'sequence': seq})

    df = pd.DataFrame(seqs)
    
    # filter to only keep the references wanted 
    if names:
        df = df[df['name'].isin(names)]

    # output to new fasta file
    df['name'] = '>' + df['name'].astype(str)

    if outuput_dir is None: 
        outuput_dir = './references'
    os.makedirs(outuput_dir, exist_ok=True)


    if multi_fasta == True:
        df.to_csv(os.path.join(outuput_dir, 'multi_reference.fasta'), sep='\n', header=False, index=False)
    else:
        for index, row in df.iterrows():
            filename = row['name'].replace('>', '') + '.fasta'
            file_path = os.path.join(outuput_dir, filename)
            with open(file_path, 'w') as f:
                f.write(row['name'] + '\n' + row['sequence'])


    return df


def get_read_lists(dir_path, R1_identifier=None, R2_identifier=None):
    """ get the list of the read files separated into read1_list and read2_list. 
    Params: 
        dir_path (str/PathLike) -- path to parent directory containing all reads 
        R1_identifier (str, optional): string contained in all read_1 file names (default R1)
        R2_identifier (str, optional): string contained in all read_2 file names. (default R2) 
    Returns: 
        read_1 paths (lst), read_2 paths (lst) 
    """
    path = dir_path
    if R1_identifier is None:
        r1 = 'R1'
    else:
        r1 = str(R1_identifier).upper().strip()
    if R2_identifier is None:
        r2 = 'R2'
    else:
        r2 = str(R2_identifier).upper().strip()
    #get paths to all the files
    files = glob.glob(os.path.join(path, '**/*.fastq'), recursive=True)
    if len(files) == 0:
        files = glob.glob(os.path.join(path, '**/*.gz'), recursive=True)
    if len(files) == 0:
        raise Exception("couldn't find any .fastq or .gz files")
    read1_list = []
    read2_list = []

    for file in files:
        if r1 in file.upper():
            read1_list.append(file)
        elif r2 in file.upper():
            read2_list.append(file)
        else:
            raise Exception("couln't seperate paths into read lists, try different identifiers")
    read1_list.sort()
    read2_list.sort()
    
    return read1_list, read2_list


def read_calls(calls_dir, position=None, write_csv=True, separate=False):
    """Reads variant calls from multiple vcf format or compressed gz files and returns a DataFrame object of the variant calls for the specified position. If no variant identified for that position, 
    the variant column will be null. 
    Arugments:
        calls_dir (str/Path) -- Path to the directory for vcf files
        position (int, Optional) -- position of the variants to keep. If None, all will be kept (Default None)
        write_csv(bool, Optional) -- write the variant calls to a csv (Default True)
        separate(bool, Optional) -- if True, will write a csv for variant calls for each vcf file in the directory
    Return:
        pandas Dataframe with variants from specified position
    """    
    if not position is None:
        try:
            int(position)
        except ValueError:
            raise Exception('position must be an integer')

    # get the paths to the vcf files
    vcf_files = glob.glob(os.path.join(calls_dir, '*.vcf'))
    vcf_files.extend(glob.glob(os.path.join(calls_dir, '*.gz')))
    if len(vcf_files) == 0:
        raise Exception("No VCF files found")

    if separate is True:
        # output for all files
        out_dir = './variant_calls'
        for file in vcf_files:
            basename = os.path.basename(file).split('.')[0]
            df = allel.vcf_to_dataframe(file)
            if not position is None:
                df = df[df['POS']==position]
            df.drop_duplicates(subset=['CHROM', 'ALT_1'], inplace=True)
            df.dropna(axis=1, how='all', inplace=True)
            df.to_csv(os.path.join(out_dir, basename + '.csv'), index=False)
    
    df = pd.concat([allel.vcf_to_dataframe(i) for i in vcf_files])
    
    # format and drop duplicates
    if not position is None:
        df = df[df['POS'] == position]
    df.drop_duplicates(subset=['CHROM', 'ALT_1'], inplace=True)
    df.dropna(axis=1, how='all', inplace=True)
    

    if write_csv is True:
        filename = os.path.join(calls_dir, f'variants_pos_{position}.csv')
        df.to_csv(filename, index=False)
    
    return df 


def quality_filter(fastq_path_read1, fastq_path_read2, window_size=4, quality_score=15, output_dir=None):
    """ Trim reads using a sliding window and minimum quality score. Uses fastp trimming from 5' end, assuming paired-end reads that are g-zip compressed
    Params:
        fastq_path_read1 (str) -- path to read 1 fastq compressed file
        fastq_path_read2 (str) -- path to read 2 fastq compressed file 
        window_size (int, optional) -- size of the sliding window used for trimming (default 4)
        quality_score (int, optional) -- minimum acceptable quality score for sliding window (default 15)
    Returns:
        None, trimmed fastq files will be written in working_dir/filtered_reads/
        
        """
    if output_dir is None:
        out_dir = './filtered_reads'
    else:
        out_dir = os.path.join(output_dir, 'filtered_reads')
    os.makedirs(out_dir, exist_ok=True)

    report_dir = os.path.join(output_dir, 'fastp_reports')
    os.makedirs(report_dir, exist_ok=True)

    base1 = os.path.basename(fastq_path_read1)
    base2 = os.path.basename(fastq_path_read2)

    out1 = os.path.join(out_dir, base1)
    out2 = os.path.join(out_dir, base2)
    h = os.path.join(report_dir, base1.strip().upper().split('R1')[0] + '.fastp.html')
    j = os.path.join(report_dir, base1.strip().upper().split('R1')[0] + '.fastp.json')

    cmd = f'fastp -i {fastq_path_read1} -I {fastq_path_read2} -o {out1} -O {out2} -5 -W {window_size} -M {quality_score} -h {h} -j {j}'
    subprocess.call(cmd, shell=True)


def bcf_call(bam_file, reference, output_dir, variants_only=True, binary=False, compressed=False):
    """ use bcftools mpileup and bcftools call to generate a variant call file from a .bam file
    Parameters:
        bam_file (str/PathLike) -- path to the .bam file containing variants
        reference (str/PathLike) -- path to the .fasta reference file
        output_dir (str/PathLike) -- path to the output directory
        variant_only (bool, optional) -- if True, will show variant calls for only positions with valid variants, otherwise will show all (Default True)
        binary (bool, optional) -- if True, output file will be bcf binary (default False)
        compressed (bool, optional) -- if True, output file will be g compressed (default False)
    Return: 
        None. Variant call file will be written in output_dir/variant_calls
    """
    if output_dir is None:
        out_dir = 'variant_calls'
    else:
        out_dir = os.path.join(output_dir, 'variant_calls')
    os.makedirs(out_dir, exist_ok=True)

    if variants_only is True:
        variants = 'mv'
    else:
        variants = 'm'
    
    if binary == True:
        if compressed == True:
            out_type = 'b'
            extension = '.bcf.gz'
        else:
            out_type = 'u'
            extension = '.bcf'
    else:
        if compressed == True:
            out_type = 'z'
            extension = '.vcf.gz'
        else:
            out_type = 'v'
            extension = '.vcf'
    output_file = os.path.join(out_dir, os.path.basename(bam_file) + extension).replace('.bam', '')

    cmd = f'bcftools mpileup -Ou -f {reference} {bam_file} | bcftools call -{variants} -O{out_type} -o {output_file}'
    subprocess.call(cmd, shell=True)


def fastqc(read1, read2, output_dir=None):
    """ Use FastQC to generate quality information reports for NGS reads
    Params: 
        read1 (str/PathLike) -- path to NGS first read (fastq or .gz file)
        read2 (str/PathLike) -- path to second read (fastq or .gz file)
        output_dir (str/PathLike, optional) -- directory to store QC files. If not specified, a folder will be created in current working directory (default None) 
    Return:
        None. QC files will be written in specified directory
        """
    if output_dir is None:
        fastqc_dir = './fastQC'
    else:
        fastqc_dir = os.path.join(output_dir, 'fastQC')
    #make fastqc directory if it doesn't already exist
    os.makedirs(fastqc_dir, exist_ok=True)

    cmd = f'fastqc "{read1}" "{read2}" -o "{fastqc_dir}"'
    subprocess.call(cmd, shell=True)


def make_index(target, reference, output_dir=None):
    """ Create an index for the reads
    Parameters: 
        target (str) -- name of the target or read
        reference (str/PathLike) -- path to reference file
        output_dir (str/PathLike) -- path to dir where files will be written. If none passed, a directory will be created in current working directory
    Return:
        None. Files will be written in index_dir
    """
    #make index directory if necessary
    if output_dir is None:
        index_dir = './index'
    else: 
        index_dir = os.path.join(output_dir, 'index')

    os.makedirs(index_dir, exist_ok=True)

    # use bowtie2-build to make the index files
    cmd = f'bowtie2-build "{reference}" "{index_dir}/{target}"'
    subprocess.call(cmd, shell=True)


def align(read1, read2, target, output_dir=None):
    """ Use bowtie2 -x to create .sam file
    Parameters: 
        read1 (str/PathLike) -- path to read1 file
        read2 (str/PathLike) -- path to read2 file
        target (str) -- name of the allele, SNP, or sample
        output_dir (str/PathLike, optional) -- path to the output directory
    Return:
        (str) path to .sam file. .sam file will be created in output_dir/alignment
    """
    # make directory if not there
    if output_dir is None:
        alignment_dir = 'alignment'
    else: 
        alignment_dir = os.path.join(output_dir, 'alignment')
    if output_dir is None:
        index_dir = 'index'
    else:
        index_dir = os.path.join(output_dir, 'index')
    if not os.path.isdir(index_dir):
        raise Exception('No index directory. Cannot complete alignment without the index')

    os.makedirs(alignment_dir, exist_ok=True)
    sam_path = os.path.join(alignment_dir, f'{target}.sam')
    
    cmd = f'bowtie2 -x "{index_dir}/{target}" -1 "{read1}" -2 "{read2}" -S "{sam_path}"'
    subprocess.call(cmd, shell=True)
    return sam_path


def make_bam(target, sam_path, output_dir=None):
    """ 
    Uses samtools view and samtools sort to create binary .bam file from .sam file
    Parameters:
        target (str) -- name of the target gene, SNP, or sample
        bam_path (str/PathLike) -- path where .bam file will be written
        sam_path (str/PathLike) -- path to the .sam file used to create the .bam file
        output_dir (str/PathLike, optional) -- path to the directory where outputs will be written
    Return:
        (str) path to the .bam file. .bam file will be written in output_dir/alignment 
    """
    if output_dir is None:
        alignment_dir = 'alignment'
    else:
        alignment_dir = os.path.join(output_dir, 'alignment')
    
    os.makedirs(alignment_dir, exist_ok=True)
    
    if not os.path.isfile(sam_path):
        raise Exception(f'{sam_path} is Not a valid sam file')
    bam_path = os.path.join(alignment_dir, f'{target}.bam')
        
    cmd = f'samtools view -S -b "{sam_path}" | samtools sort -o "{bam_path}"'
    subprocess.call(cmd, shell=True)
    return bam_path


def calculate_depth(bam_path, output_dir=None):
    """ 
    Calculates the depth from a .bam file using samtools depth 
    Parameters:
        bam_path (str/PathLike) -- path to the .bam file 
        output_dir (str/PathLike) -- path to the output directory where files will be written
    Return:
        None. Depth files will be written in output_dir/depth 
    """
    if not os.path.isfile(bam_path):
        raise Exception(f'{bam_path} is not a valid path to a .bam file')
    if output_dir is None: 
        depth_dir = 'depth'
    else:
        depth_dir = os.path.join(output_dir, 'depth')

    os.makedirs(depth_dir, exist_ok=True)

    basename = os.path.basename(bam_path).split('.')[0]
    cmd = f'samtools depth -m 1000000 -a "{bam_path}" > "{os.path.join(depth_dir, f"{basename}.dep")}"'
    subprocess.call(cmd, shell=True)
    

def nonzero_median(x):
    return np.median(x[x != 0])


def nonzero_mean(x):
    return np.mean(x[x != 0])


def graph_depths(depth_dir, output_dir=None):

    if output_dir is None:
        plots_dir = 'plots'
    else:
        plots_dir = os.path.join(output_dir, 'plots')

    os.makedirs(plots_dir, exist_ok=True)

    for depth_path in glob.glob(os.path.join(depth_dir, '*.dep')):
        basename = os.path.basename(depth_path).split('.')[0]
        depth = pd.read_csv(depth_path, sep='\t', header=None)
        depth.columns = ['gene', 'position', 'depth']
        depth_group = depth.groupby('gene')
        stats = depth_group.agg({'depth': [np.median, np.mean, nonzero_median, nonzero_mean]})
        stats.columns = stats.columns.get_level_values(1)
        stats = stats.reset_index()
        stats_melt = pd.melt(stats, id_vars='gene')
        print(stats.sort_values('nonzero_mean', ascending=False)['gene'].unique())

        fig = px.bar(stats.sort_values('nonzero_mean', ascending=False), x='gene',
                    y='nonzero_mean', title='Mean depth ' + basename,
                    labels={'nonzero_mean': 'Mean coverage'})
        fig.show()
        fig.write_html(os.path.join(plots_dir, basename + '_MeanDepth.html'))

        fig = px.bar(stats.sort_values('nonzero_mean', ascending=False), x='gene',
                    y='nonzero_mean', log_y=True, title='Mean depth, log scale ' + basename,
                    labels={'nonzero_mean': 'Mean coverage'})
        fig.show()
        fig.write_html(os.path.join(plots_dir, basename + '_MeanDepthLog.html'))

        sorted_genes = stats.sort_values('nonzero_mean', ascending=False)['gene'].unique()
        fig = make_subplots(rows=4, cols=5, subplot_titles=[gene for gene in sorted_genes])
        plot_positions = [
            (1, 1), (1, 2), (1, 3), (1, 4), (1, 5),
            (2, 1), (2, 2), (2, 3), (2, 4), (2, 5),
            (3, 1), (3, 2), (3, 3), (3, 4), (3, 5),
            (4, 1), (4, 2), (4, 3), (4, 4), (4, 5)
        ]
        for gene, position in zip(sorted_genes, plot_positions):
            x = depth.loc[depth['gene'] == gene, 'position']
            y = depth.loc[depth['gene'] == gene, 'depth']
            fig.add_scatter(x=x, y=y, row=position[0], col=position[1],
                            showlegend=False, legendgroup='a', name='')
        fig.update_layout(title='Coverage depth ' + basename)
        fig.show()
        fig.write_html(os.path.join(plots_dir, basename + '_DepthByPosition.html'))


def workflow(read1, read2, reference, target, output_dir=None):
    """ Handles the main workflow of the creation of depth file. Called within main() for each read in the read lists. 
    Params:
        read11
         """
    fastqc(read1, read2, output_dir)
    make_index(target, reference, output_dir)
    sam_path = align(read1, read2, target, output_dir)
    bam_path = make_bam(target, sam_path, output_dir)
    calculate_depth(bam_path, output_dir)


def main(read1_list, read2_list, reference, variants_only=True, targets=None, output_dir=None, filter_reads=True, var_binary=False, var_compressed=False):
    """ Generate all necessary files to calculate and graph the depths of NGS reads
    Params:
        read1_list (lst) -- list of paths for all read1 files - can handle single path
        read2_list (lst) -- list of paths for all read2 files - can handle single path
        reference (str/PathLike) -- path to multi-fasta reference file - Use get_references func to generate references from the fasta file
        variants_only (bool, Optional) -- if True, variant call will only include variants, False will return all (Default True)
        targets: (lst, Optional) -- list of the targets - must be the same length as read lists. If None is passed, targets will be inferred from path names (Defualt None)
        output_dir (str/PathLike, optional) -- path to directory where subdirectories and output files will be written (Default None)
        filter_reads (bool, optional) -- filter reads using quality_filter and fastp. Default settings for fastp will be used, to change default settings, call quality_filter seperately before running main() (Default True)
        var_binary (bool, optional) -- If true, variant call file will be written in binary bcf. If false, vcf file will be written (Default False)
        var_compressed (bool, optional) -- If True, variant call file will be compressed to vcf.gz or bcf.gz (Default False)

    Return: 
        None. Depth graphs will show in console, and be saved to 'plots' folder in working directory. 
         """
    if output_dir is None:
        output_dir = './'
    os.makedirs(output_dir, exist_ok=True)

    # trim reads to remove poor quality regions 
    if filter_reads == True:
        filtered_reads_dir = os.path.join(output_dir, 'filtered_reads')
        os.makedirs(filtered_reads_dir, exist_ok=True)
        for read1, read2 in tqdm(zip(read1_list, read2_list)):
            quality_filter(read1, read2, output_dir=output_dir)

        read1_list, read2_list = get_read_lists(filtered_reads_dir)



    if type(read1_list) is list:
        if len(read1_list) != len(read2_list):
            raise Exception('The read1 and read2 lists are not the same length')

        if targets is not None: 
            if type(targets) != list:
                raise Exception('targets must be a list')
            if len(targets) != len(read1_list):
                raise Exception('The length of "targets" must be the same as the length of the read lists')
        else:
            #generate list of targets from the file names
            targets = []
            for file in read1_list:
                targets.append(os.path.basename(file).upper().split('R1')[0])

        for (read1, read2, target) in tqdm(zip(read1_list, read2_list, targets)):
            workflow(read1, read2, reference, target, output_dir)

    elif type(read1_list) is str:
        read1 = read1_list
        if type(read2_list) is not str:
            raise Exception('invalid dtype for read2_list')
        else:
            read2 = read2_list
        if targets is None:
            target = os.path.basename(read1_list).split('R1')[0]
        elif type(targets) is not str:
            raise Exception('invalid dtype for targets')
        else:
            target = targets
        workflow(read1, read2, reference, target, output_dir)
    else:
        raise Exception('invalid dtype passed for read_lists. The lists should be paths to fastq or .gz files')
    
    depth_dir = os.path.join(output_dir, 'depth')
    graph_depths(depth_dir, output_dir=output_dir)

    # variant call
    bam_files = glob.glob(os.path.join(output_dir, 'alignment/*.bam'))
    for file in bam_files:
        bcf_call(file, reference, output_dir=output_dir, variants_only=variants_only, binary=var_binary, compressed=var_compressed)
