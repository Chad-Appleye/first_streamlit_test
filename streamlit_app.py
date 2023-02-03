import streamlit as st
import pandas as pd
import Bio
from Bio import AlignIO
import io
from itertools import groupby

import design_analysis_tools_master.primer_mismatch as pm
import codx_biotools_master.oligotools as ot
from codx_biotools_master.tools import replace_all, slidingWindow, remove_dir, get_sheet_names, rand_subset, parent_dir, flatten
from Bio.SeqIO.FastaIO import SimpleFastaParser

# ************Functions**************

# def get_fasta_stats(alignment):
#   num_records = 0
#   for _, seq in SimpleFastaParser(alignment):
#     num_records += 1
#     record_len = len(seq)
#     return num_records, record_len

# def parse_alignment(alignment):
#   num_records, record_len = get_fasta_stats(alignment)
#   seq_array = np.empty((num_records, record_len), dtype=str)
  
#   name_ls = []
#   seq_ls = []
#   idx = 0
#   for name, seq in SimpleFastaParser(alignment):
#     name_ls.append(name)
#     seq_ls.append(seq)
#     idx += 1
#   return name_ls, seq_ls

def fasta_iter(fasta_name):
  """
  modified from Brent Pedersen
  Correct Way To Parse A Fasta File In Python
  given a fasta file. yield tuples of header, sequence
  """
  "first open the file outside "
  # fh = open(fasta_name)

  # ditch the boolean (x[0]) and just keep the header or sequence since
  # we know they alternate.
  faiter = (x[1] for x in groupby(fasta_name, lambda line: line[0] == ">"))

  for header in faiter:
      # drop the ">"
      headerStr = header.__next__()[1:].strip()

      # join all sequence lines to one.
      seq = "".join(s.strip() for s in faiter.__next__())

      yield (headerStr, seq)

  fiter = fasta_iter('testset.fas')

  for ff in fiter:
      headerStr, seq = ff
      print(headerStr)


###### Working on creating the Pseudo genome. Trying to figure out if I need to use any of the biopython tools

# def pseudo_genome(alignment):
#     # convert alignment to numpy array
#     num_records, record_len =  pm._get_fasta_stats
    
    
    
#     with open(self.alignment_path) as f:
#         self.num_records, self.record_len = self._get_fasta_stats(f)
#         print(f'Calculating mismatches for {self.num_records} records of length {self.record_len}')


#     with open(self.alignment_path) as f:
#         print('\nanalyzing sequences...')
#         name_ls = []
#         self.seq_array = np.empty((self.num_records, self.record_len), dtype=str)
#         idx = 0
#         for name, seq in tqdm(SimpleFastaParser(f), total=self.num_records):
#             name_ls.append(name)
#             self.seq_array[idx, :] = list(find_insertion_deletion(seq))
#             idx += 1

#     pseudo_genome = []
#     print('\nanalyzing genome..') 
#     for i in tqdm(np.arange(self.record_len)):
#         pseudo_genome.append(self._common_base([el for el in self.seq_array[:, i] if el in ['A', 'T', 'G', 'C']]))

#     pseudo_genome = ''.join(pseudo_genome)

#     return pseudo_genome


################ Building App/GUI Elements ###############################
add_sidebar = st.sidebar.selectbox('Bioinformatics Tools', ('Alignments', 'Entropy Visualization',
                                                            'CoPrimer Selection Algorithm'))

# Alignments
if add_sidebar == 'Alignments':
    st.write("Alignments with Viral MSA")
    
    alignment_file = st.file_uploader(label='Alignment File', help="Upload a FASTA alignment file")
    if alignment_file is not None:
      byte_str = alignment_file.read()
      text_obj = byte_str.decode('UTF-8')
      alignment = AlignIO.read(io.StringIO(text_obj),"fasta")
      
      st.write(type(alignment))
      for item in alignment:
        st.write(item)
       
      
      name_list, seq_list = fasta_iter(alignment)
      st.write(name_list)
      st.write(seq_list)
      
      
      
      
      
# Primers
    primer_file = st.file_uploader(label='CoPrimer/Primer File', help="Upload a CoPrimer prediction file")
    if primer_file is not None:
      if primer_file.name.endswith('.csv') or primer_file.name.endswith('.txt'):
        raw_primer_seq = pd.read_csv(primer_file)
      if primer_file.name.endswith('.xlsx'):
        raw_primer_seq = pd.read_excel(primer_file, engine='openpyxl')
      
      lin_coprimer_list = []
      
      for index, row in raw_primer_seq.iterrows():
        coprimer = ot.clean_sequence(row['Sequence'])
        gap = row['Gap']
        lin_coprimer_list.append(ot.linearize_coprimer(coprimer, gap))
        
      raw_primer_seq['LinearizedSeq'] = lin_coprimer_list
      
      forwards, reverses = ot.split_forwards_reverses(raw_primer_seq)
      forwards = forwards[['OligoName', 'Sequence', 'Gap','LinearizedSeq']].copy(deep = True)
      reverses = reverses[['OligoName', 'Sequence', 'Gap','LinearizedSeq']].copy(deep = True)
      
      
      st.write(forwards)
      st.write(reverses)
      
    
    if alignment_file is not None and primer_file is not None:
      kwargs = pm.mismatch_args(alignment_file,primer_file)
      visual = pm.MismatchVisualization(**kwargs)
      st.write("It might have worked")
        
    
# Entropy Visualization
elif add_sidebar == 'Entropy Visualization':
    st.write("Entropy Visualized")

# CoPrimer Selection
elif add_sidebar == 'CoPrimer Selection Algorithm':
    st.write("CoPrimers Selected")
