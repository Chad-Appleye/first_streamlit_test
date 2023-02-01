
import streamlit as st
import pandas as pd
import Bio
from Bio import AlignIO
import io

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

def read_primer_file(primer_file):
  df = ot.coprimers_to_dataframe

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
        
#       name_list, seq_list = parse_alignment(alignment)
#       st.write(name_list)
#       st.write(seq_list)
      
    primer_file = st.file_uploader(label='CoPrimer/Primer File', help="Upload a CoPrimer prediction file")
    if primer_file is not None:
      raw_primer_seq = pd.read_excel(primer_file, engine='openpyxl')
      primer_seq = raw_primer_seq[['Target', 'Sequence','OligoName', 'Gap']]
    st.write(raw_primer_seq)
    st.write(primer_seq)
      
      
      
      
      
    
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
