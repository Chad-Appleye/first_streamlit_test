
import streamlit as st
import pandas as pd
import Bio
from Bio import AlignIO
import io

import design_analysis_tools_master.primer_mismatch as pm
import codx_biotools_master.oligotools
from codx_biotools_master.tools import replace_all, slidingWindow, remove_dir, get_sheet_names, rand_subset, parent_dir, flatten

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
      st.write(alignment)
  
    primer_file = st.file_uploader(label='CoPrimer/Primer File', help="Upload a CoPrimer prediction file")
    if primer_file is not None:
      primer_seq = pd.read_csv(primer_file)
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
