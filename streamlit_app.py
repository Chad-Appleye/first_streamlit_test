py -m pip install -r requirements.txt

import streamlit as st
import pandas as pd
import bio

add_sidebar = st.sidebar.selectbox('Bioinformatics Tools', ('Alignments', 'Entropy Visualization',
                                                            'CoPrimer Selection Algorithm'))

# Alignments
if add_sidebar == 'Alignments':
    st.write("Alignments with Viral MSA")

    alignment_file = st.file_uploader(label='Alignment File', help="Upload a FASTA alignment file")
    primer_file = st.file_uploader(label='CoPrimer/Primer File', help="Upload a CoPrimer prediction file")
    
    
# Entropy Visualization
elif add_sidebar == 'Entropy Visualization':
    st.write("Entropy Visualized")

# CoPrimer Selection
elif add_sidebar == 'CoPrimer Selection Algorithm':
    st.write("CoPrimers Selected")
