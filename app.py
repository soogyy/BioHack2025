import pandas as pd
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import streamlit as st

st.title("🧬 DNA Sequence Alignment and Cancer Risk Assesment 🧬")

DNA_sequence_1 = st.text_area("Please enter your DNA sequence: ")

uploaded_file = st.file_uploader("Upload cancer genes database (Excel file)", type=["xlsx"])
if uploaded_file:
    database = pd.read_excel(uploaded_file)
else:
    database = pd.read_excel("/Users/m._.nguyen/Documents/code/BioHack 2025/cancer genes.xlsx")

DNA_chars = set("AaCcGgTt")
RNA_chars = set("AaCcGgUu")
AA_chars = set("ARNDCQEGHILKMPSTWYV")

#DNA sequence
def validate_DNA_sequence(seq):
    seq = seq.upper().strip()
    if all(char in DNA_chars for char in seq):
        st.success("Searching for similar sequences...")
        return seq
    elif all(char in AA_chars for char in seq) and all(char not in RNA_chars for char in seq):
        st.error("❌ An amino acid sequence was entered. Please enter a DNA sequence.")
        return None
    elif any(char in RNA_chars for char in seq):
        st.error("❌ An RNA sequence was entered. Please enter a DNA sequence.")
        return None
    else: 
        st.error("❌ Invalid input. Please enter a DNA sequence.")
        return None

#scoring criteria
match_score = 2
mismatch_score = -1
gap_open_penalty = -1
gap_extend_penalty = -1

# Variables to keep track of the best match
best_match_percentage = 0
best_alignment_result = None
best_gene_name = None
best_condition = None
best_risk_chance = None

if DNA_sequence_1:
    with st.spinner("Aligning your sequence..."):
        for idx, row in database.iterrows():
            DNA_sequence_2 = str(row["DNA_seq"])
            alignments = pairwise2.align.localms(DNA_sequence_1, DNA_sequence_2, match_score, mismatch_score, gap_open_penalty, gap_extend_penalty)

            if alignments:
                best_alignment = alignments[0]
                aligned_seq_1, aligned_seq_2, score, begin, end = best_alignment

                # Count the number of matches (identical bases)
                matches = sum(1 for a, b in zip(aligned_seq_1, aligned_seq_2) if a == b and a != '-')
                
                # Calculate the alignment length (excluding gaps)
                alignment_length = sum(1 for char in aligned_seq_1 if char != '-')
                
                # Calculate the percentage of matches
                match_percentage = (matches / alignment_length) * 100

                # Update best match
                if match_percentage > best_match_percentage:
                    best_match_percentage = match_percentage
                    best_alignment_result = best_alignment
                    best_gene_name = row["Gene"]
                    best_condition = row["condition"]
                    best_risk_chance = row["Risk"]

if best_match_percentage >= 90:
    st.success("✅ Match found!")
    st.text(f"Gene associated with entered sequence: {best_gene_name}")
    
    # Risk assessment logic
    if best_risk_chance >= 4:
        st.text("Uh oh...")
        st.error(f"⚠️ You are at a HIGH risk of developing **{best_condition}**.")
    elif 3 >= best_risk_chance >= 2:
        st.text("Hmmm...")
        st.warning(f"You are at a moderate risk of developing **{best_condition}**")
    else:
        st.success("Good news!")
        st.success(f"✅ You have little to no risk of developing **{best_condition}**")
else:
    st.error("No match found in database :(")