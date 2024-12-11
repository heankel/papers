import pandas as pd
from collections import Counter
import random

#CURRENT FUNCTIONAL CODE
#input file contains 4 columns with titles:
#"protein_name"
#"protein_sequence"
#"aa_to_substitute"
#"aa_to_use_as_substitutes"

#row two for the purpose of our study looked like this:
#wt_TFE3_N
#MSHAAEPARDGVEASAEGPRAVFVLLEERRPADSAQLLSLNSLLPESGIVADIELENVLDPDSFYELKSQPLPLRSSLPISLQATPATPATLSASSSAGGSRTPAMSSSSSSRVLLRQQLMRAQAQEQERRERREQAAAAPFPSPAPASPAISVVGVSAGGHTLSRPPPAQVPREVLKVQTHLENPTRYHLQQARRQQVKQYLSTTLGPKLASQALTPPPGPASAQPLPAPEAAHTTGPTGSAPNSPMALLTIGSSSEKEIDDVIDEIISLESSYNDEMLSYLPGGTTGLQLPST
#AILV
#GFKF

#Note that to increase the possibilities of GOFs, some of the AILV amino acids were scrambled beforehand in certain iteractions of this code


def generate_mutant(sequence, substitutions):
    mutant_sequence = sequence[:]  # Make a copy of the sequence
    
    # Calculate the number of substitutions for each amino acid
    substitutions_count = {aa: sequence.count(aa) for aa in substitutions}
    
    for original_aa, substitute_aa in substitutions.items():
        # Calculate the number of substitutions to make for the current amino acid
        num_substitutions = substitutions_count[original_aa] // 2
        print(substitutions_count)
        print(num_substitutions)
        
        # Randomly select indices to start substitutions
        start_indices = random.sample(range(len(sequence)), num_substitutions)
        

        for index in start_indices:
            chosen_substitute = random.choice(substitute_aa)
            
            # Replace only the first occurrence at or after the selected index
            index_of_original = mutant_sequence.find(original_aa, index)
            if index_of_original != -1:
                mutant_sequence = mutant_sequence[:index_of_original] + chosen_substitute + mutant_sequence[index_of_original + 1:]
    
    return mutant_sequence

def calculate_aa_percentage(sequence):
    aa_count = Counter(sequence)
    total_aa = sum(aa_count.values())
    aa_percentage = {aa: count / total_aa * 100 for aa, count in aa_count.items()}
    return aa_percentage

# Read the Excel file
input_file = '/Users/S183959/Desktop/oncofusions_analysis/input_mutagenesis_0616.xlsx'
df = pd.read_excel(input_file)

# Get all unique amino acids from protein sequences
all_amino_acids = set(''.join(df['protein_sequence']))

# Get all unique amino acids used as substitutes
all_substitute_aa = set(df['aa_to_use_as_substitutes'].sum())

# Combine all unique amino acids
all_amino_acids.update(all_substitute_aa)

# Initialize aa_percentages dictionary with all possible amino acids
aa_percentages = {aa: [] for aa in all_amino_acids}

# Iterate over rows
mutant_sequences = []
for index, row in df.iterrows():
    protein_sequence = row['protein_sequence']
    substitutions = dict(zip(row['aa_to_substitute'], row['aa_to_use_as_substitutes']))

    # Generate mutants
    for i in range(30):  # Generate 20 mutants per row
        mutant_sequence = generate_mutant(protein_sequence, substitutions)
        mutant_sequences.append(mutant_sequence)

        # Calculate amino acid percentages
        aa_percentage = calculate_aa_percentage(mutant_sequence)

        # Fill in missing amino acids with 0 percentage
        for aa in all_amino_acids:
            if aa not in aa_percentage:
                aa_percentage[aa] = 0

        # Append percentages to aa_percentages
        for aa, percentage in aa_percentage.items():
            aa_percentages[aa].append(percentage)

# Write mutant sequences and amino acid percentages to Excel file
output_file = '/Users/S183959/Desktop/oncofusions_analysis/240616_output_mutagenesis_GOF_TFE3_AILV-GFKF.xlsx'
with pd.ExcelWriter(output_file) as writer:
    # Write mutant sequences
    pd.DataFrame({'Mutant Sequence': mutant_sequences}).to_excel(writer, sheet_name='Mutant Sequences', index=False)
    
    # Write amino acid percentages
    aa_percent_df = pd.DataFrame(aa_percentages)
    aa_percent_df.index.name = 'Mutant Index'
    aa_percent_df.to_excel(writer, sheet_name='Amino Acid Percentages')

print("Mutant sequences and amino acid percentages saved to", output_file)

