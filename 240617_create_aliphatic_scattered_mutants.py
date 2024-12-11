import pandas as pd
from localcider.sequenceParameters import SequenceParameters

# Function to calculate omega parameters
def calculate_omega_params(sequence, aliphatic_residues=['A', 'V', 'L', 'I']):
    S = SequenceParameters(sequence)
    omega_aliphatic = S.get_kappa_X(aliphatic_residues)
    return omega_aliphatic

# Input protein sequence
#PRCC_fusion as an example
protein_sequence = "MSLVAYASSDESEPDEAEPEPEEEEAVAPTSGPALGGLFASLPAPKGPALLPPPPQMLAPAFPPPLLLPPPTGDPRLQPPPPLPFGLGGFPPPPGVSPAEAAGVGEGLGLGLPSPRGPGLNLPPPIGGAGPPLGLPKPKKRKEPVKIAAPELHKGDSDSEEDEPTKKKTILQGSSEGTGLSALLPQPKNLTVKETNRLLLPHAFSRKPSDGSPDTKPSRLASKTKTSSLAPVVGTTTTTPSPSAIKAAAKSAALQVTKQITQEEDDSDEEVAPENFFSLPEKAEPPGVEPYPYPIPTVPEELPPGTEPEPAFQDDAANAPLEFKMAAGSSGAPWMPKPGDDYSYNQFSTYGDANAAGAYYQDYYSGGYYPAQDPALVPPQEIAPDASFIDDEAFKRLQGKRNRGREEINFVEIKGDDQLSGAQQWMTKSLTEEKTMKSFSKKKGEQPTGQQRRKHQITYLIHQAKERELELKNTWSENKLSRRQTQAKYGF"

# Step 1: Count A, I, L, V residues
aliphatic_residues = ['A', 'I', 'L', 'V']
count_aliphatic = sum(protein_sequence.count(res) for res in aliphatic_residues)
total_length = len(protein_sequence)

# Step 2: Calculate scattering ratio
percent_aliphatic = count_aliphatic / total_length
scattering_ratio = int(round(1 / percent_aliphatic)) if percent_aliphatic > 0 else 1

# Step 3: Remove and store aliphatic residues in order
aliphatic_removed = []
modified_sequence = []

for aa in protein_sequence:
    if aa in aliphatic_residues:
        aliphatic_removed.append(aa)
    else:
        modified_sequence.append(aa)

# Step 4: Scatter aliphatic residues back into the sequence
index = 0
for aa in aliphatic_removed:
    modified_sequence.insert(index, aa)
    index += scattering_ratio

# Join list into final protein sequence string
final_sequence = ''.join(modified_sequence)

# Step 5: Calculate omega aliphatic values
initial_omega_aliphatic = calculate_omega_params(protein_sequence)
final_omega_aliphatic = calculate_omega_params(final_sequence)

# Step 6: Output the final modified protein sequence
print("Initial Protein Sequence:")
print(protein_sequence)
print("\nModified Protein Sequence:")
print(final_sequence)
print("\nInitial Omega Aliphatic:", initial_omega_aliphatic)
print("Final Omega Aliphatic:", final_omega_aliphatic)


