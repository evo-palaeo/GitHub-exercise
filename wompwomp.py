
import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# -----------------------------
# Read sequences
# -----------------------------

file_path = "input.fas"
sequences = list(SeqIO.parse(file_path, "fasta"))

# -----------------------------
# Translate DNA sequences to proteins (longest ORF)
# -----------------------------

aa_sequences = []

for record in sequences:

    dna_seq = record.seq

    best_orf = ""
    best_frame = None

    # Check all three forward reading frames
    for frame in range(3):

        protein = dna_seq[frame:].translate(to_stop=False)

        # Split at stop codons
        fragments = str(protein).split("*")

        # Find longest fragment
        longest_in_frame = max(fragments, key=len, default="")

        # Update longest ORF
        if len(longest_in_frame) > len(best_orf):
            best_orf = longest_in_frame
            best_frame = frame

    # Store longest ORF as a Seq object
    best_protein = Seq(best_orf)

    # Extract organism name
    start = record.description.find("[organism=")

    if start != -1:
        start += len("[organism=")
        end = record.description.find("]", start)
        latin_name = record.description[start:end].replace(" ", "_")
    else:
        latin_name = record.id

    aa_record = SeqRecord(
        best_protein,
        id=record.id,
        name=latin_name,
        description=record.description,
        annotations={
            "type": "longest ORF",
            "frame": best_frame
        }
    )

    aa_sequences.append(aa_record)

# -----------------------------
# Write translated sequences
# -----------------------------

SeqIO.write(aa_sequences, "translated.fas", "fasta")

# -----------------------------
# Protein statistics
# -----------------------------

results = []

ambiguous_aas = set("XBZJUO*")

for record in aa_sequences:

    seq_str = str(record.seq)

    # Skip empty or ambiguous protein sequences
    if not seq_str or any(aa in ambiguous_aas for aa in seq_str):
        continue

    prot = ProteinAnalysis(seq_str)

    ambiguous_count = sum(aa in ambiguous_aas for aa in seq_str)

    aa_freqs = prot.get_amino_acids_percent()

    row = {
        "name": record.name,
        "id": record.id,
        "length": len(seq_str),
        "molecular_weight": prot.molecular_weight(),
        "isoelectric_point": prot.isoelectric_point(),
        "gravy": prot.gravy(),
        "ambiguous_count": ambiguous_count
    }

    row.update(aa_freqs)

    results.append(row)

# Create DataFrame
df = pd.DataFrame(results)

print(df.head())

# Save statistics
df.to_csv("protein_statistics.csv", index=False)

# -----------------------------
# END
# -----------------------------
