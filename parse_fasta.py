from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys

'''
The structure of the command is: python3 parse_fasta.py fasta_file output_file.
Takes an aligned fasta file as input, and remove all the gaps for Human species and corresponding positions 
for all other species.
'''

# The first argument will be the input fasta file, and the second argument will be the output file.
fasta_file = sys.argv[1]
output_file = sys.argv[2]

# Read all records
records = list(SeqIO.parse(fasta_file, "fasta"))

# Identify the Homo_sapiens (case insensitive)
homo_sapiens_record = None
for record in records:
    if record.id.upper() == "HOMO_SAPIENS":
        homo_sapiens_record = record
        break

if homo_sapiens_record is None:
    print("HOMO_SAPIENS not found in the fasta file.")
    sys.exit(1)

# Index of gap positions in Homo_sapiens sequence
gap_positions = [i for i, aa in enumerate(str(homo_sapiens_record.seq)) if aa == "-"]

# Remove gaps from Homo_sapiens sequence
homo_sapiens_record.seq = Seq(
    "".join([aa for i, aa in enumerate(str(homo_sapiens_record.seq)) if i not in gap_positions]))

# Create a new list of records
new_records = []

# Remove corresponding positions in other species
for record in records:
    if record.id != homo_sapiens_record.id:
        new_seq = "".join([aa for i, aa in enumerate(str(record.seq)) if i not in gap_positions])
        new_record = SeqRecord(Seq(new_seq), id=record.id, description="")
        new_records.append(new_record)
    else:
        new_records.append(homo_sapiens_record)

# Write the gapped sequences to the output fasta file
SeqIO.write(new_records, output_file, "fasta")
