import os
import sys
import subprocess
from Bio import Phylo
import click

'''
The structure of the command is: "python3 concat_fasta.py protein1.fa protein2.fa [taxa_number/NA]".

This modified script takes two user-defined raw protein fasta files as input, outputs the aligned fasta file,
both single and concatenated, the shared taxa list as a CSV file, and the pruned phylogenetic tree Newick file.

If taxa number equals to 108, then the 108 mammalian taxa phylogenetic tree was used to determine the final output tree.
If this argument is ignored, then the new 202 taxa master tree was used.
'''


def format_taxa_name(taxa_name):
    taxa_name = taxa_name.lower()
    taxa_name = taxa_name.replace(" ", "_")
    taxa_name = taxa_name.capitalize()
    return taxa_name


def make_seq_dict(fasta_lines):
    seq_dict = {}
    current_seq = ""
    for i in range(len(fasta_lines)):
        if fasta_lines[i][0] == ">":
            current_seq = fasta_lines[i][1:-1]
            seq_dict[current_seq] = ""
        else:
            seq_dict[current_seq] += fasta_lines[i][:-1]
    return seq_dict


def split_by_60(string):
    count = 0
    current_string = ""
    lines = []
    for i in string:
        count += 1
        current_string += i
        if count == 60:
            lines += [current_string]
            count = 0
            current_string = ""
    if count != 0:
        lines += [current_string]
    final_string = ""
    for i in lines:
        final_string += i
        final_string += "\n"
    return final_string


def align_fasta(fastafile, outputfile, threads=5):
    """
    Aligns a file.

    :param threads: Number of threads
    :param fastafile: The input path
    :param outputfile: The output path
    """
    assert os.path.exists(fastafile), f"Input file {fastafile} does not exist."

    # Call MAFFT
    cmd = f"mafft --maxiterate 1000 --localpair --anysymbol --thread {threads} --out {outputfile} {fastafile}"
    subprocess.call(cmd, shell=True)

    assert os.path.exists(outputfile), f"Output file {outputfile} was not created."


# Takes two command line arguments, one for each protein
protein1 = sys.argv[1]
protein2 = sys.argv[2]

# Output file names
aligned_protein1 = f"aligned_{protein1}"
aligned_protein2 = f"aligned_{protein2}"

# Align the sequences
align_fasta(protein1, aligned_protein1)
align_fasta(protein2, aligned_protein2)

Protein1_fasta = open(aligned_protein1, "r").readlines()
Protein2_fasta = open(aligned_protein2, "r").readlines()

# Load the master tree, allowing for two phylo master tree options
if len(sys.argv) > 3 and sys.argv[3] == "108":
    master_tree = Phylo.read("finished_mam_timetree.nwk", "newick")
else:
    master_tree = Phylo.read("202_taxa_set_tree.nwk", "newick")

Protein1_dict = make_seq_dict(Protein1_fasta)
Protein2_dict = make_seq_dict(Protein2_fasta)

species_list = [i for i in Protein1_dict.keys() if i in Protein2_dict.keys()]

concat_dict = {}

for i in species_list:
    concat_dict[i] = Protein1_dict[i] + Protein2_dict[i]

out = ""
taxa_list = ""
for i in concat_dict.keys():
    out += ">"
    out += i
    out += "\n"
    out += split_by_60(concat_dict[i])
    taxa_list += i.replace("_", " ")
    taxa_list += "\n"

# New, Get the base names of the input fasta files (without the '.fa' extension)
base_name1 = os.path.splitext(protein1)[0]
base_name2 = os.path.splitext(protein2)[0]

concat_fasta = f"{base_name1}_{base_name2}.fa"
csv = f"{base_name1}_{base_name2}_taxa_list.csv"

with open(concat_fasta, "w") as f:
    f.write(out)

with open(csv, "w") as l:
    l.write(taxa_list)

# Read the taxa names from the output csv file
shared_taxa_names = []
with open(csv, "r") as file:
    reader = file.readlines()
    for row in reader:
        shared_taxa_names += [format_taxa_name(row[:-1])]

all_taxa = []
for clade in master_tree.find_clades():
    if clade.is_terminal():
        all_taxa += [format_taxa_name(clade.name)]

taxa_to_prune = []
for taxa in all_taxa:
    if not (taxa in shared_taxa_names):
        taxa_to_prune += [taxa]


for taxa in taxa_to_prune:
    if len(sys.argv) > 3 and sys.argv[3] == "108":
        taxa = taxa.upper()
    else:
        master_tree.prune(taxa)

shared_phylo_tree = f"shared_taxa_{base_name1}_{base_name2}.nwk"

Phylo.write(master_tree, shared_phylo_tree, "newick")
