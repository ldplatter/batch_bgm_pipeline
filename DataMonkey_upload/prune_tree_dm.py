import os
import argparse
from Bio import Phylo

'''

This code is used to prepare for datamonkey version of BGM.
It first prunes the complete phylogeny based on protien fasta file of input, then attach 
the pruned tree at the end of the fasta file.
The output file is ready for uploading to datamonkey website.
Execute the code as:
python3 prune_tree.py path_to_protein.fasta path_to_complete_tree.nwk
'''


def make_seq_dict(fasta):
    """Convert fasta to dictionary."""
    seq_dict = {}
    header = ""
    for line in fasta:
        if line.startswith(">"):
            header = line.strip().replace(">", "")
            seq_dict[header] = ""
        else:
            seq_dict[header] += line.strip()
    return seq_dict


def format_taxa_name(name):
    """Format taxa name."""
    return name.replace(" ", "_")


import os
from Bio import Phylo

def prune_tree(protein_fasta_path, complete_phylo_path):
    # Read the protein fasta file
    base_name = os.path.splitext(os.path.basename(protein_fasta_path))[0]
    with open(protein_fasta_path, "r") as f:
        protein_fasta = f.readlines()

    # Convert fasta to dictionary
    protein_dict = make_seq_dict(protein_fasta)

    # Get the list of taxa from the protein fasta
    protein_taxa = list(protein_dict.keys())

    # Read the complete phylogeny
    tree = Phylo.read(complete_phylo_path, "newick")

    # Get all taxa from the tree
    all_taxa = [format_taxa_name(clade.name) for clade in tree.find_clades() if clade.is_terminal()]

    # Determine which taxa to prune
    taxa_to_prune = [taxa for taxa in all_taxa if taxa not in protein_taxa]

    # Prune the tree
    for taxa in taxa_to_prune:
        tree.prune(taxa)

    # Convert the pruned tree to a string in Newick format
    tree_str = tree.format("newick")

    # Make a copy of the original fasta file and append the pruned tree
    appended_fasta_path = f"appended_{base_name}.fasta"
    with open(appended_fasta_path, "w") as f:
        for line in protein_fasta:
            f.write(line)
        f.write("\n")
        f.write(tree_str)

    print(f"Appended fasta file saved to: {appended_fasta_path}")

# Assuming you have the make_seq_dict and format_taxa_name functions defined elsewhere


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prune a phylogenetic tree based on taxa in a protein fasta file.")
    parser.add_argument("protein_fasta_path", help="Path to the protein fasta file.")
    parser.add_argument("complete_phylo_path", help="Path to the complete phylogenetic tree in Newick format.")

    args = parser.parse_args()
    prune_tree(args.protein_fasta_path, args.complete_phylo_path)
