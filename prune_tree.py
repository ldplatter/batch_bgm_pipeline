import os
import argparse
from Bio import Phylo

'''
This code is a light version that only take a protein fasta file and complete phylogeny to output a pruned tree.
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

    # Output the pruned tree
    pruned_phylo_path = f"pruned_{base_name}.nwk"
    Phylo.write(tree, pruned_phylo_path, "newick")

    print(f"Pruned tree saved to: {pruned_phylo_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prune a phylogenetic tree based on taxa in a protein fasta file.")
    parser.add_argument("protein_fasta_path", help="Path to the protein fasta file.")
    parser.add_argument("complete_phylo_path", help="Path to the complete phylogenetic tree in Newick format.")

    args = parser.parse_args()
    prune_tree(args.protein_fasta_path, args.complete_phylo_path)
