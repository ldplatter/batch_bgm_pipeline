from Bio import Phylo

master_tree = Phylo.read("202_taxa_set_tree.nwk", "newick")


def format_taxa_name(taxa_name):
    taxa_name = taxa_name.lower()
    taxa_name = taxa_name.replace(" ", "_")
    taxa_name = taxa_name.capitalize()
    return(taxa_name)


shared_taxa_names = []
with open("shared_taxa.csv", "r") as file:  # This "shared_taxa.csv" need to be changed, as the name of the csv will change based on different proteins. 
                                            # See modified concat_fasta.py
    reader = file.readlines()
    for row in reader:
        shared_taxa_names += [format_taxa_name(row[:-1])]
    

all_taxa = []
for clade in master_tree.find_clades():
    if clade.is_terminal():
        all_taxa += [format_taxa_name(clade.name)]
        

taxa_to_prune = []
for taxa in all_taxa:
    if not(taxa in shared_taxa_names):
        taxa_to_prune += [taxa]
        

for taxa in taxa_to_prune:
    master_tree.prune(taxa)
Phylo.write(master_tree, "shared_taxa_tree.nwk", "newick")

 
