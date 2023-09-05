import os
import sys
import subprocess
from Bio import Phylo, SeqIO
import glob
import pexpect
import pandas as pd
import json
from networkx import *
import matplotlib.pyplot as plt
import concurrent.futures
from copy import deepcopy

'''
The structure of the command is: "python3 run_bgm.py protein1.fa protein2.fa [taxa_number/NA]".

protein2 could be either fasta file or a directory, in which case all fasta files in that directory will be 
concatenated with protein1 and run BGM tests respectively.

This modified script takes two user-defined raw protein fasta files as input, outputs the aligned fasta file,
both single and concatenated, the shared taxa list as a CSV file, and the pruned phylogenetic tree Newick file.
If taxa number equals to 108, then the 108 mammalian taxa phylogenetic tree was used to determine the final output tree.

If this argument is ignored, then the new 202 taxa master tree was used.

The script then carry out BGM utilizing HyPhy package automatically. It creates a HyPhy interactive session and 
automatically run BGM tests using GTR model and input data (generated prior to running BGMs).

The substitution model and further parameters could be also modified by editing parameters inside the call_bgm function.
'''
batch = False
batch_directory = ""
make_graph = False
protein1_length = 0
num_threads = 10


def format_taxa_name(taxa_name):
    taxa_name = taxa_name.lower()
    taxa_name = taxa_name.replace(" ", "_")
    taxa_name = taxa_name.capitalize()
    return(taxa_name)


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


def call_bgm(f_file, tree):
    command = "hyphy -i"
    child = pexpect.spawn(command)
    child.sendline("4")  # Coevolution
    child.sendline("1")  # BGM
    child.sendline("2")  # Amino-Acid
    child.sendline(f_file)
    child.sendline("11")  # GTR
    child.sendline(tree)
    child.sendline("1")  # All
    child.sendline("100000")  # MCMC Steps to sample
    child.sendline("10000")  # MCMC Steps to burn
    child.sendline("100")  # Steps to extract from the chain sample
    child.sendline("1")  # Maximum parents allowed per node
    child.sendline("1")  # Minimum substitutions per site
    child.interact()
    print(f"{f_file} BGM completed.")


def process_files(p1, p2, tree):
    # Output file names, change to "aligned_{os.path.basename(p1)}" if using the align_fasta function.
    aligned_protein1 = f"{os.path.basename(p1)}"
    aligned_protein2 = os.path.join(os.path.dirname(p2), f"{os.path.basename(p2)}")

    Protein1_fasta = open(aligned_protein1, "r").readlines()
    Protein2_fasta = open(aligned_protein2, "r").readlines()
    Protein1_dict = make_seq_dict(Protein1_fasta)
    Protein2_dict = make_seq_dict(Protein2_fasta)

    species_list = [i for i in Protein1_dict.keys() if i in Protein2_dict.keys()]

    concat_dict = {}

    for i in species_list:
        taxa = format_taxa_name(i)
        concat_dict[taxa] = Protein1_dict[taxa] + Protein2_dict[taxa]

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
    base_name1 = os.path.splitext(os.path.basename(p1))[0]
    base_name2 = os.path.splitext(os.path.basename(p2))[0]

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
    for clade in tree.find_clades():
        if clade.is_terminal():
            all_taxa += [format_taxa_name(clade.name)]

    taxa_to_prune = []
    for taxa in all_taxa:
        if not (taxa in shared_taxa_names):
            taxa_to_prune += [taxa]

    for taxa in taxa_to_prune:
        if len(sys.argv) > 3 and sys.argv[3] == "108":
            taxa = taxa.upper()
            tree.prune(taxa)
        else:
            tree.prune(taxa)

    shared_phylo_tree = f"shared_taxa_{base_name1}_{base_name2}.nwk"
    Phylo.write(tree, shared_phylo_tree, "newick")
    if batch:
        cmd = f"mkdir {base_name1}_{base_name2}"
        subprocess.call(cmd, shell=True)
        cmd = f"mv {base_name1}_{base_name2}/ {batch_directory}"
        subprocess.call(cmd, shell=True)
        cmd = f"mv {concat_fasta} {batch_directory}{base_name1}_{base_name2}/"
        subprocess.call(cmd, shell=True)
        cmd = f"mv {csv} {batch_directory}{base_name1}_{base_name2}/"
        subprocess.call(cmd, shell=True)
        cmd = f"mv {shared_phylo_tree} {batch_directory}{base_name1}_{base_name2}/"
        subprocess.call(cmd, shell=True)

        
    if batch:
        call_bgm(f"{batch_directory}{base_name1}_{base_name2}/{concat_fasta}", f"{batch_directory}{base_name1}_{base_name2}/{shared_phylo_tree}")
        cmd = f"mv {concat_fasta}.BGM.json {batch_directory}{base_name1}_{base_name2}/"
        subprocess.call(cmd, shell=True)
        
        bgm_df = initialize_results(f"{batch_directory}{base_name1}_{base_name2}/{concat_fasta}.BGM.json")
        
        
        probable_sites = find_probable_sites(bgm_df, 0.9)
        labeled_sites = add_within_between(probable_sites, protein1_length)
        labeled_sites.to_csv(f"{base_name1}_{base_name2}_90%.csv")
        cmd = f"mv {base_name1}_{base_name2}_90%.csv {batch_directory}{base_name1}_{base_name2}/"
        subprocess.call(cmd, shell=True)
        
        density, directionality = generate_summary(labeled_sites)
        summary_file = open(f"{base_name1}_{base_name2}_summary_90%.txt", "w")
        n = summary_file.write(f"{density} {directionality}")
        summary_file.close()
        cmd = f"mv {base_name1}_{base_name2}_summary_90%.txt {batch_directory}{base_name1}_{base_name2}/"
        subprocess.call(cmd, shell=True)
        
        probable_sites = find_probable_sites(bgm_df, 0.5)
        labeled_sites = add_within_between(probable_sites, protein1_length)
        labeled_sites.to_csv(f"{base_name1}_{base_name2}_50%.csv")
        cmd = f"mv {base_name1}_{base_name2}_50%.csv {batch_directory}{base_name1}_{base_name2}/"
        subprocess.call(cmd, shell=True)
        
        density, directionality = generate_summary(labeled_sites)
        summary_file = open(f"{base_name1}_{base_name2}_summary_50%.txt", "w")
        n = summary_file.write(f"{density} {directionality}")
        summary_file.close()
        cmd = f"mv {base_name1}_{base_name2}_summary_50%.txt {batch_directory}{base_name1}_{base_name2}/"
        subprocess.call(cmd, shell=True)
        
        
        if make_graph:
            create_graph(bgm_df, 0.9)
            cmd = f"mv {protein1}_{protein2}_connected_nodes.png {batch_directory}{base_name1}_{base_name2}/"
            subprocess.call(cmd, shell=True)
            
            cmd = f"mv {protein1}_{protein2}_all_nodes.png {batch_directory}{base_name1}_{base_name2}/"
            subprocess.call(cmd, shell=True)
            
    else:
        call_bgm(concat_fasta, shared_phylo_tree)
        bgm_df = initialize_results(f"{concat_fasta}.BGM.json")
        
        probable_sites = find_probable_sites(bgm_df, 0.9)
        labeled_sites = add_within_between(probable_sites, protein1_length)
        labeled_sites.to_csv(f"{base_name1}_{base_name2}_90%.csv")
        
        density, directionality = generate_summary(labeled_sites)
        summary_file = open(f"{base_name1}_{base_name2}_summary_90%.txt", "w")
        n = summary_file.write(f"{density} {directionality}")
        summary_file.close()
        
        probable_sites = find_probable_sites(bgm_df, 0.5)
        labeled_sites = add_within_between(probable_sites, protein1_length)
        labeled_sites.to_csv(f"{base_name1}_{base_name2}_50%.csv")
        
        density, directionality = generate_summary(labeled_sites)
        summary_file = open(f"{base_name1}_{base_name2}_summary_50%.txt", "w")
        n = summary_file.write(f"{density} {directionality}")
        summary_file.close()
        
        if make_graph:
            create_graph(bgm_df, 0.9)



# Post BGM file processing is below:
def length_of_first_sequence(fasta_file):
    for record in SeqIO.parse(fasta_file, "fasta"):
        return(len(record.seq))
    
def initialize_results(file):
    f = open(file)
    data = json.load(f)
    
    site1_list = []
    site2_list = []
    site1_to_site2_list = []
    site2_to_site1_list = []
    site1_or_site2_list = []
    site1_subs_list = []
    site2_subs_list = []
    shared_subs_list = []
    
    for i in data["MLE"]["content"]:
        site1_list += [i[0]]
        site2_list += [i[1]]
        site1_to_site2_list += [i[2]]
        site2_to_site1_list += [i[3]]
        site1_or_site2_list += [i[4]]
        site1_subs_list += [i[5]]
        site2_subs_list += [i[6]]
        shared_subs_list += [i[7]]
    
    bgm_df = pd.DataFrame(columns=["Site1", "Site2", "Site1 --> Site2", "Site1 <-- Site2", "Site1 <--> Site2", "Site1 Subs", "Site2 Subs", "Shared Subs"])

    bgm_df["Site1"] = site1_list
    bgm_df["Site2"] = site2_list
    bgm_df["Site1 --> Site2"] = site1_to_site2_list
    bgm_df["Site1 <-- Site2"] = site2_to_site1_list
    bgm_df["Site1 <--> Site2"] = site1_or_site2_list
    bgm_df["Site1 Subs"] = site1_subs_list
    bgm_df["Site2 Subs"] = site2_subs_list
    bgm_df["Shared Subs"] = shared_subs_list
    return(bgm_df)

def find_probable_sites(bgm_df, confidence=0.9):
    probable_sites_df = bgm_df[bgm_df["Site1 <--> Site2"] >= confidence]
    return(probable_sites_df)


#Creates graphs displaying bgm data.
def create_graph(bgm_df, confidence, protein_1, protein_2):
    prob_sites = find_probable_sites(bgm_df, confidence)
    g = Graph()
    g2 = Graph()
    nodes_added = []
    for i in range(0, len(prob_sites)):
        g.add_edge(int(prob_sites.iloc[i]["Site1"]), int(prob_sites.iloc[i]["Site2"]))
        g2.add_edge(int(prob_sites.iloc[i]["Site1"]), int(prob_sites.iloc[i]["Site2"]))
        if not(prob_sites.iloc[i]["Site1"] in nodes_added):
            nodes_added += [prob_sites.iloc[i]["Site1"]]
        if not(prob_sites.iloc[i]["Site2"] in nodes_added):
            nodes_added += [prob_sites.iloc[i]["Site2"]]
    for i in range(1,int(bgm_df["Site1"][bgm_df.index.values[-1]]) + 1):
        if not(i in nodes_added):
            g.add_node(i)
    
    nodes_with_edges = [node for node, degree in g.degree() if degree > 0]
    nodes_with_multi_edges = [node for node, degree in g.degree() if degree > 1]
    node_colors = ["red" if node in nodes_with_multi_edges else "red" if node in nodes_with_edges else "lightblue" for node in g.nodes()]
    
    pos = spring_layout(g, k=.3, scale=.1, center=(1,1), iterations=20)
    pos_nodes_with_edges = {node:position for (node, position) in pos.items() if node in nodes_with_edges}
    
    draw_networkx(g2, pos=pos_nodes_with_edges, node_color="red", node_size=150, font_size=6)
    plt.savefig(f"{protein1}_{protein2}_connected_nodes.png")
    
    draw_networkx(g, pos=pos, node_color=node_colors, node_size=150, font_size=6)
    plt.savefig(f"{protein1}_{protein2}_all_nodes.png")
    
    
# Determines if interactions are within one protein or between the two proteins, as well as determines direction of interaction.
def add_within_between(bgm_data, protein1_length):
    interaction_types = []
    direction_types = []
    for entry_ID in bgm_data.index:
        print
        if int(bgm_data["Site1"][entry_ID]) <= protein1_length and int(bgm_data["Site2"][entry_ID]) > protein1_length:
            interaction_types += ["Between"]
        else:
            interaction_types += ["Within"]
        if (float(bgm_data["Site1 --> Site2"][entry_ID])/float(bgm_data["Site1 <--> Site2"][entry_ID])) > 0.8:
            direction_types += ["1->2"]
        elif (float(bgm_data["Site1 <-- Site2"][entry_ID])/float(bgm_data["Site1 <--> Site2"][entry_ID])) > 0.8:
            direction_types += ["2->1"]
        else:
            direction_types += ["Bi"]
    bgm_data["Identity"] = interaction_types
    bgm_data["Direction"] = direction_types
    return(bgm_data)

# Generates summary data: Density and Directionality
def generate_summary(bgm_data):
    density = (len(bgm_data[bgm_data["Identity"] == "Between"]) / protein1_length)
    directionality = (len(bgm_data[bgm_data["Direction"] == "1->2"]) / len(bgm_data))
    return(density, directionality)


if __name__ == "__main__":
    protein1 = sys.argv[1]
    protein2 = sys.argv[2]
    protein1_length = length_of_first_sequence(protein1)
# Load the master tree, allowing for two phylo master tree options
    if len(sys.argv) > 3 and sys.argv[3] == "108":
        master_tree = Phylo.read("finished_mam_timetree.nwk", "newick")
    else:
        master_tree = Phylo.read("202_taxa_set_tree.nwk", "newick")

    if os.path.isdir(protein2):
        batch = True
        batch_directory = protein2
        protein2_files = [protein2 + filename for filename in os.listdir(protein2) if filename[-3:] == ".fa"]
        
        num_workers = min(num_threads, len(protein2_files))
        with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
            futures = [executor.submit(process_files, protein1, filename, deepcopy(master_tree)) for filename in protein2_files]
            print(f"Batch Started.")

    else:
        process_files(protein1, protein2, master_tree)
