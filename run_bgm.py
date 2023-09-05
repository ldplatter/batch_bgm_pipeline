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

Lastly, the script outputs summary statistics of the BGM test results and the network topology based on the input data,
using the networkx package.
'''
batch = False
batch_directory = ""
make_graph = True  # Default: True
protein1_length = 0


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


def call_bgm(f_file, tree):
    command = "hyphy -i"
    child = pexpect.spawn(command)
    child.sendline("4")  # Co-evolution
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


def write_summary_and_csv(bgm_df, protein_length, threshold, base_name1, base_name2, directory=None):
    # Write the summary files to multiple csv files and move them to new directories
    probable_sites = find_probable_sites(bgm_df, threshold)
    labeled_sites = add_within_between(probable_sites, protein_length)
    csv_name = f"{base_name1}_{base_name2}_{int(threshold * 100)}_probable_sites.csv"
    labeled_sites.to_csv(csv_name)

    density, directionality_for, directionality_rev = generate_summary(labeled_sites)
    summary_name = f"{base_name1}_{base_name2}_summary.txt"
    with open(summary_name, "w") as summary_file:
        # Write the nutshell information between target proteins.
        summary_file.write("Proportion of interactions between {} that are in-between: {}\n\n"
                           .format(f"{base_name1}_{base_name2}", density))
        summary_file.write("Percentage of directional interactions from {} to its protein partner (direction of "
                           "conditional dependence): {}%\n\n"
                           .format(base_name1, round(directionality_for, 5) * 100))
        summary_file.write("Percentage of directional interactions from protein partner to {} (direction of "
                           "conditional dependence): {}%\n\n"
                           .format(base_name1, round(directionality_rev, 5) * 100))

    if directory:
        if not os.path.exists(f"{directory}/{base_name1}_{base_name2}"):
            make_dir = f"mkdir {directory}/{base_name1}_{base_name2}"
            subprocess.call(make_dir, shell=True)

        cmd_csv = f"mv {csv_name} {directory}/{base_name1}_{base_name2}"
        cmd_summary = f"mv {summary_name} {directory}/{base_name1}_{base_name2}"
        subprocess.call(cmd_csv, shell=True)
        subprocess.call(cmd_summary, shell=True)


def process_files(p1, p2, tree):
    # Output file names, change to "aligned_{os.path.basename(p1)}" if using the align_fasta function.
    global bgm_file_path
    aligned_protein1 = f"{os.path.basename(p1)}"
    aligned_protein2 = os.path.join(os.path.dirname(p2), f"{os.path.basename(p2)}")

    protein1_fasta = open(aligned_protein1, "r").readlines()
    protein2_fasta = open(aligned_protein2, "r").readlines()
    protein1_dict = make_seq_dict(protein1_fasta)
    protein2_dict = make_seq_dict(protein2_fasta)

    species_list = [i for i in protein1_dict.keys() if i in protein2_dict.keys()]

    concat_dict = {}

    for i in species_list:
        concat_dict[i] = protein1_dict[i] + protein2_dict[i]

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

    base_dir = f"{batch_directory}/" if batch else ""

    # Calculate BGM
    if batch:
        # Building full paths
        concat_fasta_path = os.path.join(concat_fasta)
        shared_phylo_tree_path = os.path.join(shared_phylo_tree)
        bgm_file_path = os.path.join(base_dir, f"{concat_fasta}.BGM.json")

        call_bgm(concat_fasta_path, shared_phylo_tree_path)
        cmd = f"mv {concat_fasta}.BGM.json {base_dir}"
        subprocess.call(cmd, shell=True)
    else:
        call_bgm(concat_fasta, shared_phylo_tree)

    bgm_df = initialize_results(bgm_file_path)

    # Process for both 0.9 and 0.5 threshold values
    for threshold in [0.9, 0.5]:
        write_summary_and_csv(bgm_df, protein1_length, threshold, base_name1, base_name2, base_dir)

    if make_graph:
        create_graph(bgm_df, 0.9)

        if batch:
            cmd_connected = f"mv {protein1}_{protein2}_connected_nodes.png {base_dir}/"
            cmd_all = f"mv {protein1}_{protein2}_all_nodes.png {base_dir}/"
            subprocess.call(cmd_connected, shell=True)
            subprocess.call(cmd_all, shell=True)


# Post BGM file processing is below:
def length_of_first_sequence(fasta_file):
    for record in SeqIO.parse(fasta_file, "fasta"):
        return len(record.seq)


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

    bgm_df = pd.DataFrame(
        columns=["Site1", "Site2", "Site1 --> Site2", "Site1 <-- Site2", "Site1 <--> Site2", "Site1 Subs", "Site2 Subs",
                 "Shared Subs"])

    bgm_df["Site1"] = site1_list
    bgm_df["Site2"] = site2_list
    bgm_df["Site1 --> Site2"] = site1_to_site2_list
    bgm_df["Site1 <-- Site2"] = site2_to_site1_list
    bgm_df["Site1 <--> Site2"] = site1_or_site2_list
    bgm_df["Site1 Subs"] = site1_subs_list
    bgm_df["Site2 Subs"] = site2_subs_list
    bgm_df["Shared Subs"] = shared_subs_list
    return bgm_df


def find_probable_sites(bgm_df, confidence=0.9):
    probable_sites_df = bgm_df[bgm_df["Site1 <--> Site2"] >= confidence]
    return probable_sites_df


# Creates graphs displaying bgm data.
def create_graph(bgm_df, confidence):
    prob_sites = find_probable_sites(bgm_df, confidence)
    g = Graph()
    g2 = Graph()
    nodes_added = []
    for i in range(0, len(prob_sites)):
        g.add_edge(int(prob_sites.iloc[i]["Site1"]), int(prob_sites.iloc[i]["Site2"]))
        g2.add_edge(int(prob_sites.iloc[i]["Site1"]), int(prob_sites.iloc[i]["Site2"]))
        if not (prob_sites.iloc[i]["Site1"] in nodes_added):
            nodes_added += [prob_sites.iloc[i]["Site1"]]
        if not (prob_sites.iloc[i]["Site2"] in nodes_added):
            nodes_added += [prob_sites.iloc[i]["Site2"]]
    for i in range(1, int(bgm_df["Site1"][bgm_df.index.values[-1]]) + 1):
        if not (i in nodes_added):
            g.add_node(i)

    nodes_with_edges = [node for node, degree in g.degree() if degree > 0]
    nodes_with_multi_edges = [node for node, degree in g.degree() if degree > 1]
    node_colors = ["red" if node in nodes_with_edges else "lightblue" if node not in nodes_with_multi_edges else "red"
                   for node in g.nodes()]

    pos = spring_layout(g, k=.3, scale=.1, center=(1, 1), iterations=20)
    pos_nodes_with_edges = {node: position for (node, position) in pos.items() if node in nodes_with_edges}

    draw_networkx(g2, pos=pos_nodes_with_edges, node_color="red", node_size=150, font_size=6)
    plt.savefig(f"{protein1}_{protein2}_connected_nodes.png")

    draw_networkx(g, pos=pos, node_color=node_colors, node_size=150, font_size=6)
    plt.savefig(f"{protein1}_{protein2}_all_nodes.png")


# Determines if interactions are within one protein or between the two proteins, as well as determines direction of
# interaction.
def add_within_between(bgm_data, protein1_length):
    # Create a copy of the dataframe to avoid modifying the original dataframe/view
    bgm_data_copy = bgm_data.copy()

    interaction_types = []
    direction_types = []

    for entry_ID in bgm_data_copy.index:
        print()
        if int(bgm_data_copy["Site1"][entry_ID]) <= protein1_length < int(bgm_data_copy["Site2"][entry_ID]):
            interaction_types.append("Between")
        else:
            interaction_types.append("Within")

        if (float(bgm_data_copy["Site1 --> Site2"][entry_ID]) / float(
                bgm_data_copy["Site1 <--> Site2"][entry_ID])) > 0.8:
            direction_types.append("1->2")
        elif (float(bgm_data_copy["Site1 <-- Site2"][entry_ID]) / float(
                bgm_data_copy["Site1 <--> Site2"][entry_ID])) > 0.8:
            direction_types.append("2->1")
        else:
            direction_types.append("Bi")

    bgm_data_copy["Identity"] = interaction_types
    bgm_data_copy["Direction"] = direction_types

    return bgm_data_copy


# Generates summary data: Density and Directionality
def generate_summary(bgm_data):
    density = (len(bgm_data[bgm_data["Identity"] == "Between"]) / protein1_length)
    directionality_forward = (len(bgm_data[bgm_data["Direction"] == "1->2"]) / len(bgm_data))
    directionality_reverse = (len(bgm_data[bgm_data["Direction"] == "2->1"]) / len(bgm_data))
    return density, directionality_forward, directionality_reverse


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
        for filename in os.listdir(protein2):
            protein2_files = glob.glob(os.path.join(protein2, '*.fa*'))
            for protein2 in protein2_files:
                # Make a deep copy of the master tree for each protein2
                process_files(protein1, protein2, master_tree)
    else:
        process_files(protein1, protein2, master_tree)
