import pandas as pd
import json
from networkx import *
import matplotlib.pyplot as plt
import sys
from Bio import SeqIO

'''
A script to process json file(s) generated by BGM tests.

Required arguments:
python3 format_bgm_json_result.py --file [filename] --out [out_filename] --protein1 [protein1.fa]

'--file [filename]' is required as the file name of json file containing BGM data.
'--out [out_filename]' is required as the root-names for output files.
'--protein1 [protein1.fa]' is required as the fasta file for the first protein in the concatenated file.

Optional arguments:
Use '--graph [True/False]' to determine if a graph will be drawn, default is false.
Use '--threshold [number from [0,1]]' to determine the level of filtering interactions with a certain level of
conditional dependence, default is 0.9.
'''

help_text = "Arguments:\n--graph True/False\n\t Default False, determines if a graph will be drawn.\
            \n\n--threshold num\n\t Default 0.9, confidence of conditional dependence.\
            \n\n--file filename\n\t Required, filename of json file containing bgm data.\
            \n\n--out filename\n\t Root name for output files.\
            \n\n--protein1 protein1.fa\n\t Required, fasta file for the first protein in the concatenated sequence."

help_flag = False
make_graph = False
threshold = 0.9
filename = ""
out = "bgm_output"
protein1_length = 0


# Returns length of the first sequence in a fasta file.
def length_of_first_sequence(fasta_file):
    for record in SeqIO.parse(fasta_file, "fasta"):
        return len(record.seq)


# Determines what to do based on command line input.
for i in range(len(sys.argv)):
    if sys.argv[i][:2] == "--":
        if sys.argv[i] == "--graph":
            if not (sys.argv[i + 1] == "True" or sys.argv[i + 1] == "False"):
                raise ValueError("True or False are the only accepted input for the '--graph' flag.")
            else:
                if sys.argv[i + 1] == "True":
                    make_graph = True
        elif sys.argv[i] == "--threshold":
            if 1 > float(sys.argv[i + 1]) > 0:
                threshold = float(sys.argv[i + 1])
            else:
                raise ValueError("Input after '--threshold' flag must be between 0 and 1 non-inclusive.")
        elif sys.argv[i] == "--file":
            filename = sys.argv[i + 1]
        elif sys.argv[i] == "--help":
            help_flag = True
            print(help_text)
        elif sys.argv[i] == "--out":
            out = sys.argv[i + 1]
        elif sys.argv[i] == "--protein1":
            base_protein = sys.argv[i+1]
            protein1_length = length_of_first_sequence(base_protein)
        else:
            raise ValueError("Unknown argument: {}".format(sys.argv[i]))


# Loads bgm result json file and creates a new dataframe containing only the results of the bgm.
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

    bgm_df = pd.DataFrame(columns=["Site1", "Site2", "Site1 --> Site2",
                                   "Site1 <-- Site2", "Site1 <--> Site2",
                                   "Site1 Subs", "Site2 Subs", "Shared Subs"])

    bgm_df["Site1"] = site1_list
    bgm_df["Site2"] = site2_list
    bgm_df["Site1 --> Site2"] = site1_to_site2_list
    bgm_df["Site1 <-- Site2"] = site2_to_site1_list
    bgm_df["Site1 <--> Site2"] = site1_or_site2_list
    bgm_df["Site1 Subs"] = site1_subs_list
    bgm_df["Site2 Subs"] = site2_subs_list
    bgm_df["Shared Subs"] = shared_subs_list

    return bgm_df


# Finds all sites in the dataframe above some confidence of conditional dependence.
def find_probable_sites(bgm_df, confidence=0.9, print_format=False):
    probable_sites_df = bgm_df[bgm_df["Site1 <--> Site2"] >= confidence]
    if print_format:
        return probable_sites_df[["Site1", "Site2", "Site1 <--> Site2", "Shared Subs"]]
    return probable_sites_df


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
    print(density)
    directionality_forward = (len(bgm_data[bgm_data["Direction"] == "1->2"]) / len(bgm_data))
    directionality_reverse = (len(bgm_data[bgm_data["Direction"] == "2->1"]) / len(bgm_data))
    return density, directionality_forward, directionality_reverse


# Creates graphs displaying bgm data.
def create_graph(bgm_df, confidence=0.9):
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
    node_colors = ["red" if node in nodes_with_multi_edges
                   else "red" if node in nodes_with_edges else "lightblue" for node in g.nodes()]

    pos = spring_layout(g, k=.3, scale=.1, center=(1, 1), iterations=20)
    pos_nodes_with_edges = {node: position for (node, position) in pos.items() if node in nodes_with_edges}

    draw_networkx(g2, pos=pos_nodes_with_edges, node_color="red", node_size=150, font_size=6)
    plt.savefig(out + "_connected_nodes.png", dpi=400)

    draw_networkx(g, pos=pos, node_color=node_colors, node_size=150, font_size=6)
    plt.savefig(out + "_all_nodes.png", dpi=400)


# Run the program.
if not help_flag:
    bgm_df = initialize_results(filename)
    probable_sites = find_probable_sites(bgm_df, threshold)
    labeled_sites = add_within_between(probable_sites, protein1_length)
    labeled_sites.to_csv(out + "_probable_sites.csv")
    density, directionality_for, directionality_rev = generate_summary(labeled_sites)
    summary_name = out + "_summary.txt"
    protein_name = base_protein.split('.fa')[0]  # This gets the string before '.fa'
    with open(summary_name, "w") as summary_file:
        # Write the nutshell information between target proteins.
        summary_file.write("Proportion of interactions between {} that are in-between: {}\n\n"
                           .format(out, density))
        summary_file.write("Percentage of directional interactions from {} to its protein partner (direction of "
                           "conditional dependence): {}%\n\n"
                           .format(protein_name, round(directionality_for, 5) * 100))
        summary_file.write("Percentage of directional interactions from protein partner to {} (direction of "
                           "conditional dependence): {}%\n\n"
                           .format(protein_name, round(directionality_rev, 5) * 100))

    summary_file.close()

# Create and save the graphs.
if make_graph:
    # Ensure bgm_df is defined when needed
    if 'bgm_df' not in locals():
        bgm_df = initialize_results(filename)
    create_graph(bgm_df, threshold)
