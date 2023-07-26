import pandas as pd
import json
from networkx import *
import matplotlib.pyplot as plt
import sys

help_text = "Arguments:\n--graph True/False\n\tDefault False, determines if a graph will be drawn.\
            \n\n--threshold num\n\t Default 0.9, confidence of conditional dependence.\
            \n\n--file filename\n\t Required, filename of json file containing bgm data.\
            \n\n--out filename\n\t Root name for output files."

help_flag = False
make_graph = False
threshold = 0.9
filename = ""
out = "bgm_output"


# Determines what to do based on command line input.
for i in range(len(sys.argv)):
    if sys.argv[i][:2] == "--":
        if sys.argv[i] == "--graph":
            if not(sys.argv[i+1] == "True" or sys.argv[i+1] == "False"):
                raise ValueError("True or False are the only accepted input for the '--graph' flag.")
            else:
                if sys.argv[i+1] == "True":
                    make_graph = True
        elif sys.argv[i] == "--threshold":
            if 1 > float(sys.argv[i + 1]) > 0:
                threshold = float(sys.argv[i+1])
            else:
                raise ValueError("Input after '--threshold' flag must be between 0 and 1 non-inclusive.")
        elif sys.argv[i] == "--file":
            filename = sys.argv[i+1]
        elif sys.argv[i] == "--help":
            help_flag = True
            print(help_text)
        elif sys.argv[i] == "--out":
            out = sys.argv[i+1]
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


# Creates graphs displaying bgm data.
def create_graph(bgm_df, confidence=0.9):
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
    for i in range(1, int(bgm_df["Site1"][bgm_df.index.values[-1]]) + 1):
        if not(i in nodes_added):
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

    probable_sites.to_csv(out + "_probable_sites.csv")

# Create and save the graphs.
if make_graph:
    # Ensure bgm_df is defined when needed
    if 'bgm_df' not in locals():
        bgm_df = initialize_results(filename)
    create_graph(bgm_df, threshold)
