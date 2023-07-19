from Bio import SeqIO

Protein1_fasta = open("59272_ACE2.fa", "r").readlines()

Protein2_fasta = open("84873_ADGRG7.fa", "r").readlines()


def make_seq_dict(fasta_lines):
    seq_dict = {}
    current_seq = ""
    for i in range(len(fasta_lines)):
        if fasta_lines[i][0] == ">":
            current_seq = fasta_lines[i][1:-1]
            seq_dict[current_seq] = ""
        else:
            seq_dict[current_seq] += fasta_lines[i][:-1]
    return(seq_dict)
            

Protein1_dict = make_seq_dict(Protein1_fasta)
Protein2_dict = make_seq_dict(Protein2_fasta)

species_list = [i for i in Protein1_dict.keys() if i in Protein2_dict.keys()]

concat_dict = {}

for i in species_list:
    concat_dict[i] = Protein1_dict[i] + Protein2_dict[i]

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
    return(final_string)


out = ""
taxa_list = ""
for i in concat_dict.keys():
    out += ">"
    out += i
    out += "\n"
    out += split_by_60(concat_dict[i])
    taxa_list += i.replace("_", " ")
    taxa_list += "\n"
	

with open("ACE2_ADGRG7.fa", "w") as f:
    f.write(out)

with open("ACE2_ADGRG7_taxa_list.csv", "w") as l:
    l.write(taxa_list)
    
    
