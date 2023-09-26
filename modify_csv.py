import sys
import pandas as pd
from Bio.PDB import PDBParser
import numpy

def parse_pdb(pdb_file):
    """Parse a PDB file and return a Biopython PDB object."""
    parser = PDBParser()
    structure = parser.get_structure('model', pdb_file)
    return structure

def calculate_c_alpha_distance(structure, modified_site1, modified_site2):
    """Calculate the C-alpha distance between two residues in a Biopython PDB structure."""
    if modified_site1 != '-' and modified_site2 != '-':
        modified_site1 = int(modified_site1)
        modified_site2 = int(modified_site2)

    try:
        residue1 = structure[0]['A'][modified_site1]
        residue2 = structure[0]['A'][modified_site2]

        c_alpha1 = residue1['CA']
        c_alpha2 = residue2['CA']

        distance = numpy.linalg.norm(c_alpha1.coord - c_alpha2.coord)
        return distance
    except KeyError:
        return None

def add_molecular_distances(df, structure):
    """Add a new column to a DataFrame called `molecular_distances` and populate it with the calculated C-alpha distances."""

    df['molecular_distances'] = df[['modified_site1', 'modified_site2']].apply(lambda x: calculate_c_alpha_distance(structure, x[0], x[1]), axis=1)

    # Convert `molecular_distances` to numerics and handle None values.
    df['molecular_distances'] = pd.to_numeric(df['molecular_distances'], errors='coerce').fillna('-')


def parse_fasta(fasta_file):
    """Parse a FASTA file and return a dictionary of sequences."""
    sequences = {}
    with open(fasta_file, 'r') as f:
        current_key = None
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                current_key = line[1:]
                sequences[current_key] = ""
            else:
                sequences[current_key] += line
    return sequences

def calculate_modified_site(site, sequence):
    """Calculate the modified site value based on the number of gaps."""

    # Remove all the gaps in human sequence and then count total number of characters.
    homo_sapiens_seq = sequence.replace('-', '')
    total_number_of_characters = len(homo_sapiens_seq)
    print("Total number of characters:{}".format(total_number_of_characters))

    # If the modified columns site number greater than that value, then record it as '-'.
    if site > total_number_of_characters:
        return '-'

    # Return the modified site value.
    return site

def check_original_site(modified_site, sequence):
    """Check the original site #. if that human sequence residue is '-', change that location's residue back to '-'."""
    if modified_site == '-':
        return '-'

    modified_site = int(modified_site)
    if sequence[modified_site - 1] == '-':
        return '-'

    return modified_site

def generate_chimera_x_command(df):
    # Filter out rows where either modified_site1 or modified_site2 is '-'.
    filtered_df = df[(df['modified_site1'] != '-') & (df['modified_site2'] != '-')]

    # Convert the remaining modified_site1 and modified_site2 values to a list.
    sites_list = filtered_df['modified_site1'].tolist() + filtered_df['modified_site2'].tolist()

    # Construct the Chimera X command using the list of sites.
    sites_str = ', '.join(map(str, sites_list))
    command = f"select /A: {sites_str}"

    return command

def main(csv_file, fasta_file, pdb_file):
    # Load the CSV
    df = pd.read_csv(csv_file)

    # Parse the FASTA file
    sequences = parse_fasta(fasta_file)
    homo_sapiens_seq = sequences.get('Homo_sapiens', '')

    # Calculate modified_site1 and modified_site2
    df['modified_site1'] = df['Site1'].apply(lambda x: calculate_modified_site(x, homo_sapiens_seq))
    df['modified_site2'] = df['Site2'].apply(lambda x: calculate_modified_site(x, homo_sapiens_seq))

    # Check the original site #. if that human sequence residue is '-', change that location's residue back to '-'.
    df['modified_site1'] = df['modified_site1'].apply(lambda x: check_original_site(x, homo_sapiens_seq))
    df['modified_site2'] = df['modified_site2'].apply(lambda x: check_original_site(x, homo_sapiens_seq))

    # Generate Chimera X command
    command = generate_chimera_x_command(df)
    print(command)

    # Parse the PDB file
    structure = parse_pdb(pdb_file)

    # Add the molecular_distances column
    add_molecular_distances(df, structure)

    # Save the modified CSV
    df.to_csv(csv_file, index=False)
    print(f"Modified CSV saved to {csv_file}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 modify_csv.py <path_to_csv> <path_to_fasta> <path_to_pdb>")
        sys.exit(1)
    csv_file = sys.argv[1]
    fasta_file = sys.argv[2]
    pdb_file = sys.argv[3]
    main(csv_file, fasta_file, pdb_file)
