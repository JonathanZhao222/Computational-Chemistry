from Bio import PDB
import numpy as np
import matplotlib.pyplot as plt

def calculate_distance(atom1, atom2):
    """Calculate the Euclidean distance between two atoms."""
    return np.linalg.norm(atom1 - atom2)

def compare_distances(pdb_file, reference_residue_number, target_residue_numbers):
    """Compare the distances of target residues to the reference residue in a PDB file."""
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    reference_atom = None
    distances = {}

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[1] == reference_residue_number:
                    reference_atom = residue['CA']  # CA atom is commonly used for distance comparisons

    if reference_atom is None:
        raise ValueError(f"The reference residue {reference_residue_number} was not found in the PDB file.")

    for target_residue_number in target_residue_numbers:
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.id[1] == target_residue_number:
                        target_atom = residue['CA']
                        distance = calculate_distance(reference_atom.coord, target_atom.coord)
                        distances[target_residue_number] = distance

    return distances

def min_distances(pdb_file, reference_residue_number):

    target_residue_numbers = []  # Change to the desired target residue numbers
    for i in range(639):
        target_residue_numbers.append(i+63)

    try:
        distances = compare_distances(pdb_file, reference_residue_number, target_residue_numbers)
        distances_list = []
        for target_residue_number, distance in distances.items():
            print(f"The distance between residue {reference_residue_number} and residue {target_residue_number} is {distance:.2f} Ångstroms.")
            distances_list.append(distance)
    except ValueError as e:
        print(e)

    target_residue_numbers = target_residue_numbers[:648]

    print(target_residue_numbers)
    print(distances_list)

    #NEED TO FIND 5 CLOSEST AMINO ACID RESIDUES
    
def polarity():
    pass