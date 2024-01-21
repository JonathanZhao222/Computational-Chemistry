from Bio.PDB import PDBParser

# Replace 'your_pdb_file.pdb' with your actual PDB file
parser = PDBParser()
structure = parser.get_structure("my_structure", "pdb/1gyc.pdb")

metal_ions = []

for model in structure:
    for chain in model:
        for residue in chain:
            for atom in residue:
                # Check for atoms typically associated with metals
                if atom.element == "CA" or atom.element == "FE" or atom.element == "MG" or atom.element == "CU":
                    metal_ions.append(residue)
                    break  # If any atom in the residue is a metal, add the whole residue and break

print("Metal ions found:", len(metal_ions))
for ion in metal_ions:
    print(f"Metal ion: {ion.get_full_id()}")  # Printing the metal ions found



