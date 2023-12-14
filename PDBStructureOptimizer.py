import requests
import prody
from Bio import PDB
import pdbfixer
from simtk.openmm.app import PDBFile
import pyrosetta
from pyrosetta.rosetta.core.pose import Pose
from pyrosetta.rosetta.core.chemical import ChemicalManager
from pyrosetta import create_score_function
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
from pyrosetta.rosetta.core.scoring import ScoreType

def download_pdb(pdb_id):
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)
    if response.status_code == 200:
        pdb_filename = f"{pdb_id}.pdb"
        with open(pdb_filename, 'w') as file:
            file.write(response.text)
        return pdb_filename
    else:
        raise ValueError("Could not download PDB file")

def check_and_fix_pdb(pdb_file):
    fixer = pdbfixer.PDBFixer(filename=pdb_file)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixed_pdb_file = pdb_file.replace('.pdb', '_fixed.pdb')
    with open(fixed_pdb_file, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    return fixed_pdb_file

def check_quality_of_structure(pdb_file):
    structure = prody.parsePDB(pdb_file)
    bfactors = structure.getBetas()
    if bfactors is not None and bfactors.mean() > 60:
        return False
    return True

def create_pdb_with_pyrosetta(sequence, output_pdb_filename):
    pyrosetta.init()

    aa_1_to_3 = {
        "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP",
        "C": "CYS", "E": "GLU", "Q": "GLN", "G": "GLY",
        "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
        "M": "MET", "F": "PHE", "P": "PRO", "S": "SER",
        "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL"
    }

    pose = Pose()
    residue_type_set = ChemicalManager.get_instance().residue_type_set("fa_standard")

    for aa in sequence:
        three_letter_code = aa_1_to_3.get(aa, "")
        if not three_letter_code:
            raise ValueError(f"Invalid amino acid '{aa}' in sequence.")
        res_type = residue_type_set.name_map(three_letter_code)
        residue = pyrosetta.rosetta.core.conformation.ResidueFactory.create_residue(res_type)
        pose.append_residue_by_bond(residue)

    scorefxn = create_score_function("ref2015")
    mm = MoveMap()
    mm.set_bb(True)
    mm.set_chi(True)
    min_mover = MinMover()
    min_mover.score_function(scorefxn)
    min_mover.movemap(mm)
    min_mover.min_type("lbfgs_armijo_nonmonotone")

    num_minimization_steps = 100
    for _ in range(num_minimization_steps):
        min_mover.apply(pose)

    pose.dump_pdb(output_pdb_filename)
    print(f"Designed peptide saved to {output_pdb_filename}")

def extract_sequence_from_pdb(pdb_file):
    parser = PDB.PDBParser()
    structure = parser.get_structure('X', pdb_file)
    sequence = ''
    for model in structure:
        for chain in model:
            for residue in chain:
                if PDB.is_aa(residue):
                    sequence += PDB.Polypeptide.three_to_one(residue.get_resname())
    return sequence

# Main execution
pdb_id = '1m2z'  # Replace with the actual PDB ID
try:
    pdb_file = download_pdb(pdb_id)
    fixed_pdb = check_and_fix_pdb(pdb_file)
    if check_quality_of_structure(fixed_pdb):
        print(f"The structure in {fixed_pdb} is satisfactory.")
    else:
        print("Quality of structure is not satisfactory. Generating a new structure using PyRosetta.")
        sequence = extract_sequence_from_pdb(fixed_pdb)
        output_pdb_filename = 'regenerated_peptide.pdb'
        create_pdb_with_pyrosetta(sequence, output_pdb_filename)
except Exception as e:
    print(f"Error: {e}")

