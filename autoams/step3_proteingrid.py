import MDAnalysis as mda
import numpy as np
import random
from Bio.PDB import PDBParser, PDBIO, Model, Structure
__all__ = ["PeptideGrid"]

class PeptideGrid():
    '''docstring!'''
    def __init__(self, pdb_file, distance=0, num_peptides=16, minim_sep=20):
        self.molecule = mda.Universe(str(pdb_file))
        self.normal_vector = np.array([0.0, 0.0, 1.0])
        self.distance_from_origin = distance
        self.replicates = num_peptides
        self.minimum_separation = minim_sep
        self.molecules_in_plane = []

    def generate_grid(self):
        ''' generates the grid of peptides in the horizontal axis '''

        for i in range(self.replicates):
            valid_position = False

            while not valid_position:
                # Generate random angles to determine the direction
                theta = random.uniform(0, 2 * np.pi)  # Random angle in radians

                # Calculate the offsets in X and Y directions
                x_offset = random.uniform(20.0, 100.0) * np.cos(theta)  # Random separation, at least 20 nm
                y_offset = random.uniform(20.0, 100.0) * np.sin(theta)  # Random separation, at least 20 nm

                # Calculate the corresponding Z coordinate to ensure they are in the same plane
                z_offset = -(self.normal_vector[0] * x_offset + self.normal_vector[
                    1] * y_offset + self.distance_from_origin) / self.normal_vector[2]

                # Create a copy of the original molecule and position it
                new_molecule = self.molecule.copy()
                new_molecule.atoms.positions += [x_offset, y_offset, z_offset]

                # Check if this molecule's position is at least 20 nm from every other molecule
                min_distance = float('inf')

                for mol in self.molecules_in_plane:
                    distances = np.linalg.norm(new_molecule.atoms.positions - mol.atoms.positions, axis=1)
                    min_distance = min(min_distance, np.min(distances))

                if min_distance >= self.minimum_separation:
                    self.molecules_in_plane.append(new_molecule)
                    valid_position = True

    def save_grid(self):
        '''Saves the grid of molecules in the horizontal plane'''
        # Create a combined Universe
        combined_universe = mda.Merge(*[mol.atoms for mol in self.molecules_in_plane])

        # Save the combined Universe as a single PDB file
        combined_universe.atoms.write('all_molecules_combined.pdb')
        print("Saved all molecules as all_molecules_combined.pdb")
