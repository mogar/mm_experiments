"""
Carbon Bond Creator for SAMSON
Adds bonds between carbon atoms within a specified distance threshold
"""

from samson import *
import math

import utilities

class CarbonBondCreator:
    """
    Creates bonds between carbon atoms that are within a specified distance.

    Parameters:
    - distance_threshold: Maximum distance (in angstroms) for bond creation
    - bond_order: Bond order for created bonds (default 1.0)
    """

    def __init__(self, distance_threshold=3.0, bond_order=1.0):
        self.distance_threshold = distance_threshold
        self.bond_order = bond_order
        self.bonds_created = 0

    def get_selected_carbon_atoms(self):
        """Get all selected carbon atoms from the active document"""
        # Get all selected atoms
        selected_atoms = SAMSON.getNodes('node.type atom and node.selected')

        # Filter for carbon atoms only
        carbon_atoms = []
        for atom in selected_atoms:
            if atom.elementType == SBElement.Carbon:
                carbon_atoms.append(atom)

        return carbon_atoms

    def create_bonds(self):
        """
        Main method to create bonds between nearby carbon atoms.
        Returns the number of bonds created.
        """
        print(f"Starting bond creation with distance threshold: {self.distance_threshold} Å")

        # Get selected carbon atoms
        carbon_atoms = self.get_selected_carbon_atoms()
        num_atoms = len(carbon_atoms)

        print(f"Found {num_atoms} selected carbon atoms")

        if num_atoms == 0:
            print("No carbon atoms selected. Please select atoms first.")
            return 0

        # Get the structural model (parent of atoms)
        structural_model = None
        if num_atoms > 0:
            # Navigate up to find the structural model
            parent = carbon_atoms[0].getParent()
            while parent is not None:
                if isinstance(parent, SBStructuralModel):
                    structural_model = parent
                    break
                parent = parent.getParent()

        if structural_model is None:
            print("Warning: Could not find structural model. Using active document.")
            structural_model = SAMSON.getActiveDocument()

        # Begin undoable operation
        SAMSON.beginHolding("Create carbon bonds")

        self.bonds_created = 0
        bonds_to_create = []

        # Iterate through all pairs of carbon atoms
        for i in range(num_atoms):
            atom1 = carbon_atoms[i]

            for j in range(i + 1, num_atoms):
                atom2 = carbon_atoms[j]

                # Check if already bonded
                if utilities.atoms_already_bonded(atom1, atom2):
                    continue

                # Calculate distance
                distance = utilities.calculate_distance(atom1, atom2)

                # Create bond if within threshold
                if distance <= self.distance_threshold:
                    bonds_to_create.append((atom1, atom2, distance))

        # Create all the bonds
        for atom1, atom2, distance in bonds_to_create:
            bond = SBBond(atom1, atom2, self.bond_order)
            SAMSON.hold(bond)
            bond.create()
            structural_model.addChild(bond)
            self.bonds_created += 1

            # Print progress for every 100 bonds
            if self.bonds_created % 100 == 0:
                print(f"Created {self.bonds_created} bonds...")

        # End undoable operation
        SAMSON.endHolding()

        print(f"\nBond creation complete!")
        print(f"Total bonds created: {self.bonds_created}")
        print(f"Checked {num_atoms * (num_atoms - 1) // 2} atom pairs")

        return self.bonds_created


# Example usage
if __name__ == "__main__":
    distance_string = input("Distance threshold for bonding (Angstroms, default 1.6): ")
    distance_threshold = float(distance_string) if distance_string else 1.6

    bond_creator = CarbonBondCreator(distance_threshold=distance_threshold)
    bond_creator.create_bonds()
