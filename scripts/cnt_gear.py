"""
Carbon Nanotube Gear Generator for SAMSON
Generates armchair (n,n) nanotubes with lonsdaleite teeth for molecular gears
"""

import math
from samson import *

import armchair_cnt
import lonsdaleite_tooth
import utilities

class CNTGearGenerator:
    """
    Generates carbon nanotube-based molecular gears with lonsdaleite teeth.

    Parameters:
    - n: Chiral index for armchair nanotube (M=0 means armchair: (n,n))
    - z_length: Length of nanotube in angstroms
    - num_teeth: Number of gear teeth around the circumference
    """

    def __init__(self, n=6, z_length=100.0, num_teeth=6, tooth_height=10):
        self.n = n
        self.z_length = z_length
        self.num_teeth = num_teeth
        self.tooth_height = tooth_height

        # Carbon-carbon bond length in graphene/nanotubes (angstroms)
        self.acc = 1.413

    def bond_tooth_to_cnt(self, tooth_atoms, tooth_grid, tooth_grid_structure, cnt_atoms, structural_model):
        """Bond tooth base atoms to nearest CNT surface atoms"""
        SAMSON.beginHolding("Bond tooth to CNT")

        (num_x_cells, num_y_cells, num_z_cells) = tooth_grid_structure

        # only seek to bond the first layer
        for z_idx in range(num_z_cells + 1):
            for ring_idx in range(2 * num_x_cells + 1):
                tooth_atom = tooth_grid.get((ring_idx, 0, z_idx))
                if tooth_atom is None:
                    continue

                min_dist = float('inf')
                closest_cnt_atom = None

                for cnt_atom in cnt_atoms:
                    if cnt_atom is None:
                        continue
                    dist = utilities.calculate_distance(tooth_atom, cnt_atom)
                    if dist < min_dist:
                        min_dist = dist
                        closest_cnt_atom = cnt_atom

                # Create bond if reasonably close
                if closest_cnt_atom and min_dist < (self.acc * 1.1):
                    bond = SBBond(tooth_atom, closest_cnt_atom, 1.0)
                    SAMSON.hold(bond)
                    bond.create()
                    structural_model.addChild(bond)

        SAMSON.endHolding()

    def attach_teeth_to_cnt(self, structural_model, cnt_atoms, radius, tooth_height):
        """Attach lonsdaleite teeth to the nanotube at specified positions"""

        SAMSON.beginHolding("Add gear teeth")
        translation_vector = samson.SBPhysicalVector3(SBQuantity.angstrom(0),
                                                        SBQuantity.angstrom(radius + self.acc),
                                                        SBQuantity.angstrom(0))

        teeth = []
        for tooth_idx in range(self.num_teeth):
            angle_rad = 2 * math.pi * tooth_idx / self.num_teeth
            c, s = math.cos(angle_rad), math.sin(angle_rad)

            tooth_width = 6 # TODO: calculate?
            ingot = lonsdaleite_tooth.LonsdaleiteTooth(tooth_width, tooth_height, self.z_length, taper_x=True, taper_z=False)
            ingot_model, ingot_atoms, ingot_grid, ingot_cell_structure = ingot.generate_tooth()
            # remove first and last column of tool atoms
            (num_x_cells, num_y_cells, num_z_cells) = ingot_cell_structure
            for z_idx in range(num_z_cells + 1):
                left_atom = ingot_grid.pop((0, num_y_cells - 1, z_idx), None)
                if left_atom in ingot_atoms:
                    ingot_atoms.remove(left_atom)
                    utilities.remove_atom_and_bonds(left_atom)
                right_atom = ingot_grid.pop((2 * num_x_cells, num_y_cells - 1, z_idx), None)
                if right_atom in ingot_atoms:
                    ingot_atoms.remove(right_atom)
                    utilities.remove_atom_and_bonds(right_atom)
            # translate/rotate molecule atom by atom
            for atom in ingot_atoms:
                position = atom.getPosition()
                position += translation_vector
                x, y, z = position[0], position[1], position[2]
                new_position = SBPosition3(c*x - s*y, s*x + c*y, z)
                atom.setPosition(new_position)

            # Find and bond closest CNT atoms to tooth base
            self.bond_tooth_to_cnt(ingot_atoms, ingot_grid, ingot_cell_structure, cnt_atoms, structural_model)
            teeth.append(ingot_model)

        SAMSON.endHolding()
        return teeth

    def generate_gear(self):
        """Main method to generate complete CNT gear"""
        print(f"Generating CNT gear: n={self.n}, length={self.z_length}Å, teeth={self.num_teeth}")

        if self.n % self.num_teeth != 0:
            print("number of teeth is not a factor of CNT parameter N")
            print("Tooth alignment may not be optimal")

        SAMSON.beginHolding("Create Structural Model")
        structural_model = SBStructuralModel()
        structural_model.name = f"gear"
        structural_model.create()
        SAMSON.endHolding()

        # Generate base CNT structure
        cnt = armchair_cnt.CNTGenerator(self.n, self.z_length)
        cnt_molecule, cnt_atoms, radius = cnt.generate_armchair_cnt()

        SAMSON.beginHolding("Add CNT to model")
        structural_model.addChild(cnt_molecule)
        SAMSON.endHolding()

        # Attach teeth
        teeth = self.attach_teeth_to_cnt(structural_model, cnt_atoms, radius, self.tooth_height)

        SAMSON.beginHolding("Add teeth to model")
        for tooth in teeth:
            structural_model.addChild(tooth)
        SAMSON.endHolding()

        utilities.remove_dangling_carbons(structural_model)

        # Add to document
        document = SAMSON.getActiveDocument()
        SAMSON.beginHolding(f"Create CNT Gear n={self.n}")
        SAMSON.hold(structural_model)
        document.addChild(structural_model)
        SAMSON.endHolding()

        print("CNT gear generation complete!")
        return structural_model


# Example usage
if __name__ == "__main__":
    n = int(input("Gearshaft CNT N: "))
    z_length = float(input("Gear length (Angstroms): "))
    num_teeth = int(input("num teeth: "))
    tooth_height = float(input("Tooth height (Angstroms): "))

    generator = CNTGearGenerator(n, z_length, num_teeth, tooth_height)
    gear = generator.generate_gear()
