"""
Lonsdaleite gear rack generation script for SAMSON.
"""

import math
import operator
from samson import *

import guiderail
import lonsdaleite_ingot
import utilities

class GearRackGenerator:
    """Class to generate a gear rack from lonsdaleite ingot."""

    def __init__(self, length, num_clips, tooth_width, tooth_height):
        self.length = length
        self.num_clips = num_clips
        self.tooth_width = tooth_width
        self.tooth_height = tooth_height

        # Lonsdaleite crystal dimensions
        self.lons_y_step = 2.16  # angstroms between rows in y direction
        self.lons_z_step = 2.1373  # angstroms between layers in z direction
        self.lons_x_step = 2.513/2 # angstroms between atoms in x direction

        # For now, dimensions of guiderail are hardcoded
        self.anchor_width = 20
        self.anchor_height = 8
        self.narrow_width = 10
        self.narrow_height = 16
        self.rail_width = 20
        self.rail_height = 10
        self.steps_per_tooth = int(self.tooth_width / self.lons_z_step) * 2 + 3 # + 3 to ensure enough space for hydrogens
        self.num_teeth = round(self.length / (self.steps_per_tooth * self.lons_z_step))

    def bond_tooth_to_anchor(self, tooth_grid, tooth_cell_structure, anchor_grid, anchor_cell_structure, anchor_offset, structural_model):
        """Bond tooth base atoms to nearest anchor surface atoms"""
        SAMSON.beginHolding("Bond tooth to anchor")

        (tooth_x_cells, tooth_y_cells, tooth_z_cells) = tooth_cell_structure
        (anchor_x_cells, anchor_y_cells, anchor_z_cells) = anchor_cell_structure

        # only seek to bond the 0 layer
        for z_idx in range(tooth_cell_structure[2] + 1):
            for ring_idx in range(2 * tooth_cell_structure[0] + 1):
                tooth_face_atom = tooth_grid.get((ring_idx, 0, z_idx))
                if tooth_face_atom is None:
                    continue

                # The tooth was rotated 180 degrees to face the anchor
                anchor_ring_idx = 2 * tooth_cell_structure[0] - ring_idx
                anchor_face_atom = anchor_grid.get((anchor_ring_idx + anchor_offset[0], anchor_offset[1], z_idx + anchor_offset[2]))
                if anchor_face_atom is None:
                    continue

                dist = utilities.calculate_distance(tooth_face_atom, anchor_face_atom)
                if dist < (self.lons_z_step * 1.1):
                    min_dist = dist
                else:
                    continue

                bond = SBBond(tooth_face_atom, anchor_face_atom, 1.0)
                SAMSON.hold(bond)
                bond.create()
                structural_model.addChild(bond)


        SAMSON.endHolding()

    def generate(self):
        """Generate the gear rack structure."""
        # Generate the guiderail first
        rail_generator = guiderail.GuideRailGenerator(
            anchor_width=self.anchor_width,
            anchor_height=self.anchor_height,
            rail_width=self.rail_width,
            rail_height=self.rail_height,
            narrow_width=self.narrow_width,
            narrow_height=self.narrow_height,
            length=self.length,
        )
        
        rail, rail_anchor_grid, rail_anchor_cell_structure = rail_generator.generate_rail()

        SAMSON.beginHolding("Add gear teeth")

        angle_rad = math.pi # flip them around to face the anchor
        c, s = math.cos(angle_rad), math.sin(angle_rad)

        # Generate the gear teeth clips
        print(f"Generating {self.num_teeth} rack teeth...")
        for tooth_idx in range(self.num_teeth):
            ingot = lonsdaleite_ingot.LonsdaleiteIngot(self.anchor_width, tooth_height, tooth_width)
            ingot_model, ingot_atoms, ingot_grid, ingot_cell_structure = ingot.generate_ingot()
            # remove first and last column of tool atoms
            (num_x_cells, num_y_cells, num_z_cells) = ingot_cell_structure
            for x_idx in range(2*num_x_cells + 1):
                left_atom = ingot_grid.pop((x_idx, num_y_cells - 1, 0), None)
                if left_atom in ingot_atoms:
                    ingot_atoms.remove(left_atom)
                    utilities.remove_atom_and_bonds(left_atom)
                right_atom = ingot_grid.pop((2 * num_x_cells, num_y_cells - 1, num_z_cells - 1), None)
                if right_atom in ingot_atoms:
                    ingot_atoms.remove(right_atom)
                    utilities.remove_atom_and_bonds(right_atom)
                if x_idx % 2 == 0:
                    left_atom_neighbor = ingot_grid.pop((x_idx, num_y_cells - 2, 0), None)
                    if left_atom_neighbor in ingot_atoms:
                        ingot_atoms.remove(left_atom_neighbor)
                        utilities.remove_atom_and_bonds(left_atom_neighbor)
                if x_idx % 2 == 1:
                    right_atom_neighbor = ingot_grid.pop((x_idx, num_y_cells - 1, num_z_cells), None)
                    if right_atom_neighbor in ingot_atoms:
                        ingot_atoms.remove(right_atom_neighbor)
                        utilities.remove_atom_and_bonds(right_atom_neighbor)
            
            top_y = utilities.get_y_of_lonsdaleite(ingot_grid, ingot_cell_structure, comparison=operator.gt)
            # translate/rotate molecule atom by atom
            z_translation_steps = tooth_idx * self.steps_per_tooth
            x_offset = 0
            if z_translation_steps % 2 == 1:
                x_offset = 1
            translation_vector = SBPhysicalVector3(SBQuantity.angstrom(x_offset * self.lons_x_step),
                                         -SBQuantity.angstrom(self.lons_y_step),
                                         SBQuantity.angstrom(z_translation_steps * self.lons_z_step))
            for atom in ingot_atoms:
                position = atom.getPosition()
                x, y, z = position[0], position[1], position[2]
                x, y, z = c*x - s*y, s*x + c*y, z
                new_position = SBPosition3(x, y, z)
                new_position += translation_vector
                atom.setPosition(new_position)

            # Find and bond closest CNT atoms to tooth base
            anchor_offset = (x_offset, 0, z_translation_steps)
            self.bond_tooth_to_anchor(ingot_grid, ingot_cell_structure, rail_anchor_grid, rail_anchor_cell_structure, anchor_offset, rail)

            rail.addChild(ingot_model)

        SAMSON.endHolding()

        clips = []
        for clip_idx in range(self.num_clips):
            clip = rail_generator.generate_clip()
            clip.translate(SBPhysicalVector3(SBQuantity.angstrom(0),
                                     SBQuantity.angstrom(0),
                                     SBQuantity.angstrom(clip_idx * 40)))
            clips.append(clip)
        return rail, clips

        

# Example usage
if __name__ == "__main__":
    length = float(input("Guiderail length (Angstroms): "))
    num_clips = int(input("Number of clips: "))
    tooth_height = float(input("Tooth height (Angstroms): "))
    tooth_width = float(input("Tooth width (Angstroms): "))


    generator = GearRackGenerator(length=length,
                                 num_clips=num_clips,
                                 tooth_width=tooth_width,
                                 tooth_height=tooth_height)
    (rail, clips) = generator.generate()
