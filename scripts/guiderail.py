"""
Lonsdaleite Molecular Guiderail Generator for SAMSON
"""

import math
from samson import *

import lonsdaleite_ingot

class GuideRailGenerator:
    """
    Generates lonsdaleite molecular guiderails.

    Parameters:
    - anchor_width: Width of the anchor ingot in angstroms
    - anchor_height: Height of the anchor ingot in angstroms
    - rail_width: Width of the rail ingot in angstroms
    - rail_height: Height of the rail ingot in angstroms
    - narrow_width: Width of the narrow ingot in angstroms
    - narrow_height: Height of the narrow ingot in angstroms
    - length: Length of the guiderail in angstroms
    """

    def __init__(self, anchor_width, anchor_height, rail_width, rail_height, narrow_width, narrow_height, length):
        self.anchor_width = anchor_width
        self.anchor_height = anchor_height
        self.rail_width = rail_width
        self.rail_height = rail_height
        self.narrow_width = narrow_width
        self.narrow_height = narrow_height
        self.length = length

        self.lons_z_adjust = 0.3038
        self.lons_z_step = SBQuantity.angstrom(1.5737 + self.lons_z_adjust)

    def calculate_distance(self, atom1, atom2):
        """Calculate distance between two atoms"""
        x1, y1, z1 = atom1.getX().angstrom.value, atom1.getY().angstrom.value, atom1.getZ().angstrom.value
        x2, y2, z2 = atom2.getX().angstrom.value, atom2.getY().angstrom.value, atom2.getZ().angstrom.value
        return math.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)

    def get_max_y_of_lonsdaleite(self, ingot_grid, ingot_cell_structure):
        """Get maximum y-coordinate of a lonsdaleite ingot"""
        (x_cells, y_cells, z_cells) = ingot_cell_structure
        max_y = None

        # we only look at the first 9 atoms on the top face
        # if the ingot is smaller than that, we just look at whatever is there
        # since the ingot's face has two separate heights, and we may be missing atoms
        # near the corners, checking 9 positions should be sufficient
        for x_idx in range(3):
            for z_idx in range(3):
                atom = ingot_grid.get((x_idx, y_cells - 1, z_idx))
                if atom is not None:
                    position = atom.getPosition()
                    if max_y is None or position.y > max_y:
                        max_y = position.y

        if max_y is None:
            raise ValueError("Could not determine max y of lonsdaleite ingot")
        return max_y

    def translate_atoms(self, atoms, translation_vector):
        """Translate a list of atoms by a given vector"""
        for atom in atoms:
            position = atom.getPosition()
            position += translation_vector
            atom.setPosition(position)

    def bond_ingots(self, bottom_ingot, bottom_grid, top_ingot, top_grid, structural_model):
        """Bond lonsdalite ingots along their faces"""
        SAMSON.beginHolding("Bond ingots")

        (bottom_x_cells, bottom_y_cells, bottom_z_cells) = bottom_grid
        (top_x_cells, top_y_cells, top_z_cells) = top_grid
        bottom_y_cells -= 1    # adjust to last layer index (y only)
        top_y_cells -= 1     # adjust to last layer index (y only)

        # Find and bond closest atoms between bottom and top ingots
        # central atoms (x = top_*_cells/2) along both faces align (by assumption)
        # z-length is the same for both ingots (by assumption)

        wide_on_top = bottom_x_cells <= top_x_cells

        narrow_ingot  = bottom_ingot   if wide_on_top else top_ingot
        narrow_grid   = bottom_grid    if wide_on_top else top_grid
        narrow_face_y = bottom_y_cells if wide_on_top else 0
        wide_ingot    = top_ingot      if wide_on_top else bottom_ingot
        wide_grid     = top_grid       if wide_on_top else bottom_grid
        wide_face_y   = 0              if wide_on_top else bottom_y_cells

        # only seek to bond the first layer
        for z_idx in range(narrow_grid[2] + 1):
            for ring_idx in range(2 * narrow_grid[0] + 1):
                narrow_face_atom = narrow_ingot.get((ring_idx, narrow_face_y, z_idx))
                if narrow_face_atom is None:
                    continue

                wide_ring_idx = ring_idx + (2*wide_grid[0] - 2*narrow_grid[0])/2
                wide_face_atom = wide_ingot.get((wide_ring_idx, wide_face_y, z_idx))
                if wide_face_atom is None:
                    continue

                dist = self.calculate_distance(narrow_face_atom, wide_face_atom)
                if dist < (self.lons_z_step.angstrom.value * 1.1):
                    min_dist = dist
                    closest_cnt_atom = wide_face_atom
                else:
                    continue

                bond = SBBond(narrow_face_atom, wide_face_atom, 1.0)
                SAMSON.hold(bond)
                bond.create()
                structural_model.addChild(bond)

        SAMSON.endHolding()

    def layers_are_offset(self, bottom_cell_structure, top_cell_structure):
        """Check if two lonsdaleite ingots have offset layers"""
        (bottom_x_cells, bottom_y_cells, bottom_z_cells) = bottom_cell_structure
        (top_x_cells, top_y_cells, top_z_cells) = top_cell_structure

        bottom_y_layer_mod = (-1) if (bottom_y_cells - 1) % 2 == 0 else 1
        bottom_z_layer_mod = (-1) if (bottom_x_cells % 2 == 0) else 1

        top_y_layer_mod = -1
        top_z_layer_mod = (-1) if (top_x_cells % 2 == 0) else 1


        # if 
        if (bottom_y_layer_mod * bottom_z_layer_mod) == (top_y_layer_mod * top_z_layer_mod):
            return True
        return False

    def generate_rail(self):
        """Main method to generate complete guide rail"""
        print(f"Generating guide rail: length={self.length}Å")

        SAMSON.beginHolding("Create Structural Model")
        structural_model = SBStructuralModel()
        structural_model.name = f"guiderail L={self.length}Å"
        structural_model.create()
        SAMSON.endHolding()

        anchor_ingot = lonsdaleite_ingot.LonsdaleiteIngot(self.anchor_width, self.anchor_height, self.length)
        anchor_model, anchor_atoms, anchor_grid, anchor_cell_structure = anchor_ingot.generate_ingot()
        anchor_model.name = "anchor ingot"
        anchor_top = self.get_max_y_of_lonsdaleite(anchor_grid, anchor_cell_structure)

        narrow_ingot = lonsdaleite_ingot.LonsdaleiteIngot(self.narrow_width, self.narrow_height, self.length)
        narrow_model, narrow_atoms, narrow_grid, narrow_cell_structure = narrow_ingot.generate_ingot()
        narrow_model.name = "narrow ingot"
        narrow_bottom_y = anchor_top + self.lons_z_step
        shift_narrow = self.layers_are_offset(anchor_cell_structure, narrow_cell_structure)
        if shift_narrow:
            raise ValueError("Anchor and narrow ingots have offset layers; please adjust dimensions.")
        self.translate_atoms(narrow_atoms, samson.SBPhysicalVector3(SBQuantity.angstrom(0),
                                                                 SBQuantity.angstrom(narrow_bottom_y),
                                                                 SBQuantity.angstrom(0)))
        narrow_top_y = self.get_max_y_of_lonsdaleite(narrow_grid, narrow_cell_structure)

        rail_ingot = lonsdaleite_ingot.LonsdaleiteIngot(self.rail_width, self.rail_height, self.length)
        rail_model, rail_atoms, rail_grid, rail_cell_structure = rail_ingot.generate_ingot()
        rail_model.name = "rail ingot"
        rail_bottom_y = narrow_top_y + self.lons_z_step
        shift_rail = self.layers_are_offset(narrow_cell_structure, rail_cell_structure)
        if shift_rail:
            raise ValueError("Narrow and rail ingots have offset layers; please adjust dimensions.")
        self.translate_atoms(rail_atoms, samson.SBPhysicalVector3(SBQuantity.angstrom(0),
                                                                 SBQuantity.angstrom(rail_bottom_y),
                                                                 SBQuantity.angstrom(0)))
        self.bond_ingots(anchor_grid, anchor_cell_structure, narrow_grid, narrow_cell_structure, structural_model)
        self.bond_ingots(narrow_grid, narrow_cell_structure, rail_grid, rail_cell_structure, structural_model)

        SAMSON.beginHolding("Add ingots to model")
        structural_model.addChild(anchor_model)
        structural_model.addChild(rail_model)
        structural_model.addChild(narrow_model)
        SAMSON.endHolding()

        # Add to document
        document = SAMSON.getActiveDocument()
        SAMSON.beginHolding(f"Create guide rail L={self.length}Å")
        SAMSON.hold(structural_model)
        document.addChild(structural_model)
        SAMSON.endHolding()

        print("Guide rail generation complete!")
        return structural_model

# Example usage
if __name__ == "__main__":
    length = float(input("Guiderail length (Angstroms): "))
    anchor_width = 30 #float(input("Anchor ingot width (Angstroms): "))
    anchor_height = 16 #float(input("Anchor ingot height (Angstroms): "))
    narrow_width = 10 #float(input("Narrow ingot width (Angstroms): "))
    narrow_height = 16 #float(input("Narrow ingot height (Angstroms): "))
    rail_width = 20 #float(input("Rail ingot width (Angstroms): "))
    rail_height = 10 #float(input("Rail ingot height (Angstroms): "))

    generator = GuideRailGenerator(anchor_width, anchor_height, rail_width, rail_height, narrow_width, narrow_height, length)
    rail = generator.generate_rail()