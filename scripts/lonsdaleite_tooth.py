"""
Lonsdaleite Tooth Structure Creator for SAMSON
Generates a single tooth structure made of lonsdaleite. Tooth structure is semi-triangular in
cross-section, with a flat base and pointed top. Tapering is achieved as uniformly as possible
along the height of the tooth, subject to the constraint that the tooth remains symmetric
and the lonsdaleite unit cells are preserved.
"""

import math
from samson import *

import utilities
import lonsdaleite_ingot

class LonsdaleiteTooth:
    """
    Generates a single lonsdaleite tooth structure.

    Parameters:
    - tooth_width: Width of the tooth at its base (angstroms)
    - tooth_height: Height of the tooth (angstroms)
    - tooth_length: Length of the tooth along the z-axis (angstroms)
    """

    def __init__(self, tooth_width, tooth_height, tooth_length=12.0):
        self.tooth_width = tooth_width
        self.tooth_height = tooth_height
        self.tooth_length = tooth_length

        # Carbon-carbon bond length in lonsdaleite (angstroms)
        self.acc = 1.59
        
        # bond lengths (2D)
        self.y_adjust = 0.575/2
        self.y_step = 1.59 +self.y_adjust*2
        self.x_step = 2.513/2
        self.z_adjust = 0.3038
        self.z_step = 1.5297 + self.z_adjust*2

    def generate_tooth(self):
        """
        Generate a single lonsdaleite tooth structure.
        Lonsdaleite is hexagonal diamond with ABAB stacking.
        Tooth tapers from self.tooth_width at base to a point at the top.
        """
        atoms = []
        atom_grid = {}

        molecular_model = SBMolecule()
        molecular_model.name = f"Lonsdaleite Tooth (w={self.tooth_width}, h={self.tooth_height}, l={self.tooth_length})"
        molecular_model.create()

        axial_length = self.tooth_length

        # Number of layers for lonsdaleite structure
        num_layers = int(self.tooth_height / self.y_step) + 1
        num_x_cells = int(self.tooth_width / (self.x_step*2))
        num_z_cells = int(self.tooth_length / self.z_step)

        print(f"cell size: width={num_x_cells}, height={num_layers}, length={num_z_cells}")

        # Create lonsdaleite hexagonal layers
        for layer in range(num_layers):
            layer_y = layer * self.y_step
            layer_width = self.tooth_width * (1 - ((layer - 1) / num_layers))
            
            layer_x_cells = 2*int(layer_width / (self.x_step*2))
            ignore_ring_x_low = (2*num_x_cells - layer_x_cells) // 2
            ignore_ring_x_high = 2*num_x_cells - ignore_ring_x_low

            for z_idx in range(num_z_cells + 1):
                start_x = -num_x_cells * self.x_step
                for ring_idx in range(2 * num_x_cells + 1):
                    if ring_idx < ignore_ring_x_low or ring_idx > ignore_ring_x_high:
                        continue  # skip atoms outside tapered width

                    if (z_idx == 0):
                        if (ring_idx == 0) or (ring_idx == 2*num_x_cells):
                            continue  # skip edge atoms at z=0

                    y_layer_mod = (-1) if layer % 2 == 0 else 1
                    z_layer_mod = (-1) if (ring_idx % 2 == z_idx % 2) else 1

                    atom_x = start_x + ring_idx * self.x_step
                    atom_y = layer_y + y_layer_mod * z_layer_mod * self.y_adjust
                    atom_z = self.z_adjust + z_idx * self.z_step + z_layer_mod * self.z_adjust

                    atom = SBAtom(
                        SBElement.Carbon,
                        SBQuantity.angstrom(atom_x),
                        SBQuantity.angstrom(atom_y),
                        SBQuantity.angstrom(atom_z)
                    )
                    atoms.append(atom)
                    atom_grid[(ring_idx, layer, z_idx)] = atom
                    SAMSON.hold(atom)
                    atom.create()
                    molecular_model.addChild(atom)

        cell_structure = (num_x_cells, num_layers, num_z_cells)
        lonsdaleite_ingot.create_lonsdaleite_bonds(molecular_model, atoms, atom_grid, cell_structure)
        return molecular_model, atoms, atom_grid, cell_structure

if __name__ == "__main__":
    x_width = float(input("Tooth width (Angstroms): "))
    y_height = float(input("Tooth height (Angstroms): "))
    z_length = float(input("Tooth length (Angstroms): "))
    generator = LonsdaleiteTooth(x_width, y_height, z_length)
    molecular_model, _, _, _ = generator.generate_tooth()

    SAMSON.beginHolding("Create Lonsdaleite Tooth Model")
    structural_model = SBStructuralModel()
    structural_model.name = f"Lonsdaleite Tooth (w={x_width}, h={y_height}, l={z_length})"
    structural_model.create()
    structural_model.addChild(molecular_model)
    SAMSON.endHolding()

    # Add to document
    document = SAMSON.getActiveDocument()
    SAMSON.beginHolding(f"Create Lonsdaleite Tooth width={x_width}, height={y_height}, length={z_length}")
    SAMSON.hold(structural_model)
    document.addChild(structural_model)
    SAMSON.endHolding()
