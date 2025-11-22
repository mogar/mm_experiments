"""
Lonsdaleite Ingot Generator for SAMSON
"""

import math
from samson import *

class LonsdaleiteIngot():
    """
    Generates lonsdaleite ingots.

    Parameters:
    - n: Chiral index for armchair nanotube (M=0 means armchair: (n,n))
    - z_length: Length of nanotube in angstroms
    - num_teeth: Number of gear teeth around the circumference
    """

    def __init__(self, x_width, y_height, z_length=100.0):
        self.x_width = x_width
        self.y_height = y_height
        self.z_length = z_length

        # Carbon-carbon bond length in graphene/nanotubes (angstroms)
        self.acc = 1.413

    def generate_ingot(self, structural_model = None):
        """
        Generate a single lonsdaleite ingot structure.
        Lonsdaleite is hexagonal diamond with ABAB stacking.
        """
        atoms = []
        atom_grid = {}

        # Create structural model
        if structural_model is None:
            structural_model = SBStructuralModel()
            structural_model.name = f"Lonsdaleite Ingot (w={self.x_width}, h={self.y_height}, l={self.z_length})"
            structural_model.create()

        axial_length = self.z_length

        x_step_size = self.acc*math.sqrt(3)/2

        # Number of layers for lonsdaleite structure
        layer_adjust = self.acc/6
        layer_height = self.acc + 2 * layer_adjust
        num_layers = int(self.y_height / layer_height)
        num_x_cells = int(self.x_width / (x_step_size*2))
        num_z_cells = int(self.z_length / (self.acc * 1.5))

        print(f"cell size: width={num_x_cells}, height={num_layers}, length={num_z_cells}")

        # Create lonsdaleite hexagonal layers
        for layer in range(num_layers):
            for z_idx in range(num_z_cells + 1):
                start_x = -num_x_cells * x_step_size
                for ring_idx in range(2 * num_x_cells + 1):
                    if (z_idx == 0):
                        if (ring_idx == 0) or (ring_idx == 2*num_x_cells):
                            # skip first and last carbons of first row because they'll dangle
                            continue

                    x = start_x + ring_idx * x_step_size
                    y = layer * layer_height
                    z = self.acc/4 + z_idx * self.acc * 1.5

                    # modify y up or down based on grid
                    y_layer_mod = (-1) if layer % 2 == 0 else 1
                    z_layer_mod = (-1) if (ring_idx % 2 == z_idx % 2) else 1

                    point_x = x
                    point_y = y + y_layer_mod * z_layer_mod * layer_adjust
                    point_z = z + z_layer_mod * self.acc/4

                    atom = SBAtom(
                        SBElement.Carbon,
                        SBQuantity.angstrom(point_x),
                        SBQuantity.angstrom(point_y),
                        SBQuantity.angstrom(point_z)
                    )
                    atoms.append(atom)
                    atom_grid[(ring_idx, layer, z_idx)] = atom
                    SAMSON.hold(atom)
                    atom.create()
                    structural_model.addChild(atom)

        cell_structure = (num_x_cells, num_layers, num_z_cells)
        self.create_ingot_bonds(atoms, atom_grid, cell_structure, structural_model)
        return structural_model, atoms, atom_grid, cell_structure

    def create_ingot_bonds(self, ingot_atoms, ingot_grid, cell_structure, structural_model):
        """Create bonds within ingot structure"""
        SAMSON.beginHolding("Create ingot bonds")

        (num_x_cells, num_y_cells, num_z_cells) = cell_structure

        for y_idx in range(num_y_cells):
            for z_idx in range(num_z_cells + 1):
                for ring_idx in range(2 * num_x_cells + 1):
                    current_atom = ingot_grid.get((ring_idx, y_idx, z_idx))

                    if current_atom is None:
                        continue

                    # Bond to next atom in ring (circumferential)
                    next_ring_idx = ring_idx + 1
                    if next_ring_idx < (2*num_x_cells + 1):
                        next_atom = ingot_grid.get((next_ring_idx, y_idx, z_idx))

                        if next_atom is not None:
                            bond = SBBond(current_atom, next_atom, 1.0)
                            SAMSON.hold(bond)
                            bond.create()
                            structural_model.addChild(bond)

                    # Bond to atoms in next z-layer (axial)
                    if (ring_idx % 2 != z_idx % 2) and (z_idx < num_z_cells):
                        ring_offset = 0 #(-1) if (z_idx % 2 == 0) else 1
                        bonding_ring_idx = ring_idx + ring_offset
                        # Connect to corresponding atom in next layer
                        next_z_atom = ingot_grid.get((bonding_ring_idx, y_idx, z_idx + 1))
                        if next_z_atom is not None:
                            bond = SBBond(current_atom, next_z_atom, 1.0)
                            SAMSON.hold(bond)
                            bond.create()
                            structural_model.addChild(bond)

                    # layer-to-layer bonds for y
                    y_layer_mod = (-1) if y_idx % 2 == 0 else 1
                    z_layer_mod = (-1) if (ring_idx % 2 == z_idx % 2) else 1
                    if y_layer_mod * z_layer_mod > 0:
                        next_y_atom = ingot_grid.get((ring_idx, y_idx + 1, z_idx))
                        if next_y_atom is not None:
                            bond = SBBond(current_atom, next_y_atom, 1.0)
                            SAMSON.hold(bond)
                            bond.create()
                            structural_model.addChild(bond)

        SAMSON.endHolding()

if __name__ == "__main__":
    x_width = float(input("Ingot width (Angstroms): "))
    y_height = float(input("Ingot height (Angstroms): "))
    z_length = float(input("Ingot length (Angstroms): "))
    generator = LonsdaleiteIngot(x_width, y_height, z_length)
    structural_model, _, _, _ = generator.generate_ingot()

    # Add to document
    document = SAMSON.getActiveDocument()
    SAMSON.beginHolding(f"Create Lonsdaleite Ingot width={x_width}, height={y_height}, length={z_length}")
    SAMSON.hold(structural_model)
    document.addChild(structural_model)
    SAMSON.endHolding()
