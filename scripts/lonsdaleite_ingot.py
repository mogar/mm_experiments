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

        # Carbon-carbon bond length in lonsdaleite (angstroms)
        self.acc = 1.59
        
        # bond lengths (2D)
        # planar y-vertical bonds = 1.59 angstroms
        # planar y-vertical distance between columns = 0.575 angstroms
        # planar y-vertical distance between lower and upper points of hexagon = bond length + 2 * 0.575 = 2.74 angstroms
        self.y_adjust = 0.575/2
        self.y_step = 1.59 + 0.575
        # planar x-horizontal distance between columns = 2.513 angstroms
        self.x_step = 2.513/2
        # planar z-axial distance between layers = 2.909 angstroms
        # planar z-axial bond length (ignoring off-axis length) = 1.4397 angstroms
        # planar z-axial distance between lower and upper points of hexagon = bond length + 2 * 0.73465 = 2.909 angstroms
        self.z_step = 1.4397 + 0.73465
        self.z_adjust = 0.73465/2

    def generate_ingot(self):
        """
        Generate a single lonsdaleite ingot structure.
        Lonsdaleite is hexagonal diamond with ABAB stacking.
        """
        atoms = []
        atom_grid = {}

        molecular_model = SBMolecule()
        molecular_model.name = f"Lonsdaleite Ingot (w={self.x_width}, h={self.y_height}, l={self.z_length})"
        molecular_model.create()

        axial_length = self.z_length

        # Number of layers for lonsdaleite structure
        num_layers = int(self.y_height / self.y_step) + 1
        num_x_cells = int(self.x_width / (self.x_step*2))
        num_z_cells = int(self.z_length / self.z_step)

        print(f"cell size: width={num_x_cells}, height={num_layers}, length={num_z_cells}")

        # Create lonsdaleite hexagonal layers
        for layer in range(num_layers):
            for z_idx in range(num_z_cells + 1):
                start_x = -num_x_cells * self.x_step
                for ring_idx in range(2 * num_x_cells + 1):
                    if (z_idx == 0):
                        if (ring_idx == 0) or (ring_idx == 2*num_x_cells):
                            # skip first and last carbons of first row because they'll dangle
                            continue

                    x = start_x + ring_idx * self.x_step
                    y = layer * self.y_step
                    z = self.z_adjust + z_idx * self.z_step

                    # modify y up or down based on grid
                    y_layer_mod = (-1) if layer % 2 == 0 else 1
                    z_layer_mod = (-1) if (ring_idx % 2 == z_idx % 2) else 1

                    point_x = x
                    point_y = y + y_layer_mod * z_layer_mod * self.y_adjust
                    point_z = z + z_layer_mod * self.z_adjust

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
                    molecular_model.addChild(atom)

        cell_structure = (num_x_cells, num_layers, num_z_cells)
        self.create_ingot_bonds(molecular_model, atoms, atom_grid, cell_structure)
        return molecular_model, atoms, atom_grid, cell_structure

    def create_ingot_bonds(self, molecular_model, ingot_atoms, ingot_grid, cell_structure):
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
                            molecular_model.addChild(bond)

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
                            molecular_model.addChild(bond)

                    # layer-to-layer bonds for y
                    y_layer_mod = (-1) if y_idx % 2 == 0 else 1
                    z_layer_mod = (-1) if (ring_idx % 2 == z_idx % 2) else 1
                    if y_layer_mod * z_layer_mod > 0:
                        next_y_atom = ingot_grid.get((ring_idx, y_idx + 1, z_idx))
                        if next_y_atom is not None:
                            bond = SBBond(current_atom, next_y_atom, 1.0)
                            SAMSON.hold(bond)
                            bond.create()
                            molecular_model.addChild(bond)

        SAMSON.endHolding()

if __name__ == "__main__":
    x_width = float(input("Ingot width (Angstroms): "))
    y_height = float(input("Ingot height (Angstroms): "))
    z_length = float(input("Ingot length (Angstroms): "))
    generator = LonsdaleiteIngot(x_width, y_height, z_length)
    molecular_model, _, _, _ = generator.generate_ingot()

    SAMSON.beginHolding("Create Lonsdaleite Ingot Model")
    structural_model = SBStructuralModel()
    structural_model.name = f"Lonsdaleite Ingot (w={x_width}, h={y_height}, l={z_length})"
    structural_model.create()
    structural_model.addChild(molecular_model)
    SAMSON.endHolding()

    # Add to document
    document = SAMSON.getActiveDocument()
    SAMSON.beginHolding(f"Create Lonsdaleite Ingot width={x_width}, height={y_height}, length={z_length}")
    SAMSON.hold(structural_model)
    document.addChild(structural_model)
    SAMSON.endHolding()
