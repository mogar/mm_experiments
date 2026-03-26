"""
Lonsdaleite Bar With evenly spaced hexagonal sockets.
"""

import math
from samson import *

import lonsdaleite_ingot
import utilities

class BarWithHexSockets:
    def __init__(self, bar_length, bar_width, bar_height, num_sockets):
        self.bar_length = bar_length
        self.bar_width = bar_width
        self.bar_height = bar_height
        self.socket_radius = 20.5 # angstroms
        self.num_sockets = num_sockets

    def excavate_socket(self, bar_model, socket_position, hex_radius):
        """Excavate a hexagonal socket at the specified position.

        - Determines geometry of a 2D hexagon centered at socket_position with radius hex_radius
        - Iterates through all atoms once, checks if x,z are inside the hexagon
        - Removes all atoms inside the hexagon and their bonds
        """
        if isinstance(socket_position, (tuple, list)) and len(socket_position) == 3:
            cx, cy, cz = socket_position
        elif hasattr(socket_position, 'x') and hasattr(socket_position, 'y') and hasattr(socket_position, 'z'):
            cx = socket_position.x.angstrom.value
            cy = socket_position.y.angstrom.value
            cz = socket_position.z.angstrom.value
        else:
            raise TypeError("socket_position must be tuple/list of 3 floats or position-like object")

        # Find the closest atom within the cylinder to use as the socket center
        cylinder_radius = 1.5  # angstroms
        center_x_min = cx - cylinder_radius
        cylinder_x_max = cx + cylinder_radius
        cylinder_z_min = cz - cylinder_radius
        cylinder_z_max = cz + cylinder_radius


        possible_center_xs = set(bar_model.getNodes(f'node.type atom and a.x {center_x_min}A:{cylinder_x_max}A'))
        possible_center_zs = set(bar_model.getNodes(f'node.type atom and a.z {cylinder_z_min}A:{cylinder_z_max}A'))
        possible_center_atoms = list(possible_center_xs & possible_center_zs)

        closest_atom = None
        min_distance = float('inf')

        for atom in possible_center_atoms:
            ax = atom.getX().angstrom.value
            az = atom.getZ().angstrom.value
            dx = ax - cx
            dz = az - cz
            distance = math.sqrt(dx * dx + dz * dz)

            if distance < min_distance:
                min_distance = distance
                closest_atom = atom

        # Use the closest atom's position as the socket center
        if closest_atom is not None:
            cx = closest_atom.getX().angstrom.value
            cz = closest_atom.getZ().angstrom.value

        bounding_box_radius = hex_radius * 1.25 # add some padding to ensure we get all atoms in the hexagon
        x_min = cx - bounding_box_radius
        x_max = cx + bounding_box_radius
        z_min = cz - bounding_box_radius
        z_max = cz + bounding_box_radius

        xs = set(bar_model.getNodes(f'node.type atom and a.x {x_min}A:{x_max}A'))
        zs = set(bar_model.getNodes(f'node.type atom and a.z {z_min}A:{z_max}A'))
        all_atoms = list(xs & zs)

        if not all_atoms:
            return

        atoms_to_remove = []
        for atom in all_atoms:
            dx = atom.getX().angstrom.value - cx
            dz = atom.getZ().angstrom.value - cz

            # Check if the atom's x,z position is inside the hexagon
            hex_r = utilities.hexagon_distance(dz, dx) # reverse dz, dx because hexagon_distance assumes hexagon is oriented with flat sides on top/bottom, but we want it rotated 90 degrees
            if hex_r <= hex_radius:
                atoms_to_remove.append(atom)

        # Remove all atoms inside the hexagon
        for atom in atoms_to_remove:
            utilities.remove_atom_and_bonds(atom)


    def generate(self):
        # Generate the main bar
        print(f"Bar(length={self.bar_length}, width={self.bar_width}, height={self.bar_height})")
        
        SAMSON.beginHolding("Create Structural Model")
        structural_model = SBStructuralModel()
        structural_model.name = f"bar with sockets"
        structural_model.create()
        SAMSON.endHolding()

        # Calculate the spacing between sockets
        spacing = (self.bar_length - 2 * self.socket_radius) / (self.num_sockets)
        start_x = -(self.num_sockets - 1) * spacing / 2.0

        # Generate bar
        bar_ingot = lonsdaleite_ingot.LonsdaleiteIngot(self.bar_length, self.bar_height, self.bar_width)
        bar_model, bar_atoms, bar_grid, bar_cell_structure = bar_ingot.generate_ingot()
        bar_model.name = "bar ingot"

        # Generate the sockets and excavate each one from the bar model
        for i in range(self.num_sockets):
            x_position = start_x + i * spacing
            socket_center = (x_position, self.bar_height / 2.0, self.bar_width / 2.0)
            print(f"HexSocket(radius={self.socket_radius}, position=({socket_center[0]}, {socket_center[1]}, {socket_center[2]}))")
            self.excavate_socket(bar_model, socket_center, hex_radius=self.socket_radius)

        SAMSON.beginHolding("Add ingots to model")
        structural_model.addChild(bar_model)
        SAMSON.endHolding()

        # Add to document
        document = SAMSON.getActiveDocument()
        SAMSON.beginHolding(f"Create guide clip L={self.bar_length}Å")
        SAMSON.hold(structural_model)
        document.addChild(structural_model)
        SAMSON.endHolding()

        return bar_model


# Example usage
if __name__ == "__main__":
    length = float(input("Guiderail length (Angstroms): "))
    num_sockets = int(input("Number of sockets: "))
    bar_width = 48 #float(input("Anchor ingot width (Angstroms): "))
    bar_height = 15 #float(input("Anchor ingot height (Angstroms): "))

    generator = BarWithHexSockets(length, bar_width, bar_height, num_sockets)
    stock = generator.generate()