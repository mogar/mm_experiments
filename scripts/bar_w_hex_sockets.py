"""
Lonsdaleite Bar With evenly spaced hexagonal sockets.
"""

from samson import *

import lonsdaleite_ingot
import utilities

class BarWithHexSockets:
    def __init__(self, bar_length, bar_width, bar_height, num_sockets):
        self.bar_length = bar_length
        self.bar_width = bar_width
        self.bar_height = bar_height
        self.socket_radius = 11 #socket_radius
        self.num_sockets = num_sockets

    def excavate_socket(self, bar_model, socket_position):
        """Excavate a hexagonal socket at the specified position."""
        # This is a placeholder for the actual excavation logic.
        # You would need to implement the logic to remove atoms from the bar_model
        # based on the geometry of the hexagonal socket and its position.
        pass

    def generate(self):
        # Generate the main bar
        bar = f"Bar(length={self.bar_length}, width={self.bar_width}, height={self.bar_height})"
        
        # Calculate the spacing between sockets
        spacing = (self.bar_length - 2 * self.socket_radius) / (self.num_sockets - 1)
        
        # Generate bar
        bar_ingot = lonsdaleite_ingot.LonsdaleiteIngot(self.bar_height, self.bar_width, self.bar_length)
        bar_model, bar_atoms, bar_grid, bar_cell_structure = bar_ingot.generate_ingot()
        bar_model.name = "bar ingot"

        # Generate the sockets
        sockets = []
        for i in range(self.num_sockets):
            x_position = self.socket_radius + i * spacing
            socket = f"HexSocket(radius={self.socket_radius}, position=({x_position}, {self.bar_width / 2}, {self.bar_height / 2}))"
            sockets.append(socket)
        
        return bar, sockets




# Example usage
if __name__ == "__main__":
    length = float(input("Guiderail length (Angstroms): "))
    bar_width = 20 #float(input("Anchor ingot width (Angstroms): "))
    bar_height = 8 #float(input("Anchor ingot height (Angstroms): "))

    generator = BarWithHexSockets(length, bar_width, bar_height, 4)
    stock, _ = generator.generate()