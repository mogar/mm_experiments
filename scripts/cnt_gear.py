"""
Carbon Nanotube Gear Generator for SAMSON
Generates armchair (n,n) nanotubes with lonsdaleite teeth for molecular gears
"""

import math
from samson import *

class CNTGearGenerator:
    """
    Generates carbon nanotube-based molecular gears with lonsdaleite teeth.
    
    Parameters:
    - n: Chiral index for armchair nanotube (M=0 means armchair: (n,n))
    - z_length: Length of nanotube in angstroms
    - num_teeth: Number of gear teeth around the circumference
    """
    
    def __init__(self, n=6, z_length=100.0, num_teeth=6):
        self.n = n
        self.z_length = z_length
        self.num_teeth = num_teeth
        
        # Carbon-carbon bond length in graphene/nanotubes (angstroms)
        self.acc = 1.413
        
        # Lonsdaleite lattice parameters
        self.a_lons = 2.52  # Lattice parameter in angstroms
        self.c_lons = 4.12  # c-axis parameter for hexagonal structure
        
    def calculate_cnt_radius(self):
        """Calculate the radius of an armchair (n,n) nanotube"""
        # For armchair nanotubes (n,n), radius = n * sqrt(3) * acc / pi
        return self.n * math.sqrt(3) * self.acc / (2*math.pi)
    
    def generate_armchair_cnt(self):
        """
        Generate an armchair carbon nanotube using SAMSON's coordinate generation.
        For an armchair (n,n) nanotube, we create the structure programmatically.
        """
        radius = self.calculate_cnt_radius()
        
        # Create structural model
        structural_model = SBStructuralModel()
        structural_model.create()
        
        # Number of unit cells along z-axis
        num_z_cells = int(self.z_length / (self.acc * 1.5))
        
        atoms = []
        atom_grid = {}  # For tracking atoms by (ring_index, z_index) for bonding
        
        # Generate atoms
        for z_idx in range(num_z_cells + 1):            
            for ring_idx in range(2 * self.n):
                theta = 2 * math.pi * ring_idx / (2 * self.n)
                
                z = z_idx * self.acc * 1.5

                # Offset every other z-layer by half a ring position
                if z_idx % 2 == 1:
                    theta += 2 * math.pi / (2 * self.n)

                # Offset every other atom on-axis 
                if ring_idx % 2 == 1:
                    # adjust by cos(60)*acc (= acc/2)
                    z = z + self.acc/2
                
                x = radius * math.cos(theta)
                y = radius * math.sin(theta)
                
                # Create carbon atom
                atom = SBAtom(
                    SBElement.Carbon,
                    SBQuantity.angstrom(x),
                    SBQuantity.angstrom(y),
                    SBQuantity.angstrom(z)
                )
                
                atoms.append(atom)
                atom_grid[(ring_idx, z_idx)] = atom
        
        # Add atoms to structural model
        SAMSON.beginHolding("Create CNT atoms")
        for atom in atoms:
            SAMSON.hold(atom)
            atom.create()
            structural_model.addChild(atom)
        SAMSON.endHolding()
        
        # Create bonds
        self.create_cnt_bonds(atoms, atom_grid, num_z_cells, structural_model)
        
        return structural_model, atoms, radius
    
    def create_cnt_bonds(self, atoms, atom_grid, num_z_cells, structural_model):
        """Create bonds for the carbon nanotube"""
        SAMSON.beginHolding("Create CNT bonds")
        
        for z_idx in range(num_z_cells + 1):
            for ring_idx in range(2 * self.n):
                current_atom = atom_grid.get((ring_idx, z_idx))
                
                if current_atom is None:
                    continue
                
                # Bond to next atom in ring (circumferential)
                next_ring_idx = (ring_idx + 1) % (2 * self.n)
                next_atom = atom_grid.get((next_ring_idx, z_idx))
                
                if next_atom is not None:
                    bond = SBBond(current_atom, next_atom, 1.0)
                    SAMSON.hold(bond)
                    bond.create()
                    structural_model.addChild(bond)
                
                # Bond to atoms in next z-layer (axial)
                if (ring_idx % 2 == 1) and (z_idx < num_z_cells):
                    ring_offset = (-1) if (z_idx % 2 == 0) else 1
                    bonding_ring_idx = (ring_idx + ring_offset) % (2 * self.n)
                    # Connect to corresponding atom in next layer
                    next_z_atom = atom_grid.get((bonding_ring_idx, z_idx + 1))
                    if next_z_atom is not None:
                        bond = SBBond(current_atom, next_z_atom, 1.0)
                        SAMSON.hold(bond)
                        bond.create()
                        structural_model.addChild(bond)
        
        SAMSON.endHolding()
    
    def calculate_distance(self, atom1, atom2):
        """Calculate distance between two atoms"""
        x1, y1, z1 = atom1.getX().angstrom.value, atom1.getY().angstrom.value, atom1.getZ().angstrom.value
        x2, y2, z2 = atom2.getX().angstrom.value, atom2.getY().angstrom.value, atom2.getZ().angstrom.value
        return math.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
    
    def generate_lonsdaleite_tooth(self, orientation_angle, radius, tooth_height):
        """
        Generate a single lonsdaleite tooth structure.
        Lonsdaleite is hexagonal diamond with ABAB stacking.
        
        orientation_angle: angle around nanotube axis
        radius: CNT radius for proper attachment
        tooth_height: Radial extension from CNT surface
        """
        atoms = []
        atom_grid = {}
 
        axial_length = self.z_length

        x_step_size = self.acc*math.sqrt(3)/2

        # TODO: use a_lons and c_lons

        # Number of layers for lonsdaleite structure
        layer_height = self.acc # TODO: make this more accurate to account for "bumps"
        num_layers = int(tooth_height / layer_height)
        num_x_cells = 3 # TODO: calculate or pass in
        num_z_cells = int(self.z_length / (self.acc * 1.5))

        # Create lonsdaleite hexagonal layers
        for layer in range(num_layers):
            # Radial distance from CNT center
            r = radius + layer * layer_height
            
            layer_y = r*math.sin(orientation_angle)
            layer_x = r*math.cos(orientation_angle)

            for z_idx in range(num_z_cells + 1):  
                start_x = -num_x_cells * x_step_size + x_step_size/2      
                for ring_idx in range(2 * num_x_cells):
                    x = start_x + ring_idx *x_step_size
                    
                    z = z_idx * self.acc * 1.5

                    # Offset every other z-layer by half a ring position
                    if z_idx % 2 == 1:
                        x += x_step_size

                    # Offset every other atom on-axis 
                    if ring_idx % 2 == 1:
                        # adjust by cos(60)*acc (= acc/2)
                        z = z + self.acc/2
                
                    # modify y up or down based on grid
                    y_layer_mod = (-1) if layer % 2 == 0 else 1
                    z_layer_mod = (-1) if (ring_idx % 2 == z_idx % 2) else 1
                    y = y_layer_mod * z_layer_mod * self.acc/6

                    s = math.sin(orientation_angle + math.pi/2)
                    c = math.cos(orientation_angle + math.pi/2)

                    point_x = layer_x + x * c + y * s
                    point_y = layer_y + x * s + y * c
                    point_z = z

                    atom = SBAtom(
                        SBElement.Carbon,
                        SBQuantity.angstrom(point_x),
                        SBQuantity.angstrom(point_y),
                        SBQuantity.angstrom(point_z)
                    )
                    atoms.append(atom)
                    atom_grid[(ring_idx, layer, z_idx)] = atom

        return atoms, atom_grid, (num_x_cells, num_layers, num_z_cells)
    
    def create_tooth_bonds(self, tooth_atoms, tooth_grid, cell_structure, structural_model):
        """Create bonds within a tooth structure"""
        SAMSON.beginHolding("Create tooth bonds")
        
        (num_x_cells, num_y_cells, num_z_cells) = cell_structure

        for y_idx in range(num_y_cells):
            for z_idx in range(num_z_cells + 1):
                for ring_idx in range(2 * num_x_cells):
                    current_atom = tooth_grid.get((ring_idx, y_idx, z_idx))
                    
                    if current_atom is None:
                        continue
                    
                    # Bond to next atom in ring (circumferential)
                    next_ring_idx = ring_idx + 1
                    if next_ring_idx < 2*num_x_cells:
                        next_atom = tooth_grid.get((next_ring_idx, y_idx, z_idx))
                    
                        if next_atom is not None:
                            bond = SBBond(current_atom, next_atom, 1.0)
                            SAMSON.hold(bond)
                            bond.create()
                            structural_model.addChild(bond)
                    
                    # Bond to atoms in next z-layer (axial)
                    if (ring_idx % 2 == 1) and (z_idx < num_z_cells):
                        ring_offset = (-1) if (z_idx % 2 == 0) else 1
                        bonding_ring_idx = ring_idx + ring_offset
                        # Connect to corresponding atom in next layer
                        next_z_atom = tooth_grid.get((bonding_ring_idx, y_idx, z_idx + 1))
                        if next_z_atom is not None:
                            bond = SBBond(current_atom, next_z_atom, 1.0)
                            SAMSON.hold(bond)
                            bond.create()
                            structural_model.addChild(bond)

                    # layer-to-layer bonds for y
                    y_layer_mod = (-1) if y_idx % 2 == 0 else 1
                    z_layer_mod = (-1) if (ring_idx % 2 == z_idx % 2) else 1
                    if y_layer_mod * z_layer_mod > 0:
                        next_y_atom = tooth_grid.get((ring_idx, y_idx + 1, z_idx))
                        if next_y_atom is not None:
                            bond = SBBond(current_atom, next_y_atom, 1.0)
                            SAMSON.hold(bond)
                            bond.create()
                            structural_model.addChild(bond)
                
        SAMSON.endHolding()
    
    def bond_tooth_to_cnt(self, tooth_atoms, tooth_grid, tooth_grid_structure, cnt_atoms, structural_model):
        """Bond tooth base atoms to nearest CNT surface atoms"""
        SAMSON.beginHolding("Bond tooth to CNT")
        
        (num_x_cells, num_y_cells, num_z_cells) = tooth_grid_structure

        # only seek to bond the first layer
        for z_idx in range(num_z_cells + 1):
            for ring_idx in range(2 * num_x_cells):
                tooth_atom = tooth_grid.get((ring_idx, 0, z_idx))
        
                min_dist = float('inf')
                closest_cnt_atom = None
                
                for cnt_atom in cnt_atoms:
                    dist = self.calculate_distance(tooth_atom, cnt_atom)
                    if dist < min_dist:
                        min_dist = dist
                        closest_cnt_atom = cnt_atom
                
                # Create bond if reasonably close
                # TODO: may want to merge instead
                if closest_cnt_atom and min_dist < 3.0:
                    bond = SBBond(tooth_atom, closest_cnt_atom, 1.0)
                    SAMSON.hold(bond)
                    bond.create()
                    structural_model.addChild(bond)
        
        SAMSON.endHolding()

    def attach_teeth_to_cnt(self, structural_model, cnt_atoms, radius, tooth_height):
        """Attach lonsdaleite teeth to the nanotube at specified positions"""
                
        # TODO: determine angle offset to line up with cnt grid
        tooth_phase_offset = 0

        SAMSON.beginHolding("Add gear teeth")
        
        for tooth_idx in range(self.num_teeth):
            angle = 2 * math.pi * tooth_idx / self.num_teeth + tooth_phase_offset
            
            # Generate tooth structure
            tooth_atoms, tooth_grid, tooth_cell_sizes = self.generate_lonsdaleite_tooth(
                angle,
                radius,
                tooth_height
            )
            
            # Add tooth atoms to structural model
            for atom in tooth_atoms:
                SAMSON.hold(atom)
                atom.create()
                structural_model.addChild(atom)
            
            # Create bonds within tooth structure
            self.create_tooth_bonds(tooth_atoms, tooth_grid, tooth_cell_sizes, structural_model)
            
            # Find and bond closest CNT atoms to tooth base
            self.bond_tooth_to_cnt(tooth_atoms, tooth_grid, tooth_cell_sizes, cnt_atoms, structural_model)
        
        SAMSON.endHolding()
    
    def generate_gear(self):
        """Main method to generate complete CNT gear"""
        print(f"Generating CNT gear: n={self.n}, length={self.z_length}Å, teeth={self.num_teeth}")
        
        # Generate base CNT structure
        structural_model, cnt_atoms, radius = self.generate_armchair_cnt()
        
        # Attach teeth
        tooth_height = 10 # angstroms
        self.attach_teeth_to_cnt(structural_model, cnt_atoms, radius, tooth_height)
        
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
    # Create a gear with n=12 (armchair CNT), 30 angstrom length, 6 teeth
    generator = CNTGearGenerator(n=12, z_length=30.0, num_teeth=6)
    gear = generator.generate_gear()
    