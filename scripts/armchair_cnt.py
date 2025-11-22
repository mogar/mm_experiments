"""
Carbon Nanotube Generator for SAMSON
Generates armchair (n,n) nanotubes
"""

import math
from samson import *

class CNTGenerator:
	"""
	Generates carbon nanotube-based molecular gears with lonsdaleite teeth.
	
	Parameters:
	- n: Chiral index for armchair nanotube (M=0 means armchair: (n,n))
	- z_length: Length of nanotube in angstroms
	"""
	
	def __init__(self, n=6, z_length=100.0):
		self.n = n
		self.z_length = z_length
		
		# Carbon-carbon bond length in graphene/nanotubes (angstroms)
		self.acc = 1.413
				
	def calculate_cnt_radius(self):
		"""Calculate the radius of an armchair (n,n) nanotube"""
		# For armchair nanotubes (n,n), radius = n * sqrt(3) * acc / pi
		return self.n * math.sqrt(3) * self.acc / (2*math.pi)
	
	def generate_armchair_cnt(self, structural_model = None):
		"""
		Generate an armchair carbon nanotube using SAMSON's coordinate generation.
		For an armchair (n,n) nanotube, we create the structure programmatically.
		"""
		radius = self.calculate_cnt_radius()
		
		SAMSON.beginHolding("Create CNT atoms")

		# Create structural model
		if structural_model is None:
			structural_model = SBStructuralModel()
			structural_model.name = f"Armchair CNT (N={self.n}, h={self.z_length})"
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
		for atom in atoms:
			SAMSON.hold(atom)
			atom.create()
			structural_model.addChild(atom)
			
		SAMSON.endHolding()
		
		# Create bonds
		self.create_cnt_bonds(structural_model, atoms, atom_grid, num_z_cells)
		
		return structural_model, atoms, radius
	
	def create_cnt_bonds(self, structural_model, atoms, atom_grid, num_z_cells):
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

if __name__ == "__main__":
	n = int(input("carbon nanotube N: "))
	z_length = int(input("carbon nanotube length (Angstroms): "))
	generator = CNTGenerator(n, z_length)
	structural_model, _, _ = generator.generate_armchair_cnt()

	# Add to document
	document = SAMSON.getActiveDocument()
	SAMSON.beginHolding(f"Create CNT n={n}")
	SAMSON.hold(structural_model)
	document.addChild(structural_model)
	SAMSON.endHolding()
	