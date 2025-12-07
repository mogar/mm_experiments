
import math, operator
from samson import *


def get_selected_carbon_atoms(self):
    """Get all selected carbon atoms from the active document"""
    # Get all selected atoms
    selected_atoms = SAMSON.getNodes('node.type atom and node.selected')

    # Filter for carbon atoms only
    carbon_atoms = []
    for atom in selected_atoms:
        if atom.elementType == SBElement.Carbon:
            carbon_atoms.append(atom)

    return carbon_atoms

def calculate_distance(atom1, atom2):
    """Calculate Euclidean distance between two atoms"""
    if atom1 is None or atom2 is None:
        # "Cannot calculate distance: one or both atoms are None"
        return float('inf')
    if not isinstance(atom1, SBAtom) or not isinstance(atom2, SBAtom):
        # "Both inputs must be SBAtom instances"
        return float('inf')
    x1, y1, z1 = atom1.getX().angstrom.value, atom1.getY().angstrom.value, atom1.getZ().angstrom.value
    x2, y2, z2 = atom2.getX().angstrom.value, atom2.getY().angstrom.value, atom2.getZ().angstrom.value
    return math.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)

def translate_atoms(atoms, translation_vector):
    """Translate a list of atoms by a given vector"""
    for atom in atoms:
        position = atom.getPosition()
        position += translation_vector
        atom.setPosition(position)

def atoms_already_bonded(atom1, atom2):
    """Check if two atoms are already bonded"""
    # Get all bonds connected to atom1
    bonds_indexer = atom1.getNodes('n.t b')

    for bond in bonds_indexer:
        # Check if this bond connects to atom2
        if bond.getOppositeAtom(atom1) == atom2:
            return True

    return False

def remove_atom_and_bonds(atom):
    parent = atom.getParent()
    # collect bonds attached to the atom
    bonds = []
    for child in parent.getChildren():
        if isinstance(child, SBBond) and (child.leftAtom is atom or child.rightAtom is atom):
            bonds.append(child)

    with SAMSON.holding("Remove atom and its bonds"):
        for b in bonds:
            b.getParent().removeChild(b)  # remove each bond
        parent.removeChild(atom)          # finally remove the atom

def remove_dangling_carbons(molecular_model):
    """Remove carbon atoms with only one bond from the molecular model"""
    atoms_to_remove = []
    for atom in molecular_model.getChildren():
        if isinstance(atom, SBAtom) and atom.elementType == SBElement.Carbon:
            # Count bonds connected to this atom
            bond_count = 0
            for child in molecular_model.getChildren():
                if isinstance(child, SBBond) and (child.leftAtom is atom or child.rightAtom is atom):
                    bond_count += 1
                    if bond_count > 1:
                        break  # No need to count further
            if bond_count <= 1:
                atoms_to_remove.append(atom)

    with SAMSON.holding("Remove dangling carbon atoms"):
        for atom in atoms_to_remove:
            remove_atom_and_bonds(atom)

def get_x_of_lonsdaleite(ingot_grid, ingot_cell_structure, comparison=operator.lt):
    """Get minimum x coordinate of lonsdaleite ingot."""
    (x_cells, y_cells, z_cells) = ingot_cell_structure
    min_x = None

    # we only look at the first 9 atoms on the face
    # if the ingot is smaller than that, we just look at whatever is there
    # since the ingot's face has two separate heights, and we may be missing atoms
    # near the corners, checking 9 positions should be sufficient
    for y_idx in range(3):
        for z_idx in range(3):
            atom = ingot_grid.get((0, y_idx, z_idx))
            if atom is not None:
                position = atom.getPosition()
                if min_x is None or comparison(position.x, min_x):
                    min_x = position.x

    if min_x is None:
        raise ValueError("Could not determine extreme x of lonsdaleite ingot")
    return min_x

def get_y_of_lonsdaleite(ingot_grid, ingot_cell_structure, comparison=operator.lt):
    """Get extremum y coordinate of lonsdaleite ingot."""
    (x_cells, y_cells, z_cells) = ingot_cell_structure
    extreme_y = None

    if operator.lt == comparison:
        y_idx = 0
    else:
        y_idx = y_cells - 1

    # we only look at the first 3 rows of atoms on the face
    # if the ingot is smaller than that, we just look at whatever is there
    # since the ingot's face has two separate heights, and we may be missing atoms
    # near the edges, checking a few rows should be sufficient
    for x_idx in range(2*x_cells + 1):
        for z_idx in range(z_cells + 1):
            atom = ingot_grid.get((x_idx, y_idx, z_idx))
            if atom is not None:
                position = atom.getPosition()
                if extreme_y is None or comparison(position.y, extreme_y):
                    extreme_y = position.y

    if extreme_y is None:
        raise ValueError("Could not determine extreme y of lonsdaleite ingot")
    return extreme_y


def lons_layers_are_offset(bottom_cell_structure, top_cell_structure):
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