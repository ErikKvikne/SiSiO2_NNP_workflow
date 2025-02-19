from scipy.spatial.distance import cdist
from ase.visualize import view
import numpy as np
from pymatgen.core import Structure, Lattice
from pymatgen.core.surface import SlabGenerator, get_symmetrically_distinct_miller_indices
from pymatgen.io.cif import CifWriter
from pymatgen.core.interface import fix_pbc
from pymatgen.io.xyz import XYZ

def fetch_asio2(filename, asio2_shift = [0,0,0]):
    """
    Fetches the amorphous SiO2 structure from an XYZ file and applies a shift.
    Atoms are shifted by a vector specified as a fraction of the lattice parameters.
    
    Parameters:
    - filename: str, the name of the XYZ file containing the structure.
    - asio2_shift: list of three floats, fractional shift vector [x_shift, y_shift, z_shift].
    
    Returns:
    - asio2_structure: pymatgen Structure object representing the shifted SiO2 structure.
    """
    asio2_molecule = XYZ.from_file(filename).molecule

    # Find the maximum coordinate value among all atomic positions
    max_position = max(max(coord) for coord in asio2_molecule.cart_coords)

    # Create a cubic lattice using the largest coordinate value
    lattice = Lattice.from_parameters(max_position, max_position, max_position, 90, 90, 90)

    # Convert to a pymatgen Structure
    asio2_structure = Structure(lattice, asio2_molecule.species, asio2_molecule.cart_coords, coords_are_cartesian=True)

    # Calculate the Cartesian shift vector
    shift_cartesian = np.array(asio2_shift) * np.array([lattice.a, lattice.b, lattice.c])

    # Apply the shift and wrap atoms within the periodic boundaries
    asio2_structure.translate_sites(
        range(len(asio2_structure)),
        shift_cartesian,
        frac_coords=False,
        to_unit_cell=True
    )

    return asio2_structure

def generate_si_slab(miller_index, min_thickness, asio2 = None):
    # Create bulk Si structure, lattice length from materials project
    si_bulk = Structure.from_spacegroup(227, Lattice.cubic(5.44370237), ["Si"], [[0, 0, 0]])

    #ls = get_symmetrically_distinct_miller_indices(si_bulk, 3)

    # Generate slab without vacuum, using the conventional cell and orthogonal lattice
    slabgen = SlabGenerator(
        si_bulk,
        miller_index,
        min_slab_size=min_thickness,
        min_vacuum_size=0,
        center_slab=False,
        primitive=False,
        lll_reduce=True
        #reorient_lattice=True  # Ensure lattice is reoriented to have orthogonal axes
    )
    slab = slabgen.get_slab()

    slab_test = fix_pbc(slab.get_orthogonal_c_slab().get_sorted_structure())

    if asio2 == None:
        scaling_matrix = [1,1,1]
    else:
        # Define the Si lattice dimensions
        si_a, si_b = slab_test.lattice.a, slab_test.lattice.b

        # Define the desired lattice dimensions from the amorphous SiO2
        asio2_a, asio2_b= asio2.lattice.a, asio2.lattice.b

        # Calculate the number of repetitions along each dimension
        repeat_a = int(round(asio2_a / si_a))
        repeat_b = int(round(asio2_b / si_b))


        # Define the scaling matrix for the supercell
        scaling_matrix = [[repeat_a, 0, 0], [0, repeat_b, 0], [0, 0, 1]]



    supercell = slab_test.copy()
    supercell.make_supercell(scaling_matrix)  # Scale the Si slab


    #view(supercell.to_ase_atoms())

    return supercell

def generate_interface(si, asio2, structure_filename, vaccum_separation = 0, vaccum_over = 0):


    # Define a combined lattice that only adjusts the bounding box in the c-axis direction
    combined_lattice = Lattice.from_parameters(
        max(si.lattice.a, asio2.lattice.a),
        max(si.lattice.a, asio2.lattice.b),
        si.lattice.c + asio2.lattice.c + vaccum_separation + vaccum_over,
        90, 90, 90
    )

    # Create the combined structure with the new lattice
    combined_structure = Structure(combined_lattice, si.species, si.cart_coords, coords_are_cartesian=True)
    for site in asio2:
        combined_structure.append(site.specie, site.coords + [0,0, si.lattice.c + vaccum_separation], coords_are_cartesian=True)

    # Output the modified structure to a CIF file
    combined_structure.to(structure_filename, fmt="cif")
    view(combined_structure.to_ase_atoms())

    print("Combined Si structure created with an elongated bounding box in the c-axis.")

def indent(si,asio2, x_int, y_int, si_indent = True, fill_indent = False):
    """
    Does not yet take PBC into concideration when looking for neighbor atoms, works fine when indent is in centre
    """

    unit_cell_length = 5.44370237
    delta = 0.001

    si_si_threshold=2.3 
    si_o_threshold=1.55 

    a_scaling = round(si.lattice.a/unit_cell_length)
    b_scaling = round(si.lattice.b/unit_cell_length)
    c_scaling = round(si.lattice.c/unit_cell_length)


    if si_indent == True:

        lower_bounds = [x_int * unit_cell_length - delta, y_int * unit_cell_length - delta, si.lattice.c - unit_cell_length - delta]
        upper_bounds = [x_int * unit_cell_length + unit_cell_length - delta, y_int * unit_cell_length + unit_cell_length - delta, si.lattice.c - delta]


        # Get Cartesian coordinates for all si atoms
        cart_coords_si = np.array(si.cart_coords)

        # Identify indices of atoms within the specified region
        indices_to_remove = [
            i for i, coord in enumerate(cart_coords_si)
            if all(lower_bounds <= coord) and all(coord <= upper_bounds)
        ]

        # Remove the identified atoms
        si.remove_sites(indices_to_remove)

        if fill_indent == True:

            lower_bounds_asio2 = [x_int * unit_cell_length - delta, y_int * unit_cell_length - delta, asio2.lattice.c - unit_cell_length - delta]
            upper_bounds_asio2 = [x_int * unit_cell_length + unit_cell_length - delta, y_int * unit_cell_length + unit_cell_length - delta, asio2.lattice.c - delta]

            # Get Cartesian coordinates for all asio2 atoms
            cart_coords_asio2 = np.array(asio2.cart_coords)

            indices_to_add = [
                i for i, coord in enumerate(cart_coords_asio2)
                if all(lower_bounds_asio2 <= coord) and all(coord <= upper_bounds_asio2)
            ]

            coords_to_add_asio2_frame = [cart_coords_asio2[i] for i in indices_to_add]
            coords_to_add = [[x,y,z - asio2.lattice.c + si.lattice.c] for x,y,z in coords_to_add_asio2_frame]
            species_to_add = [asio2.species[i] for i in indices_to_add]

            for i in range(len(coords_to_add)):

                # Compute distances to existing Si atoms
                distances = cdist([coords_to_add[i]], np.array(si.cart_coords))
                if len(distances) > 0:
                    min_distance = distances.min()

                    # Determine threshold based on atom type
                    if species_to_add[i].symbol == "Si":
                        threshold = si_si_threshold
                    elif species_to_add[i].symbol == "O":
                        threshold = si_o_threshold
                    else:
                        print(species_to_add[i].symbol)

                    # Skip appending if too close to existing atoms
                    if min_distance < threshold:
                        print(min_distance)
                        continue

                si.append(species_to_add[i], coords_to_add[i], coords_are_cartesian=True)
    
    if si_indent == False:
        lower_bounds = [x_int * unit_cell_length - delta, y_int * unit_cell_length - delta, - delta]
        upper_bounds = [x_int * unit_cell_length + unit_cell_length - delta, y_int * unit_cell_length + unit_cell_length - delta, unit_cell_length - delta]


        # Get Cartesian coordinates for all si atoms
        cart_coords_asio2 = np.array(asio2.cart_coords)

        # Identify indices of atoms within the specified region
        indices_to_remove = [
            i for i, coord in enumerate(cart_coords_asio2)
            if all(lower_bounds <= coord) and all(coord <= upper_bounds)
        ]


        # Remove the identified atoms
        asio2.remove_sites(indices_to_remove)

        if fill_indent == True:

            lower_bounds_si = [x_int * unit_cell_length - delta, y_int * unit_cell_length - delta, - delta]
            upper_bounds_si = [x_int * unit_cell_length + unit_cell_length - delta, y_int * unit_cell_length + unit_cell_length - delta, unit_cell_length - delta]

            # Get Cartesian coordinates for all asio2 atoms
            cart_coords_si = np.array(si.cart_coords)

            indices_to_add = [
                i for i, coord in enumerate(cart_coords_si)
                if all(lower_bounds_si <= coord) and all(coord <= upper_bounds_si)
            ]

            coords_to_add_si_frame = [cart_coords_si[i] for i in indices_to_add]
            coords_to_add = [[x,y,z] for x,y,z in coords_to_add_si_frame] # no change since both start at 0 (if I use vacuum this needs to be changed)
            species_to_add = [si.species[i] for i in indices_to_add] # Always Si

            for i in range(len(coords_to_add)):

                coord_to_add = coords_to_add[i]
                species_to_insert = species_to_add[i]
                print(coord_to_add)
                print(species_to_insert)

                # Append the Si atom to asio2
                asio2.append(species_to_insert, coord_to_add, coords_are_cartesian=True)

                # Compute distances between the newly added Si atom and existing asio2 atoms
                new_atom_index = len(asio2) - 1
                distances = cdist([asio2[new_atom_index].coords], np.array(asio2.cart_coords[:-1]))

                # Find indices of atoms to remove based on thresholds
                indices_to_remove = []
                for j, distance in enumerate(distances[0]):
                    existing_species = asio2.species[j]

                    if existing_species.symbol == "Si" and distance < si_si_threshold:
                        indices_to_remove.append(j)
                    elif existing_species.symbol == "O" and distance < si_o_threshold:
                        indices_to_remove.append(j)

                # Remove overlapping atoms
                if indices_to_remove:
                    asio2.remove_sites(indices_to_remove)

    return si, asio2

def main_generator(miller_index, min_thickness, filename_asio2, filename_interface, asio2_shift = [0,0,0]):
    asio2 = fetch_asio2(filename_asio2, asio2_shift = asio2_shift)
    si_slab = generate_si_slab(miller_index, min_thickness, asio2)
    generate_interface(si_slab,asio2,filename_interface)
