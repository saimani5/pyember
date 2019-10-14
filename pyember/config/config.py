"""The configurations class should have

1. 3D Array
2. List of atoms and coordinates
3. Conversion between 1 and 2
4. Periodic boundary conditions
5. Lattice Type
6. Atom Names
7. Helper functions (possibly)
8. Easier slicing methods
"""

import numpy as np


class config():


    def __init__(self):
        self.bc = boundary_conditions

    def array_to_list(self, input_array):
        """Inputs: Array of size X x Y x Z, where the values represent the type of atom (0 = no atom),
        this function converts this array to list of atomic positions with atom type"""

        list_array = []
        for i in range(input_array.shape[0]):
            for j in range(input_array.shape[1]):
                for k in range(input_array.shape[2]):
                    if input_array[i,j,k]!=0: list_array.append((i,j,k,input_array[i,j,k]))

        return list_array

    def list_to_array(self, list_array):
        """Inputs: list of atomic positions, where the last column represent the type of atom.
       This function converts this list to an [x,y,z] box"""

        list_np_array = np.array(list_array)

        # Find number of steps required
        min_x, min_y, min_z = [np.min(list_np_array[:, 0]), np.min(list_np_array[:, 1]), np.min(list_np_array[:, 2])]
        max_x, max_y, max_z = [np.max(list_np_array[:, 0]), np.max(list_np_array[:, 1]), np.max(list_np_array[:, 2])]

        # Find number of steps required, based off of unique entries in each axis
        num_x_steps = (len(np.unique(list_np_array[:, 0])))
        num_y_steps = (len(np.unique(list_np_array[:, 1])))
        num_z_steps = (len(np.unique(list_np_array[:, 2])))

        # Find multiplication factors along each axis required to convert xyz into integers
        factors = [(num_x_steps - 1) / (max_x - min_x), (num_y_steps - 1) / (max_y - min_y),
                   (num_z_steps - 1) / (max_z - min_z)]

        atoms_array = np.zeros(shape=(num_x_steps, num_y_steps, num_z_steps))

        for entry in list_array:
            x, y, z, c = entry
            #Find multiplication factors along each axis requried to convert xyz into integers
            factors = [(num_x_steps - 1) / (max_x - min_x), (num_y_steps - 1) / (max_y - min_y),
                       (num_z_steps - 1) / (max_z - min_z)]
            atoms_array[int(x * factors[0]), int(y * factors[1]), int(z * factors[2])] = c

        return atoms_array


    def make_lattice(self, latt_type = 'cubic', lat_parms):
        """This function will make and return a lattice of type latt_type, made
        with lattice parameters given in lat_parms"""

        if latt_type = 'cubic':
            lx, ly, lz = lat_parms
            latt = {}
            latt['box'] = ['cubic', lx, ly, lz]
            latt['xyzs'] = []

            # box dimensions in lattice units

            # layer number
            for iz in range(lz):
                # layer structure
                for iy in range(ly):
                    for ix in range(lx):
                        latt['xyzs'].append([ix, iy, iz,1])

        elif latt_type = 'bcc':
            lx, ly, lz = lat_parms
            latt = {}
            latt['box'] = ['bcc', lx, ly, lz]
            latt['xyzs'] = []

            # box dimensions in lattice units

            # layer number
            for iz in range(lz):
                # layer structure
                for iy in range(ly):
                    for ix in range(lx):
                        if ix + 0.5 <= (lx - 1) and iy + 0.5 <= (ly - 1) and iz + 0.5 <= (lz - 1):
                            latt['xyzs'].append([ix + 0.5, iy + 0.5, iz + 0.5, 1])
                        latt['xyzs'].append([1 * ix, 1 * iy, 1 * iz, 1])



        elif latt_type = 'fcc':
            lx, ly, lz = lat_parms

            latt = {}
            latt['nat'] = lx * ly * lz
            latt['box'] = ['fcc', 2 * lx, ly, lz]
            latt['xyzs'] = []

            # box dimensions in lattice units

            # layer number
            for iz in range(lz):
                # layer structure
                for iy in range(ly):
                    for ix in range(lx):
                        rx = 2 * ix + (iy + iz) % 2
                        latt['xyzs'].append([rx, iy, iz,1])

        return latt