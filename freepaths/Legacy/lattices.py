from random import random
from numpy import sign, zeros
from math import pi, cos, sin, sqrt, asin
import matplotlib.pyplot as plt

from parameters import *


def hole_positioning(hole_lattice_type, rectangular_hole_side_y, rectangular_hole_side_x, width, period_x, period_y):
    """This function places holes in space, depending on the lattice type,
    and returns coordinates of their centers and changes in their size (if any).
    In the hole_coordinates array, column 0 is X, column 1 is Y, column 2 is correction of the size,
    i.e. 0 is no correction, 1 means the hole will be +100% larger, -0.5 means -50% smaller.
    hole_shapes list contains shapes of hole e.g. 'rectangle', 'circle', 'triangle_up' or 'triangle_down' """

    if hole_lattice_type == 'yana_pnc':
        '''Regular square lattice of holes with secondary periodicity'''

        # Parameters:
        first_hole_coordinate = 500e-9
        number_of_periods_x = 3
        number_of_periods_y = 4
        number_of_subperiods_x = 3
        number_of_subperiods_y = 2
        space_x = 300e-9
        space_y = 300e-9

        # Filling the lattice (do not change):
        hole_coordinates = zeros((number_of_periods_x * number_of_periods_y * number_of_subperiods_x * number_of_subperiods_y
                                + number_of_subperiods_x*number_of_subperiods_y*(number_of_periods_y//2), 3))
        hole_shapes = ['circle' for x in range(hole_coordinates.shape[0])]
        hole_number = 0
        shift_y = 0
        global_shift_x = (number_of_periods_x*number_of_subperiods_x*period_x+space_x)/2
        for i in range(number_of_periods_y):
            shift_x = -(i%2)*(number_of_subperiods_x*period_x + space_x)/2 - global_shift_x
            for j in range(number_of_periods_x + i%2):
                for sub_i in range(number_of_subperiods_y):
                    for sub_j in range(number_of_subperiods_x):
                        hole_coordinates[hole_number, 0] = shift_x + sub_j * period_x
                        hole_coordinates[hole_number, 1] = shift_y + first_hole_coordinate + sub_i * period_y
                        hole_number += 1
                shift_x += number_of_periods_x*period_x + space_x
            shift_y += number_of_periods_y*period_y + space_y

    if hole_lattice_type == 'directional_source':
        '''Directional source without any holes'''
        hole_coordinates = zeros((1, 3))
        hole_shapes = ['rectangle' for x in range(hole_coordinates.shape[0])]

    if hole_lattice_type == 'probe':
        '''Point source without any holes'''
        hole_coordinates = zeros((1, 3))
        hole_shapes = ['rectangle' for x in range(hole_coordinates.shape[0])]

    if hole_lattice_type == 'square':
        '''Regular square lattice of holes'''
        first_hole_coordinate = 300e-9
        number_of_periods_x = 1
        number_of_periods_y = 7
        hole_coordinates = zeros((number_of_periods_x * number_of_periods_y, 3))
        hole_shapes = ['circle' for x in range(hole_coordinates.shape[0])]
        # hole_shapes = ['rectangle' for x in range(hole_coordinates.shape[0])]
        # hole_shapes = ['triangle_up' for x in range(hole_coordinates.shape[0])]
        hole_number = 0
        for i in range(number_of_periods_y):
            for j in range(number_of_periods_x):
                hole_coordinates[hole_number, 0] = -(number_of_periods_x - 1) * period_x / 2 + j * period_x
                hole_coordinates[hole_number, 1] = first_hole_coordinate + i * period_y
                hole_number += 1

    if hole_lattice_type == 'diamond_particle':
        '''Diamond shaped particle'''
        number_of_periods_x = 2
        number_of_periods_y = 2
        hole_coordinates = zeros((number_of_periods_x * number_of_periods_y, 3))
        hole_coordinates = zeros((number_of_periods_x * number_of_periods_y, 3))
        hole_shapes = ['triangle_up' for x in range(hole_coordinates.shape[0])]
        hole_number = 0
        for i in range(number_of_periods_y):
            for j in range(number_of_periods_x):
                hole_coordinates[hole_number, 0] = -(number_of_periods_x - 1) * period_x / 2 + j * period_x
                hole_coordinates[hole_number, 1] = period_y/2 + i * period_y
                hole_shapes[hole_number] = 'triangle_up' if i == 0 else 'triangle_down'
                hole_number += 1

    if hole_lattice_type == 'source':
        '''Lattice similar to staggered, but with a passage in the middle'''
        first_hole_coordinate = 300e-9
        number_of_periods_x = 5
        number_of_periods_y = 3
        hole_coordinates = zeros( (number_of_periods_x * number_of_periods_y * 4 + number_of_periods_y * 2, 3))
        hole_shapes = ['circle' for x in range(hole_coordinates.shape[0])]
        hole_number = 0
        for i in range(number_of_periods_y):
            for j in range(number_of_periods_x):
                # Aligned rows
                hole_coordinates[hole_number, 0] = j * period_x + period_x / 2.0
                hole_coordinates[hole_number, 1] = first_hole_coordinate + i * period_y * 2
                hole_coordinates[hole_number + 1, 0] = -j * period_x - period_x / 2.0
                hole_coordinates[hole_number + 1, 1] = first_hole_coordinate + i * period_y * 2
                # Staggered rows
                hole_coordinates[hole_number + 2, 0] = j * period_x + period_x
                hole_coordinates[ hole_number + 2, 1] = first_hole_coordinate + i * period_y * 2.0 + period_y
                hole_coordinates[hole_number + 3, 0] = -j * period_x - period_x
                hole_coordinates[ hole_number + 3, 1] = first_hole_coordinate + i * period_y * 2.0 + period_y
                hole_number += 4
        for i in range(number_of_periods_y):
            # Additional holes
            hole_coordinates[hole_number, 0] = period_x / 2
            hole_coordinates[ hole_number, 1] = first_hole_coordinate + i * period_y * 2.0 + period_y
            hole_coordinates[hole_number + 1, 0] = -period_x / 2
            hole_coordinates[ hole_number + 1, 1] = first_hole_coordinate + i * period_y * 2.0 + period_y
            hole_number += 2

    if hole_lattice_type == 'sources':
        '''Lattice similar to staggered, but with several passages'''
        first_hole_coordinate = 300e-9
        number_of_periods_x = 3
        number_of_periods_y = 3
        number_of_sources = 2
        hole_coordinates = zeros(
            (number_of_sources *
             (number_of_periods_x * number_of_periods_y * 4 +
              number_of_periods_y * 4), 3))
        hole_shapes = ['circle' for x in range(hole_coordinates.shape[0])]
        hole_number = 0
        for k in range(number_of_sources):
            for i in range(number_of_periods_y):
                for j in range(number_of_periods_x):
                    # Aligned rows
                    hole_coordinates[
                        hole_number, 0] = j * period_x + period_x / 2.0 + ( number_of_periods_x + 0.5) * k * period_x
                    hole_coordinates[ hole_number, 1] = first_hole_coordinate + i * period_y * 2
                    hole_coordinates[ hole_number + 1, 0] = -j * period_x - period_x / 2.0 - ( number_of_periods_x + 0.5) * k * period_x
                    hole_coordinates[ hole_number + 1, 1] = first_hole_coordinate + i * period_y * 2
                    # Staggered rows
                    hole_coordinates[ hole_number + 2, 0] = j * period_x + period_x + (number_of_periods_x + 0.5) * k * period_x
                    hole_coordinates[ hole_number + 2, 1] = first_hole_coordinate + i * period_y * 2.0 + period_y
                    hole_coordinates[ hole_number + 3, 0] = -j * period_x - period_x - (number_of_periods_x + 0.5) * k * period_x
                    hole_coordinates[ hole_number + 3, 1] = first_hole_coordinate + i * period_y * 2.0 + period_y
                    hole_number += 4
        for k in range(number_of_sources):
            for i in range(number_of_periods_y):
                # Additional holes
                hole_coordinates[hole_number, 0] = period_x / 2 + (number_of_periods_x + 0.5) * k * period_x
                hole_coordinates[ hole_number, 1] = first_hole_coordinate + i * period_y * 2.0 + period_y
                hole_coordinates[hole_number + 1, 0] = -period_x / 2 - (number_of_periods_x + 0.5) * k * period_x
                hole_coordinates[ hole_number + 1, 1] = first_hole_coordinate + i * period_y * 2.0 + period_y
                hole_coordinates[hole_number + 2, 0] = -number_of_periods_x * period_x - ( number_of_periods_x + 0.5) * k * period_x
                hole_coordinates[ hole_number + 2, 1] = first_hole_coordinate + i * period_y * 2.0
                hole_coordinates[hole_number + 3, 0] = number_of_periods_x * period_x + ( number_of_periods_x + 0.5) * k * period_x
                hole_coordinates[ hole_number + 3, 1] = first_hole_coordinate + i * period_y * 2.0
                hole_number += 4

    if hole_lattice_type == 'source_with_mirror':
        '''Lattice similar to staggered, but with a passage in the middle plus a mirror'''
        first_hole_coordinate = 300e-9
        #passage=100e-9
        number_of_periods_x = 5
        number_of_periods_y = 3
        hole_coordinates = zeros((number_of_periods_x * number_of_periods_y * 4 + number_of_periods_y * 2 + 1, 3))
        hole_shapes = ['circle' for x in range(hole_coordinates.shape[0])]
        hole_number = 0
        for i in range(number_of_periods_y):
            for j in range(number_of_periods_x):
                # Aligned rows
                hole_coordinates[hole_number, 0] = j * period_x + period_x / 2.0
                hole_coordinates[hole_number, 1] = first_hole_coordinate + i * period_y * 2
                hole_coordinates[hole_number + 1, 0] = -j * period_x - period_x / 2.0
                hole_coordinates[hole_number + 1, 1] = first_hole_coordinate + i * period_y * 2
                # Staggered rows
                hole_coordinates[hole_number + 2, 0] = j * period_x + period_x
                hole_coordinates[ hole_number + 2, 1] = first_hole_coordinate + i * period_y * 2.0 + period_y
                hole_coordinates[hole_number + 3, 0] = -j * period_x - period_x
                hole_coordinates[ hole_number + 3, 1] = first_hole_coordinate + i * period_y * 2.0 + period_y
                hole_number += 4

        for i in range(number_of_periods_y):
            # Additional holes
            hole_coordinates[hole_number, 0] = period_x / 2
            hole_coordinates[ hole_number, 1] = first_hole_coordinate + i * period_y * 2.0 + period_y
            hole_coordinates[hole_number + 1, 0] = -period_x / 2
            hole_coordinates[ hole_number + 1, 1] = first_hole_coordinate + i * period_y * 2.0 + period_y
            hole_number += 2

        # Mirror triangle
        hole_coordinates[hole_number, 0] = -600e-9
        hole_coordinates[hole_number, 1] = 3000e-9 - rectangular_hole_side_y / 2
        hole_shapes[hole_number] = 'triangle_down'

    if hole_lattice_type == 'slits':
        '''This lattece is array of slits from PRB 101, 115301 (2020)'''
        first_hole_coordinate = 500e-9 - 75e-9
        number_of_periods_x = 13
        number_of_periods_y = 4
        hole_coordinates = zeros((number_of_periods_x * number_of_periods_y, 3))
        hole_shapes = ['rectangle' for x in range(hole_coordinates.shape[0])]
        hole_number = 0
        for i in range(number_of_periods_y):
            for j in range(number_of_periods_x):
                hole_coordinates[hole_number, 0] = -(number_of_periods_x - 1) * period_x / 2 + j * period_x
                hole_coordinates[hole_number, 1] = first_hole_coordinate + i * period_y
                hole_number += 1

    if hole_lattice_type == 'tesla valve forward':
        hole_coordinates = zeros((13, 3))
        hole_shapes = []
        first_hole_coordinate = 500e-9
        repetitions = 4
        # Big triangles
        for i in range(repetitions):
            hole_coordinates[i, 0] = ((-1)**(i + 1)) * period_x / 2
            hole_coordinates[i, 1] = i * period_y / 2 + first_hole_coordinate
            hole_shapes.append('triangle_down')

        # Small triangles
        shift = 250e-9
        shift_y = 270e-9
        for i, j in enumerate(range(repetitions, repetitions * 2)):
            hole_coordinates[j, 0] = ((-1)**i) * period_x / 2 + ((-1)**(i + 1)) * shift
            hole_coordinates[j, 1] = first_hole_coordinate + i * 0.5 * period_y + shift_y
            hole_coordinates[j, 2] = -0.6
            hole_shapes.append('triangle_up')

        # Additional small triangle at the end
        hole_coordinates[2 * repetitions, 0] = -period_x / 2 - 10e-9
        hole_coordinates[2 * repetitions, 1] = 1.8 * period_y + first_hole_coordinate - 10e-9
        hole_coordinates[2 * repetitions, 2] = -0.4
        hole_shapes.append('triangle_down')

        # Lower rectangles
        hole_coordinates[2 * repetitions + 1, 0] = -3 * period_x / 4 - 20e-9
        hole_coordinates[2 * repetitions + 1, 1] = 0
        hole_shapes.append('rectangle')
        hole_coordinates[2 * repetitions + 2, 0] = +2 * period_x / 4 - 20e-9
        hole_coordinates[2 * repetitions + 2, 1] = 0
        hole_shapes.append('rectangle')

        # Upper rectangles
        hole_coordinates[2 * repetitions + 3, 0] = -3 * period_x / 4
        hole_coordinates[2 * repetitions + 3, 1] = period_y * 3 + 100e-9
        hole_shapes.append('rectangle')
        hole_coordinates[2 * repetitions + 4, 0] = +2 * period_x / 4
        hole_coordinates[2 * repetitions + 4, 1] = period_y * 3
        hole_shapes.append('rectangle')

    if hole_lattice_type == 'tesla valve backward':
        hole_coordinates = zeros((13, 3))
        hole_shapes = []
        first_hole_coordinate = 500e-9
        repetitions = 4
        # Big triangles
        for i in range(repetitions):
            hole_coordinates[i, 0] = ((-1) ** i) * period_x / 2
            hole_coordinates[
                i, 1] = (i + 1) * period_y / 2 + first_hole_coordinate
            hole_shapes.append('triangle_up')

        # Small triangles
        shift = 250e-9
        shift_y = 230e-9
        for i, j in enumerate(range(repetitions, repetitions * 2)):
            hole_coordinates[j, 0] = ((-1)**(i + 1)) * period_x / 2 + ((-1) ** i) * shift
            hole_coordinates[j, 1] = first_hole_coordinate + i * 0.5 * period_y + shift_y
            hole_coordinates[j, 2] = -0.6
            hole_shapes.append('triangle_down')

        # Additional small triangle at the end
        hole_coordinates[2 * repetitions, 0] = -period_x / 2 - 45e-9
        hole_coordinates[2 * repetitions, 1] = 1 * period_y / 2 + first_hole_coordinate - 220e-9
        hole_coordinates[2 * repetitions, 2] = -0.35
        hole_shapes.append('triangle_up')

        # Lower rectangles
        hole_coordinates[2 * repetitions + 1, 0] = -3 * period_x / 4 - 5e-9
        hole_coordinates[2 * repetitions + 1, 1] = 0 - 55e-9
        hole_shapes.append('rectangle')
        hole_coordinates[2 * repetitions + 2, 0] = +2 * period_x / 4 - 5e-9
        hole_coordinates[2 * repetitions + 2, 1] = 0 - 20e-9
        hole_shapes.append('rectangle')

        # Upper rectangles
        hole_coordinates[2 * repetitions + 3, 0] = -3 * period_x / 4 - 20e-9
        hole_coordinates[2 * repetitions + 3, 1] = period_y * 3
        hole_shapes.append('rectangle')
        hole_coordinates[2 * repetitions + 4, 0] = +2 * period_x / 4 - 20e-9
        hole_coordinates[2 * repetitions + 4, 1] = period_y * 3
        hole_shapes.append('rectangle')

    if hole_lattice_type == 'slit':
        '''Single slit from Hao et al.'''
        first_hole_coordinate = 500e-9
        number_of_periods_x = 2
        number_of_periods_y = 1
        hole_coordinates = zeros((number_of_periods_x * number_of_periods_y, 3))
        hole_shapes = ['rectangle' for x in range(hole_coordinates.shape[0])]
        hole_number = 0
        for i in range(number_of_periods_y):
            for j in range(number_of_periods_x):
                hole_coordinates[hole_number, 0] = -(number_of_periods_x - 1) * period_x / 2 + j * period_x
                hole_coordinates[hole_number, 1] = first_hole_coordinate + i * period_y
                hole_number += 1

    if hole_lattice_type == 'fishbone':
        '''Fishbone lattice, similar to the one in Maire et al. Scientific Report 8, 4452 (2018)'''
        first_hole_coordinate = 300e-9
        number_of_periods_x = 2
        number_of_periods_y = 4
        hole_coordinates = zeros((number_of_periods_x * number_of_periods_y, 3))
        hole_shapes = ['rectangle' for x in range(hole_coordinates.shape[0])]
        hole_number = 0
        for i in range(number_of_periods_y):
            for j in range(number_of_periods_x):
                hole_coordinates[hole_number, 0] = -(number_of_periods_x - 1) * period_x / 2 + j * period_x
                hole_coordinates[hole_number, 1] = first_hole_coordinate + i * period_y
                hole_number += 1

    if hole_lattice_type == 'staggered':
        first_hole_coordinate = 300e-9
        number_of_periods_x = 5
        number_of_periods_y = 3
        hole_coordinates = zeros((number_of_periods_x * number_of_periods_y * 2 + number_of_periods_y, 3))
        hole_shapes = ['circle' for x in range(hole_coordinates.shape[0])]
        hole_number = 0
        for i in range(number_of_periods_y):
            for j in range(number_of_periods_x):
                hole_coordinates[hole_number, 0] = -(number_of_periods_x - 1) * period_x / 2 + j * period_x
                hole_coordinates[hole_number, 1] = first_hole_coordinate + i * period_y * 2
                hole_number += 1
        for i in range(number_of_periods_y):
            for j in range(number_of_periods_x + 1):
                hole_coordinates[ hole_number, 0] = -(number_of_periods_x - 1) * period_x / 2 + j * period_x - period_x / 2
                hole_coordinates[ hole_number, 1] = first_hole_coordinate + period_y + i * period_y * 2
                hole_number += 1

    if hole_lattice_type == 'staggered_triangles_down':
        first_hole_coordinate = 600e-9
        number_of_periods_x = 5
        number_of_periods_y = 3
        hole_coordinates = zeros((number_of_periods_x * number_of_periods_y * 2 + number_of_periods_y, 3))
        hole_shapes = ['triangle_down' for x in range(hole_coordinates.shape[0])]
        hole_number = 0
        for i in range(number_of_periods_y):
            for j in range(number_of_periods_x):
                hole_coordinates[hole_number, 0] = -(number_of_periods_x - 1) * period_x / 2 + j * period_x
                hole_coordinates[ hole_number, 1] = first_hole_coordinate + i * period_y * 1.0
                hole_number += 1
        for i in range(number_of_periods_y):
            for j in range(number_of_periods_x + 1):
                hole_coordinates[ hole_number, 0] = -(number_of_periods_x - 1) * period_x / 2 + j * period_x - period_x / 2
                hole_coordinates[ hole_number, 1] = first_hole_coordinate + period_y * 0.5 + i * period_y * 1.0
                hole_number += 1

    if hole_lattice_type == 'staggered_triangles_up':
        first_hole_coordinate = 600e-9
        number_of_periods_x = 5
        number_of_periods_y = 3
        hole_coordinates = zeros((number_of_periods_x * number_of_periods_y * 2 + number_of_periods_y, 3))
        hole_shapes = ['triangle_up' for x in range(hole_coordinates.shape[0])]
        hole_number = 0
        for i in range(number_of_periods_y):
            for j in range(number_of_periods_x):
                hole_coordinates[hole_number, 0] = -(number_of_periods_x - 1) * period_x / 2 + j * period_x
                hole_coordinates[ hole_number, 1] = first_hole_coordinate + i * period_y * 1.0
                hole_number += 1
        for i in range(number_of_periods_y):
            for j in range(number_of_periods_x + 1):
                hole_coordinates[ hole_number, 0] = -(number_of_periods_x - 1) * period_x / 2 + j * period_x - period_x / 2
                hole_coordinates[ hole_number, 1] = first_hole_coordinate + period_y * 0.5 + i * period_y * 1.0
                hole_number += 1

    if hole_lattice_type == 'serpentine':
        hole_coordinates = zeros((7, 3))
        neck = 155e-9
        hole_coordinates[0, 0] = neck / 2
        hole_coordinates[0, 1] = 0.0
        hole_coordinates[1, 0] = neck / 2
        hole_coordinates[1, 1] = (rectangular_hole_side_y + 1 * neck) / 4
        hole_coordinates[2, 0] = neck / 2
        hole_coordinates[2, 1] = (rectangular_hole_side_y + 2.5 * neck) / 4
        for i in range(1, 5):
            hole_coordinates[i + 2, 0] = sign(-0.5 + i % 2) * neck / 2
            hole_coordinates[i + 2, 1] = (2 * i - 1) * (rectangular_hole_side_y) / 2 + neck * i
        hole_shapes = ['rectangle' for x in range(hole_coordinates.shape[0])]

    if hole_lattice_type == 'diode_with_wires':
        first_hole_coordinate = period_y
        number_of_periods_x = 2
        number_of_periods_y = 3
        hole_coordinates = zeros((number_of_periods_x * number_of_periods_y +
                                  number_of_periods_y + 4, 3))
        hole_shapes = ['triangle_down' for x in range(hole_coordinates.shape[0])]
        hole_number = 0
        for i in range(number_of_periods_y):  # These are large holes
            for j in range(number_of_periods_x):
                hole_coordinates[hole_number,
                                 0] = -(number_of_periods_x -
                                        1) * period_x / 2 + j * period_x
                hole_coordinates[hole_number,
                                 1] = first_hole_coordinate + i * period_y
                hole_coordinates[hole_number, 2] = 0.0
                hole_number += 1
        for i in range(number_of_periods_y):  # These are small holes
            hole_coordinates[hole_number, 0] = 0.0
            hole_coordinates[
                hole_number,
                1] = (1.1 * period_y) + i * period_y  # first_hole_coordinate
            hole_coordinates[hole_number, 2] = -0.6  #-0.5625
            hole_number += 1
        # These are rectangles that form wires at the beginning and the end
        hole_coordinates[hole_number, 0] = -width / 2
        hole_coordinates[hole_number, 1] = 0.0
        hole_shapes[hole_number] = 'rectangle'
        hole_coordinates[hole_number + 1, 0] = +width / 2
        hole_coordinates[hole_number + 1, 1] = 0.0
        hole_shapes[hole_number + 1] = 'rectangle'
        hole_coordinates[hole_number + 2, 0] = -width / 2
        hole_coordinates[hole_number + 2, 1] = period_y * 4
        hole_shapes[hole_number + 2] = 'rectangle'
        hole_coordinates[hole_number + 3, 0] = +width / 2
        hole_coordinates[hole_number + 3, 1] = period_y * 4
        hole_shapes[hole_number + 3] = 'rectangle'

    if hole_lattice_type == 'cloak':
        first_hole_coordinate = 200e-9
        total_number_of_holes = 38
        hole_in_the_center = False
        total_number_of_holes += hole_in_the_center * 1  # i.e. if there is hole in the center, than add a hole
        hole_coordinates = zeros((total_number_of_holes, 3))
        hole_shapes = ['circle' for x in range(hole_coordinates.shape[0])]
        hole_coordinates[0, 0] = 0.0
        hole_coordinates[0, 1] = first_hole_coordinate
        hole_coordinates[1, 0] = -period_x / 2
        hole_coordinates[1, 1] = first_hole_coordinate + period_y
        hole_coordinates[2, 0] = +period_x / 2
        hole_coordinates[2, 1] = first_hole_coordinate + period_y
        for i in range(3):
            hole_coordinates[3 + i, 0] = -period_x + period_x * i
            hole_coordinates[3 + i, 1] = first_hole_coordinate + period_y * 2
        for i in range(4):
            hole_coordinates[6 + i, 0] = -1.5 * period_x + period_x * i
            hole_coordinates[6 + i, 1] = first_hole_coordinate + period_y * 3
        for i in range(5):
            hole_coordinates[10 + i, 0] = -2.0 * period_x + period_x * i
            hole_coordinates[10 + i, 1] = first_hole_coordinate + period_y * 4
        hole_coordinates[15, 0] = -2.25 * period_x
        hole_coordinates[15, 1] = first_hole_coordinate + period_y * 5
        hole_coordinates[16, 0] = +2.25 * period_x
        hole_coordinates[16, 1] = first_hole_coordinate + period_y * 5
        for i in range(2):
            hole_coordinates[17 + 2 * i, 0] = -2.4 * period_x
            hole_coordinates[17 + 2 * i,
                             1] = first_hole_coordinate + period_y * (6 + i)
            hole_coordinates[18 + 2 * i, 0] = +2.4 * period_x
            hole_coordinates[18 + 2 * i,
                             1] = first_hole_coordinate + period_y * (6 + i)
        hole_coordinates[21, 0] = -2.25 * period_x
        hole_coordinates[21, 1] = first_hole_coordinate + period_y * 8
        hole_coordinates[22, 0] = +2.25 * period_x
        hole_coordinates[22, 1] = first_hole_coordinate + period_y * 8
        for i in range(5):
            hole_coordinates[23 + i, 0] = -2.0 * period_x + period_x * i
            hole_coordinates[23 + i, 1] = first_hole_coordinate + period_y * 9
        for i in range(4):
            hole_coordinates[28 + i, 0] = -1.5 * period_x + period_x * i
            hole_coordinates[28 + i, 1] = first_hole_coordinate + period_y * 10
        for i in range(3):
            hole_coordinates[32 + i, 0] = -period_x + period_x * i
            hole_coordinates[32 + i, 1] = first_hole_coordinate + period_y * 11
        hole_coordinates[35, 0] = -period_x / 2
        hole_coordinates[35, 1] = first_hole_coordinate + period_y * 12
        hole_coordinates[36, 0] = period_x / 2
        hole_coordinates[36, 1] = first_hole_coordinate + period_y * 12
        hole_coordinates[37, 0] = 0.0
        hole_coordinates[37, 1] = first_hole_coordinate + period_y * 13

        # This is the big hole in the center
        if hole_in_the_center:
            hole_coordinates[38, 0] = 0
            hole_coordinates[38, 1] = first_hole_coordinate + period_y * 6.5
            hole_coordinates[38, 2] = 3.5

    if hole_lattice_type == 'cloak_2':
        first_hole_coordinate = 300e-9
        period_x = 300e-9
        period_y = (sqrt(3) / 2) * period_x
        total_number_of_holes = 54
        hole_in_the_center = False
        total_number_of_holes += hole_in_the_center * 1  # i.e. if there is hole in the center, than add a hole
        hole_coordinates = zeros((total_number_of_holes, 3))
        hole_shapes = ['circle' for x in range(hole_coordinates.shape[0])]
        hole_coordinates[0, 0] = 0.0
        hole_coordinates[0, 1] = first_hole_coordinate
        hole_coordinates[1, 0] = -period_x / 2
        hole_coordinates[1, 1] = first_hole_coordinate + period_y
        hole_coordinates[2, 0] = +period_x / 2
        hole_coordinates[2, 1] = first_hole_coordinate + period_y
        for i in range(3):
            hole_coordinates[3 + i, 0] = -period_x + period_x * i
            hole_coordinates[3 + i, 1] = first_hole_coordinate + period_y * 2
        for i in range(4):
            hole_coordinates[6 + i, 0] = -1.5 * period_x + period_x * i
            hole_coordinates[6 + i, 1] = first_hole_coordinate + period_y * 3
        for i in range(5):
            hole_coordinates[10 + i, 0] = -2.0 * period_x + period_x * i
            hole_coordinates[10 + i, 1] = first_hole_coordinate + period_y * 4
        hole_coordinates[15, 0] = -2.25 * period_x
        hole_coordinates[15, 1] = first_hole_coordinate + period_y * 5
        hole_coordinates[16, 0] = +2.25 * period_x
        hole_coordinates[16, 1] = first_hole_coordinate + period_y * 5
        for i in range(2):
            hole_coordinates[17 + 2 * i, 0] = -2.4 * period_x
            hole_coordinates[17 + 2 * i, 1] = first_hole_coordinate + period_y * (6 + i)
            hole_coordinates[18 + 2 * i, 0] = +2.4 * period_x
            hole_coordinates[18 + 2 * i, 1] = first_hole_coordinate + period_y * (6 + i)
        hole_coordinates[21, 0] = -2.25 * period_x
        hole_coordinates[21, 1] = first_hole_coordinate + period_y * 8
        hole_coordinates[22, 0] = +2.25 * period_x
        hole_coordinates[22, 1] = first_hole_coordinate + period_y * 8
        for i in range(5):
            hole_coordinates[23 + i, 0] = -2.0 * period_x + period_x * i
            hole_coordinates[23 + i, 1] = first_hole_coordinate + period_y * 9

        # Holes at the sides
        for i in range(9):
            hole_coordinates[28 + i, 0] = 3.35 * period_x
            hole_coordinates[28 + i, 1] = first_hole_coordinate + period_x * 0.1 + period_x * i
        for i in range(9):
            hole_coordinates[37 + i, 0] = -3.35 * period_x
            hole_coordinates[37 + i, 1] = first_hole_coordinate + period_x * 0.1 + period_x * i

        # Additional holes right
        hole_coordinates[46, 0] = period_x * sqrt(3) / 2 + 1.1 * period_x / 2
        hole_coordinates[46, 1] = first_hole_coordinate + period_y / 2
        hole_coordinates[46, 2] = -0.2
        hole_coordinates[47, 0] = period_x * sqrt(3) / 2 + 2.5 * period_x / 2
        hole_coordinates[47, 1] = first_hole_coordinate + 0.65 * period_y / 2 + period_x * sqrt(3) / 2
        hole_coordinates[47, 2] = 0.35
        hole_coordinates[48, 0] = period_x * sqrt(3) / 2 + 2.9 * period_x / 2
        hole_coordinates[48, 1] = first_hole_coordinate + period_y / 2 + 2 * period_x * sqrt(3) / 2
        hole_coordinates[48, 2] = -0.2
        # Additional holes left
        hole_coordinates[49, 0] = -period_x * sqrt(3) / 2 - 1.1 * period_x / 2
        hole_coordinates[49, 1] = first_hole_coordinate + period_y / 2
        hole_coordinates[49, 2] = -0.2
        hole_coordinates[50, 0] = -period_x * sqrt(3) / 2 - 2.5 * period_x / 2
        hole_coordinates[50, 1] = first_hole_coordinate + 0.65 * period_y / 2 + period_x * sqrt(3) / 2
        hole_coordinates[50, 2] = 0.35
        hole_coordinates[51, 0] = -period_x * sqrt(3) / 2 - 2.9 * period_x / 2
        hole_coordinates[51, 1] = first_hole_coordinate + period_y / 2 + 2 * period_x * sqrt(3) / 2
        hole_coordinates[51, 2] = -0.2

        # Small holes at the beginning
        hole_coordinates[52, 0] = 3.35 * period_x - period_x
        hole_coordinates[52, 1] = first_hole_coordinate + period_x * 0.1
        hole_coordinates[52, 2] = -0.15
        hole_coordinates[53, 0] = -3.35 * period_x + period_x
        hole_coordinates[53, 1] = first_hole_coordinate + period_x * 0.1
        hole_coordinates[53, 2] = -0.15

        # Big hole in the center
        if hole_in_the_center:
            hole_coordinates[54, 0] = 0
            hole_coordinates[54, 1] = first_hole_coordinate + period_y * 6.5
            hole_coordinates[54, 2] = 3.5

    if hole_lattice_type == 'cloak_3':
        first_hole_coordinate = 200e-9
        period_x = 300e-9
        period_y = (sqrt(3) / 2) * period_x
        total_number_of_holes = 42
        hole_coordinates = zeros((total_number_of_holes, 3))
        hole_shapes = ['circle' for x in range(hole_coordinates.shape[0])]
        N = 0
        # Main triangle
        for i in range(2):
            hole_coordinates[N, 0] = -period_x / 2 + i * period_x
            hole_coordinates[N, 1] = first_hole_coordinate + period_y
            N += 1
        for i in range(3):
            hole_coordinates[N, 0] = -period_x + period_x * i
            hole_coordinates[N, 1] = first_hole_coordinate + period_y * 2
            N += 1
        for i in range(4):
            hole_coordinates[N, 0] = -1.5 * period_x + period_x * i
            hole_coordinates[N, 1] = first_hole_coordinate + period_y * 3
            N += 1
        for i in range(5):
            hole_coordinates[N, 0] = -2.0 * period_x + period_x * i
            hole_coordinates[N, 1] = first_hole_coordinate + period_y * 4
            N += 1
        for i in range(6):
            hole_coordinates[N, 0] = -2.5 * period_x + period_x * i
            hole_coordinates[N, 1] = first_hole_coordinate + period_y * 5
            N += 1

        # Holes on the sides
        for j in range(2):
            for i in range(5):
                hole_coordinates[N, 0] = sign(1 - 2 * j) * (
                    (i + 1) * period_x / 2 + period_x * cos(1.0 * pi / 6.0))
                hole_coordinates[N, 1] = first_hole_coordinate + period_y * (
                    i + 1) - period_y * sin(1.0 * pi / 6.0)
                N += 1
            for i in range(4):
                hole_coordinates[N, 0] = sign(1 - 2 * j) * (
                    (i + 3) * period_x / 2 + period_x * cos(1.0 * pi / 6.0))
                hole_coordinates[N, 1] = first_hole_coordinate + period_y * (
                    i + 1) - period_y * sin(1.0 * pi / 6.0)
                N += 1
            for i in range(2):
                hole_coordinates[N, 0] = sign(1 - 2 * j) * (
                    (i + 5) * period_x / 2 + period_x * cos(1.0 * pi / 6.0))
                hole_coordinates[N, 1] = first_hole_coordinate + period_y * (
                    i + 1) - period_y * sin(1.0 * pi / 6.0)
                N += 1

    if hole_lattice_type == 'turn':
        first_hole_coordinate = 300e-9
        turning_point = 10  # This is an integer number N in Lx=N*period_x where Lx is the distance form the center to the center of the turn
        number_of_periods_y_in_turn = 4  # Number of periods in the turning part
        number_of_periods_y_before_turn = 4  # Number of periods before the turn
        number_of_periods_x = 5  # Number of periods along x (should be odd number)
        hole_coordinates = zeros(
            (number_of_periods_x * number_of_periods_y_before_turn +
             number_of_periods_y_in_turn * number_of_periods_x + 4, 3))
        hole_shapes = ['circle' for x in range(hole_coordinates.shape[0])]
        hole_number = 0
        for Ny in range(number_of_periods_y_before_turn
                        ):  # This is the square lattice before the turn
            for Nx in range(number_of_periods_x):
                hole_coordinates[hole_number, 0] = -(number_of_periods_x - 1) * period_x / 2 + Nx * period_x
                hole_coordinates[hole_number, 1] = first_hole_coordinate + Ny * period_y
                hole_number += 1

        offset_x = -period_x * turning_point
        offset_y = first_hole_coordinate + period_y * (number_of_periods_y_before_turn - 1)
        turning_angle = 2.0 * asin(1.0 / (2.0 * turning_point))
        for Ny in range(number_of_periods_y_in_turn):  # This is the turning part
            for Nx in range(number_of_periods_x):
                hole_coordinates[hole_number, 0] = offset_x + (Nx + turning_point - (number_of_periods_x - 1) / 2) * period_x * cos(turning_angle * (Ny + 1))
                hole_coordinates[hole_number, 1] = offset_y + (Nx + turning_point - (number_of_periods_x - 1) / 2) * period_x * sin(turning_angle * (Ny + 1))
                #hole_coordinates[hole_number,2]=period_y*(-turning_point + (Nx+turning_point-(number_of_periods_x-1)/2))/(turning_point*period_x)
                hole_number += 1
        # Additional structures
        hole_coordinates[hole_number, 0] = -2.0 * period_x - rectangular_hole_side_x / 2
        hole_coordinates[hole_number, 1] = rectangular_hole_side_y / 2
        hole_shapes[hole_number] = 'rectangle'
        hole_coordinates[hole_number + 1, 0] = +2.0 * period_x + rectangular_hole_side_x / 2
        hole_coordinates[hole_number + 1, 1] = rectangular_hole_side_y / 2
        hole_shapes[hole_number + 1] = 'rectangle'
        hole_coordinates[hole_number + 2, 0] = -2.0 * period_x - rectangular_hole_side_x / 2
        hole_coordinates[hole_number + 2, 1] = 3 * rectangular_hole_side_y / 2
        hole_shapes[hole_number + 2] = 'triangle_up'
        hole_coordinates[hole_number + 3, 0] = +2.0 * period_x
        hole_coordinates[hole_number + 3, 1] = 3 * rectangular_hole_side_y / 2 - 100e-9
        hole_shapes[hole_number + 3] = 'triangle_down'

    if hole_lattice_type == 'turn90':
        '''Staggered lattice for turning 90 degree'''
        first_hole_coordinate = 1000e-9
        number_of_periods_x = 8
        number_of_periods_y = 3
        hole_coordinates = zeros((number_of_periods_x * number_of_periods_y * 2 + 1, 3))
        hole_shapes = ['circle' for x in range(hole_coordinates.shape[0])]
        N = 0
        for i in range(number_of_periods_y):
            for j in range(number_of_periods_x):
                hole_coordinates[ N, 0] = -(number_of_periods_x - 1) * period_x + j * period_x
                hole_coordinates[N, 1] = first_hole_coordinate + i * 2 * period_y
                hole_coordinates[ N + 1, 0] = -(number_of_periods_x) * period_x + period_x / 2 + j * period_x
                hole_coordinates[ N + 1, 1] = first_hole_coordinate + i * 2 * period_y + period_y
                N += 2

        hole_coordinates[N, 0] = (number_of_periods_x) * period_x / 2
        hole_coordinates[N, 1] = first_hole_coordinate / 2.0
        hole_shapes[N] = 'rectangle'

    if hole_lattice_type == 'cavity':
        '''Regular square lattice of holes'''
        first_hole_coordinate = 500e-9
        number_of_periods_x = 1
        number_of_periods_y = 8
        hole_coordinates = zeros((number_of_periods_x * number_of_periods_y + 2, 3))
        hole_coordinates = zeros((number_of_periods_x * number_of_periods_y + 2, 3))
        hole_shapes = ['circle' for x in range(hole_coordinates.shape[0])]
        # hole_shapes = ['rectangle' for x in range(hole_coordinates.shape[0])]
        # hole_shapes = ['triangle_up' for x in range(hole_coordinates.shape[0])]
        hole_number = 0
        for i in range(number_of_periods_y):
            for j in range(number_of_periods_x):
                hole_coordinates[hole_number, 0] = -(number_of_periods_x - 1) * period_x / 2 + j * period_x
                hole_coordinates[hole_number, 1] = first_hole_coordinate + i * period_y
                hole_number += 1

        # Top wall to create a cavity
        hole_coordinates[hole_number, 0] = 0
        hole_coordinates[hole_number, 1] = 3250e-9
        hole_shapes[hole_number] = 'rectangle'
        hole_number += 1

        # Bottom wall to create a cavity
        hole_coordinates[hole_number, 0] = 0
        hole_coordinates[hole_number, 1] = 0
        hole_shapes[hole_number] = 'rectangle'

    return hole_coordinates, hole_shapes


def pillar_positioning(pillar_lattice_type, period_x, period_y):
    """This function places pillars in space, depending on the lattice type,
    and returns coordinates of their centers and changes in their size (if any).
    In the pillar_coordinates array, column 0 is X, column 1 is Y, column 2 is correction of the base diameter.
    i.e. 0 is no correction, 1 means the pillar will be +100% larger, -0.5 means 50% smaller."""

    if pillar_lattice_type == 'square':
        number_of_periods_x = 5
        number_of_periods_y = 5
        first_pillar_coordinate = period_y / 2.0
        pillar_coordinates = zeros((number_of_periods_x * number_of_periods_y, 3))
        hole_number = 0
        for i in range(number_of_periods_y):
            for j in range(number_of_periods_x):
                pillar_coordinates[hole_number, 0] = -(number_of_periods_x - 1) * period_x / 2 + j * period_x
                pillar_coordinates[hole_number, 1] = first_pillar_coordinate + i * period_y
                hole_number += 1

    if pillar_lattice_type == 'black_silicon':
        number_of_periods_x = 6
        number_of_periods_y = 6
        disorder = 0.15
        first_pillar_coordinate = period_y / 2.0
        pillar_coordinates = zeros(
            (number_of_periods_x * number_of_periods_y, 3))
        hole_number = 0
        for i in range(number_of_periods_y):
            for j in range(number_of_periods_x):
                pillar_coordinates[hole_number, 0] = -(number_of_periods_x - 1) * period_x / 2 + j * period_x + disorder * (1.0 - 2.0 * random()) * period_x
                pillar_coordinates[hole_number, 1] = (first_pillar_coordinate + i * period_y) + disorder * (1.0 - 2.0 * random()) * period_y
                pillar_coordinates[hole_number, 2] = disorder * (1.0 - 2.0 * random())
                hole_number += 1

    if pillar_lattice_type == 'square_test':
        number_of_periods_x = 1
        number_of_periods_y = 2
        first_pillar_coordinate = period_y / 2.0
        pillar_coordinates = zeros(
            (number_of_periods_x * number_of_periods_y, 3))
        hole_number = 0
        for i in range(number_of_periods_y):
            for j in range(number_of_periods_x):
                pillar_coordinates[hole_number, 0] = -(number_of_periods_x - 1) * period_x / 2 + j * period_x
                pillar_coordinates[hole_number, 1] = first_pillar_coordinate + i * period_y
                hole_number += 1
    return pillar_coordinates


def print_lattice():
    """This function simply runs a test to print the lattice"""
    hole_lattice_type = 'yana_pnc'
    rectangular_hole_side_y = 250e-9
    rectangular_hole_side_x = 250e-9
    width = 2200e-9
    period_x = 300e-9
    period_y = 300e-9
    hole_coordinates, hole_shapes = hole_positioning(hole_lattice_type, rectangular_hole_side_y, rectangular_hole_side_x, width, period_x, period_y)

    plt.plot(hole_coordinates[:, 0], hole_coordinates[:, 1], 'o', markersize=15)
    # plt.axes().set_aspect('equal', 'datalim')
    # plt.savefig("lattice.pdf", dpi=300, format = 'pdf', bbox_inches="tight")
    plt.show()


# Generate the coordinates of holes and pillars:
pillar_coordinates = pillar_positioning(PILLAR_LATTICE_TYPE, PERIOD_X, PERIOD_Y)
hole_coordinates, hole_shapes = hole_positioning(HOLE_LATTICE_TYPE, RECTANGULAR_HOLE_SIDE_Y, RECTANGULAR_HOLE_SIDE_X, WIDTH, PERIOD_X, PERIOD_Y)


if __name__ == "__main__":
    print_lattice()