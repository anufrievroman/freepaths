"""Module that draws shape of the structure"""

from matplotlib.patches import Rectangle, Circle, Polygon


def draw_structure(cf):
    """Draw shape of the structure using patches from matplotlib"""

    # Overal shape as gray area:
    patches = [Rectangle((-1e6*cf.width/2, 0),
                         1e6*cf.width,
                         1e6*cf.length,
                         facecolor=cf.output_structure_color)]
    if cf.include_holes:
        # Holes as white patches:
        for index, shape in enumerate(cf.hole_shapes):
            x_coor = 1e6*cf.hole_coordinates[index, 0]
            y_coor = 1e6*cf.hole_coordinates[index, 1]
            scale_factor = cf.hole_coordinates[index, 2]

            if shape == 'circle':
                hole_radius = 1e6*(1 + scale_factor)*cf.circular_hole_diameter/2
                patch = Circle((x_coor, y_coor), hole_radius, facecolor='white')
                patches.append(patch)

            if shape == 'rectangle':
                hole_size_x = 1e6*(1 + scale_factor)*cf.rectangular_hole_side_x
                hole_size_y = 1e6*(1 + scale_factor)*cf.rectangular_hole_side_y
                patch = Rectangle((x_coor - hole_size_x/2, y_coor - hole_size_y/2),
                                  hole_size_x, hole_size_y, facecolor='white')
                patches.append(patch)

            if shape == 'triangle_up':
                hole_size_x = 1e6*(1 + scale_factor)*cf.rectangular_hole_side_x
                hole_size_y = 1e6*(1 + scale_factor)*cf.rectangular_hole_side_y
                patch = Polygon([[x_coor-hole_size_x/2, y_coor-hole_size_y/2],
                                 [x_coor+hole_size_x/2, y_coor-hole_size_y/2],
                                 [x_coor, y_coor+hole_size_y/2]], closed=True, facecolor='white')
                patches.append(patch)

            if shape == 'triangle_down':
                hole_size_x = 1e6*(1 + scale_factor)*cf.rectangular_hole_side_x
                hole_size_y = 1e6*(1 + scale_factor)*cf.rectangular_hole_side_y
                patch = Polygon([[x_coor-hole_size_x/2, y_coor+hole_size_y/2],
                                 [x_coor+hole_size_x/2, y_coor+hole_size_y/2],
                                 [x_coor, y_coor-hole_size_y/2]], closed=True, facecolor='white')
                patches.append(patch)

    return patches
