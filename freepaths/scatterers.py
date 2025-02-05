"""
Module that contains the Hole classes.
These classes contain all methods associated with scattering on holes.
"""


from math import atan
from numpy import pi, array, linspace, column_stack, vstack
from random import random
from matplotlib.patches import Rectangle, Circle, Polygon
from scipy.spatial import cKDTree

from freepaths.scattering_primitives import *
from freepaths.scattering_types import ScatteringTypes


class Hole:
    def is_inside(self, x, y, z, cf) -> bool:
        """
        Check if phonon with given coordinates traverses the boundary.
        It returns True or False depending whether x, y, z is inside the hole
        """
        pass

    def scatter(self, ph, scattering_types, x, y, z, cf):
        """
        Calculate the new direction after scattering on the hole.
        It returns ScatteringTypes object with the scattering type that occued
        """
        pass

    def get_patch(self, color_holes, cf):
        """
        Create a patch in the shape of the hole to use in the plots.
        It returns the matplotlib.patches objects like Rectangle etc.
        """
        pass


class CircularHole(Hole):
    """Shape of a circular hole"""

    def __init__(self, x=0, y=0, diameter=100e-9):
        self.x0 = x
        self.y0 = y
        self.diameter = diameter

    def is_inside(self, x, y, z, cf):
        """Check if phonon with given coordinates traverses the boundary"""
        radius = self.diameter / 2
        return (x - self.x0) ** 2 + (y - self.y0) ** 2 <= radius**2

    def scatter(self, ph, scattering_types, x, y, z, cf):
        """Calculate the new direction after scattering on the hole"""
        if y == self.y0:
            y += 1e-9  # Prevent division by zero
        tangent_theta = atan((x - self.x0) / (y - self.y0))

        # check if the phonon is travelling towards the hole
        current_distance = sqrt((self.x0 - ph.x)**2 + (self.y0 - ph.y)**2)
        next_distance = sqrt((self.x0 - x)**2 + (self.y0 - y)**2)
        if next_distance <= current_distance:
            scattering_types.holes = circle_outer_scattering(
                ph, tangent_theta, y, self.y0, cf.hole_roughness, cf
            )

    def get_patch(self, color_holes, cf):
        """Create a patch in the shape of the hole to use in the plots"""
        return Circle(
            (1e6 * self.x0, 1e6 * self.y0),
            1e6 * self.diameter / 2,
            facecolor=color_holes,
        )


class RectangularHole(Hole):
    """Shape of a rectangular hole"""

    def __init__(self, x=0, y=0, size_x=100e-9, size_y=100e-9, depth=None):
        self.x0 = x
        self.y0 = y
        self.size_x = size_x
        self.size_y = size_y
        self.depth = depth

    def is_inside(self, x, y, z, cf):
        """Check if phonon with given coordinates traverses the boundary. It also depens on in the hole is complete or partial."""
        if self.depth and z:
            return (abs(x - self.x0) <= self.size_x / 2) and (abs(y - self.y0) <= self.size_y / 2) and (z > cf.thickness/2 - self.depth)
        else:
            return (abs(x - self.x0) <= self.size_x / 2) and (abs(y - self.y0) <= self.size_y / 2)


    def scatter(self, ph, scattering_types, x, y, z, cf):
        """Calculate the new direction after scattering on the hole"""

        # If phonon arrives from below, then it's bottom scattering:
        if self.depth and ph.z < (cf.thickness/2 - self.depth):
            scattering_types.holes = in_plane_surface_scattering(ph, cf.top_roughness)
            return

        # Otherwise, calculate coordinate y of the intersection with the hole side:
        y1 = (self.y0 - y) + cos(ph.theta) * (self.size_x / 2 - abs(self.x0 - x)) / abs(sin(ph.theta))

        # Scattering on the left wall:
        if abs(y1) <= self.size_y / 2 and x < self.x0:
            scattering_types.holes = vertical_surface_left_scattering(ph, cf.hole_roughness, cf)

        # Scattering on the right wall:
        elif abs(y1) <= self.size_y / 2 and x > self.x0:
            scattering_types.holes = vertical_surface_right_scattering(ph, cf.hole_roughness, cf)

        # Scattering on the top wall:
        elif y > self.y0:
            scattering_types.holes = horizontal_surface_up_scattering(ph, cf.hole_roughness)

        # Scattering on the bottom wall:
        else:
            scattering_types.holes = horizontal_surface_down_scattering(ph, cf.hole_roughness)

    def get_patch(self, color_holes, cf):
        """Create a patch in the shape of the hole to use in the plots"""
        return Rectangle(
            (1e6 * (self.x0 - self.size_x / 2), 1e6 * (self.y0 - self.size_y / 2)),
            1e6 * self.size_x,
            1e6 * self.size_y,
            facecolor=color_holes,
        )


class TriangularUpHole(Hole):
    """Shape of a triangular hole facing up"""

    def __init__(self, x=0, y=0, size_x=100e-9, size_y=100e-9):
        self.x0 = x
        self.y0 = y
        self.size_x = size_x
        self.size_y = size_y

    def is_inside(self, x, y, z, cf):
        """Check if phonon with given coordinates traverses the boundary"""
        beta = atan(0.5 * self.size_x / self.size_y)
        return ((self.size_y / 2 + (y - self.y0) <= (self.size_x / 2 - abs(x - self.x0)) / tan(beta))
                and (abs(y - self.y0) < self.size_y / 2))

    def scatter(self, ph, scattering_types, x, y, z, cf):
        """Calculate the new direction after scattering on the hole"""

        # Scattering on the bottom wall of the triangle:
        if (ph.y < self.y0 - self.size_y / 2) and (abs(ph.theta) < pi / 2):
            scattering_types.holes = horizontal_surface_down_scattering(ph, cf.hole_roughness)
            # triangle_scattering_places.floor = scattering_types.holes

        # Scattering on the sidewalls of the triangle:
        else:
            beta = atan(0.5 * self.size_x / self.size_y)
            scattering_types.holes = inclined_surfaces_up_scattering(ph, beta, x, self.x0, cf.hole_roughness)
            # if x > self.x0:
                # triangle_scattering_places.right_wall = scattering_types.holes
            # else:
                # triangle_scattering_places.left_wall = scattering_types.holes

    def get_patch(self, color_holes, cf):
        """Create a patch in the shape of the hole to use in the plots"""
        return Polygon(
            [
                [1e6 * (self.x0 - self.size_x / 2), 1e6 * (self.y0 - self.size_y / 2)],
                [1e6 * (self.x0 + self.size_x / 2), 1e6 * (self.y0 - self.size_y / 2)],
                [1e6 * self.x0, 1e6 * (self.y0 + self.size_y / 2)],
            ],
            closed=True,
            facecolor=color_holes,
        )


class TriangularDownHole(Hole):
    """Shape of a triangular hole facing down"""

    def __init__(self, x=0, y=0, size_x=100e-9, size_y=100e-9):
        self.x0 = x
        self.y0 = y
        self.size_x = size_x
        self.size_y = size_y

    def is_inside(self, x, y, z, cf):
        """Check if phonon with given coordinates traverses the boundary"""
        beta = atan(0.5 * self.size_x / self.size_y)
        return ((self.size_y / 2 - (y - self.y0) <= (self.size_x / 2 - abs(x - self.x0)) / tan(beta))
                and (abs(y - self.y0) < self.size_y / 2))

    def scatter(self, ph, scattering_types, x, y, z, cf):
        """Calculate the new direction after scattering on the hole"""

        # Angle of the triangle:
        beta = atan(0.5 * self.size_x / self.size_y)

        # Scattering on the top wall of the triangle:
        if (ph.y > self.y0 + self.size_y / 2) and (abs(ph.theta) > pi / 2):
            scattering_types.holes = horizontal_surface_up_scattering(ph, cf.hole_roughness)

        # Scattering on the sidewalls of the triangle:
        else:
            scattering_types.holes = inclined_surfaces_down_scattering(ph, beta, x, self.x0, cf.hole_roughness)

    def get_patch(self, color_holes, cf):
        """Create a patch in the shape of the hole to use in the plots"""
        return Polygon(
            [
                [1e6 * (self.x0 - self.size_x / 2), 1e6 * (self.y0 + self.size_y / 2)],
                [1e6 * (self.x0 + self.size_x / 2), 1e6 * (self.y0 + self.size_y / 2)],
                [1e6 * self.x0, 1e6 * (self.y0 - self.size_y / 2)],
            ],
            closed=True,
            facecolor=color_holes,
        )


class TriangularDownHalfHole(Hole):
    """Shape of a half triangular hole facing down"""

    def __init__(self, x=0, y=0, size_x=100e-9, size_y=100e-9, is_right_half=True):
        self.x0 = x
        self.y0 = y
        self.size_x = size_x
        self.size_y = size_y
        self.is_right_half = is_right_half

    def is_inside(self, x, y, z, cf):
        """Check if phonon with given coordinates traverses the boundary"""
        beta = atan(0.5 * self.size_x / self.size_y)
        if self.is_right_half:
            return ((self.size_y / 2 - (y - self.y0) <= (self.size_x / 2 - abs(x - self.x0)) / tan(beta))
            and (abs(y - self.y0) < self.size_y / 2) and (x > self.x0))
        else:
           return ((self.size_y / 2 - (y - self.y0) <= (self.size_x / 2 - abs(x - self.x0)) / tan(beta))
            and (abs(y - self.y0) < self.size_y / 2) and (x < self.x0))

    def scatter(self, ph, scattering_types, x, y, z, cf):
        """Calculate the new direction after scattering on the hole"""

        # Angle of the triangle:
        beta = atan(0.5 * self.size_x / self.size_y)

        in_area = self.is_inside(x, y, z, cf)

        # Triangle oriented right:
        if self.is_right_half:
            # Scattering on the top wall of the triangle:
            if (ph.y > self.y0 + self.size_y / 2) and (abs(ph.theta) > pi / 2):
                scattering_types.holes = horizontal_surface_up_scattering(ph, cf.hole_roughness)

            # Scattering on the vertical sidewall of the triangle:
            elif ph.x < self.x0:
                scattering_types.holes = vertical_surface_left_scattering(ph, cf.hole_roughness, cf)

            # Scattering on the inclined sidewall of the triangle:
            else:
                scattering_types.holes = inclined_surfaces_down_scattering(ph, beta, x, self.x0, cf.hole_roughness)

        # Triangle oriented left:
        else:
            # Scattering on the top wall of the triangle:
            if (ph.y > self.y0 + self.size_y / 2) and (abs(ph.theta) > pi / 2):
                scattering_types.holes = horizontal_surface_up_scattering(
                    ph, cf.hole_roughness
                )

            # Scattering on the vertical sidewall of the triangle:
            elif ph.x > self.x0:
                scattering_types.holes = vertical_surface_right_scattering(ph, cf.hole_roughness, cf)

            # Scattering on the inclined sidewall of the triangle:
            else:
                scattering_types.holes = inclined_surfaces_down_scattering(ph, beta, x, self.x0, cf.hole_roughness)

    def get_patch(self, color_holes, cf):
        """Create a patch in the shape of the hole to use in the plots"""
        if self.is_right_half:
            return Polygon(
                [
                    [1e6 * (self.x0), 1e6 * (self.y0 + self.size_y / 2)],
                    [
                        1e6 * (self.x0 + self.size_x / 2),
                        1e6 * (self.y0 + self.size_y / 2),
                    ],
                    [1e6 * self.x0, 1e6 * (self.y0 - self.size_y / 2)],
                ],
                closed=True,
                facecolor=color_holes,
            )
        else:
            return Polygon(
                [
                    [
                        1e6 * (self.x0 - self.size_x / 2),
                        1e6 * (self.y0 + self.size_y / 2),
                    ],
                    [1e6 * (self.x0), 1e6 * (self.y0 + self.size_y / 2)],
                    [1e6 * self.x0, 1e6 * (self.y0 - self.size_y / 2)],
                ],
                closed=True,
                facecolor=color_holes,
            )


class TriangularUpHalfHole(Hole):
    """Shape of a half triangular hole facing up"""

    def __init__(self, x=0, y=0, size_x=100e-9, size_y=100e-9, is_right_half=True):
        self.x0 = x
        self.y0 = y
        self.size_x = size_x
        self.size_y = size_y
        self.is_right_half = is_right_half

    def is_inside(self, x, y, z, cf):

        """Check if phonon with given coordinates traverses the boundary"""
        beta = atan(0.5 * self.size_x / self.size_y)

        if self.is_right_half:
            return (( self.size_y / 2 + (y - self.y0) <= (self.size_x / 2 - abs(x - self.x0)) / tan(beta))
                    and (abs(y - self.y0) < self.size_y / 2) and (x > self.x0))
        else:
            return ((self.size_y / 2 + (y - self.y0) <= (self.size_x / 2 - abs(x - self.x0)) / tan(beta))
                    and (abs(y - self.y0) < self.size_y / 2) and (x < self.x0))

    def scatter(self, ph, scattering_types, x, y, z, cf):
        """Calculate the new direction after scattering on the hole"""

        # Angle of the triangle:
        beta = atan(0.5 * self.size_x / self.size_y)

        in_area = self.is_inside(x, y, z, cf)

        # Triangle oriented right:
        if self.is_right_half:
            # Scattering on the bottom wall of the triangle:
            if (ph.y < self.y0 - self.size_y / 2) and (abs(ph.theta) < pi / 2):
                scattering_types.holes = horizontal_surface_down_scattering(ph, cf.hole_roughness)

            # Scattering on the vertical sidewall of the triangle:
            elif ph.x < self.x0:
                scattering_types.holes = vertical_surface_left_scattering(ph, cf.hole_roughness, cf)

            # Scattering on the inclined sidewall of the triangle:
            else:
                scattering_types.holes = inclined_surfaces_up_scattering(ph, beta, x, self.x0, cf.hole_roughness)

        # Triangle oriented left:
        else:
            # Scattering on the bottom wall of the triangle:
            if (ph.y < self.y0 - self.size_y / 2) and (abs(ph.theta) < pi / 2):
                scattering_types.holes = horizontal_surface_down_scattering(ph, cf.hole_roughness)

            # Scattering on the vertical sidewall of the triangle:
            elif ph.x > self.x0:
                scattering_types.holes = vertical_surface_right_scattering(ph, cf.hole_roughness, cf)

            # Scattering on the inclined sidewall of the triangle:
            else:
                scattering_types.holes = inclined_surfaces_up_scattering(ph, beta, x, self.x0, cf.hole_roughness)

    def get_patch(self, color_holes, cf):
        """Create a patch in the shape of the hole to use in the plots"""
        if self.is_right_half:
            return Polygon(
                [
                    [1e6 * (self.x0), 1e6 * (self.y0 - self.size_y / 2)],
                    [
                        1e6 * (self.x0 + self.size_x / 2),
                        1e6 * (self.y0 - self.size_y / 2),
                    ],
                    [1e6 * self.x0, 1e6 * (self.y0 + self.size_y / 2)],
                ],
                closed=True,
                facecolor=color_holes,
            )
        else:
            return Polygon(
                [
                    [
                        1e6 * (self.x0 - self.size_x / 2),
                        1e6 * (self.y0 - self.size_y / 2),
                    ],
                    [1e6 * (self.x0), 1e6 * (self.y0 - self.size_y / 2)],
                    [1e6 * self.x0, 1e6 * (self.y0 + self.size_y / 2)],
                ],
                closed=True,
                facecolor=color_holes,
            )


class PointLineHole(Hole):
    """General shape that can be defined by a list of points"""

    def __init__(self, x=0, y=0, points=None, thickness=100e-9, rotation=0):
        # add option for rounded/angled corners?

        assert points is not None, "Please provide some points to the PointLineHole"

        # Rotate points:
        if rotation != 0 and rotation is not None:
            points = self.rotate_points(points, rotation)

        # Move points to x0, y0 position:
        self.points = array(points) + (x, y)

        # Build cKDTree for fast search:
        self.tree = cKDTree(self.points)
        self.thickness = thickness

    def is_inside(self, x, y, z, cf):
        """Check if phonon with given coordinates traverses the boundary"""
        distance, _ = self.tree.query((x, y))
        return distance < self.thickness / 2

    def scatter(self, ph, scattering_types, x, y, z, cf):
        """Calculate the new direction after scattering on the hole"""

        # First, find nearest circle and get its coordinates:
        distance, index = self.tree.query((x, y))
        x0, y0 = self.points[index]
        if y == y0:
            y += 1e-9  # Prevent division by zero
        tangent_theta = atan((x - x0) / (y - y0))

        # Check if the phonon is traveling towards the hole:
        current_distance, _ = self.tree.query((ph.x, ph.y))
        if distance <= current_distance:
            scattering_types.holes = circle_outer_scattering(ph, tangent_theta, y, y0, cf.hole_roughness, cf)

    def get_patch(self, color_holes, cf):
        """Create a patch in the shape of the hole to use in the plots"""
        return [Circle((x*1e6, y*1e6), self.thickness*1e6/2, facecolor=color_holes,) for x, y in self.points]


    def rotate_points(self, points, angle):
        rotated_points = []
        cos_theta = cos(-angle/180*pi)
        sin_theta = sin(-angle/180*pi)

        for point in points:
            x = point[0]
            y = point[1]

            # Perform rotation using rotation matrix
            new_x = x * cos_theta - y * sin_theta
            new_y = x * sin_theta + y * cos_theta

            rotated_points.append((new_x, new_y))

        return rotated_points


class FunctionLineHole(PointLineHole):
    """Create a line of holes from a mathematical function"""

    def __init__(self, x=0, y=0, thickness=60e-9, function=lambda x: sin(x*2*pi/300e-9)/2*200e-9, function_range=(-150e-9, 150e-9),
                size_x=None, size_y=None, resolution=1e-9, rotation=0):
        points = self.points_from_function(function, function_range, size_x, size_y, resolution, thickness)
        super().__init__(x, y, points, thickness, rotation)

    def points_from_function(self, function, function_range, size_x, size_y, resolution, thickness):
        """Generate the points in the function space"""
        number_of_circles = round((function_range[1] - function_range[0])/resolution if size_x is None else (size_x-thickness)/resolution)

        xs = linspace(function_range[0], function_range[1], number_of_circles)
        ys = array([0.0]*len(xs))
        for i, x in enumerate(xs):
            ys[i] = function(x)

        # normalize the point to a range of 1 and then rescale to wanted dimensions and account for the thickness of the shape
        if size_x is not None:
            xs /= function_range[1] - function_range[0]
            xs *= size_x - thickness

        if size_y is not None:
            ys /= max(ys) - min(ys)
            ys *= size_y - thickness

        return vstack((xs, ys)).T


class ParabolaTop(Hole):
    """Shape of a parabolic wall"""

    def __init__(self, tip=0, focus=0):
        self.tip = tip
        self.focus = focus

    def is_inside(self, x, y, z, cf):
        """Check if phonon with given coordinates traverses the boundary"""
        y_cept = -((cf.width / 2) ** 2) / (4 * self.focus) + self.tip
        return (y > y_cept) and (x**2 + 4 * self.focus * (y - self.tip)) >= 0

    def scatter(self, ph, scattering_types, x, y, z, cf):
        """Calculate the new direction after scattering on the hole"""

        # Calculate angle to the surface and specular scattering probability:
        normal_theta = pi * (x < 0) - atan(2 * self.focus / x)
        dot_product = cos(ph.phi) * sin(ph.theta - normal_theta)
        angle = acos(dot_product)
        p = specularity(angle, cf.side_wall_roughness, ph.wavelength)

        # Specular scattering:
        if random() < p:
            if abs(ph.theta) > pi / 2:
                ph.theta = ph.theta - 2 * normal_theta
            else:
                ph.theta = 2 * normal_theta - ph.theta
            scattering_types.walls = Scattering.SPECULAR

        # Diffuse scattering:
        else:
            scattering_types.walls = Scattering.DIFFUSE
            for _ in range(10):
                # Lambert distribution
                ph.theta = normal_theta + asin(2 * random() - 1) - pi / 2
                ph.phi = asin((asin(2 * random() - 1)) / (pi / 2))

                # Accept the angles only if they do not immediately cause new scattering:
                if no_new_scattering(ph, cf):
                    break

    def get_patch(self, color_holes, cf):
        """Create a patch in the shape of the hole to use in the plots"""
        eval_xs = linspace(-cf.width / 2, cf.width / 2, 50)
        parabola_ys = -(eval_xs**2) / (4 * self.focus) + self.tip
        parabola_points = column_stack((eval_xs, parabola_ys))
        polygon_point = vstack(
            (parabola_points, [cf.width / 2, cf.length], [-cf.width / 2, cf.length])
        )
        return Polygon(polygon_point * 1e6, closed=True, facecolor=color_holes,)


class ParabolaBottom(Hole):
    """Shape of a parabolic wall"""

    def __init__(self, tip=0, focus=0):
        self.tip = tip
        self.focus = focus

    def is_inside(self, x, y, z, cf):
        """Check if phonon with given coordinates traverses the boundary"""
        y_cept = (cf.width / 2) ** 2 / (4 * self.focus) + self.tip
        return (y < y_cept) and (x**2 - 4 * self.focus * (y - self.tip)) >= 0

    def scatter(self, ph, scattering_types, x, y, z, cf):
        """Calculate the new direction after scattering on the hole"""

        # Calculate angle to the surface and specular scattering probability:
        normal_theta = pi * (x < 0) - atan(-2 * self.focus / x)
        dot_product = cos(ph.phi) * sin(ph.theta - normal_theta)
        angle = acos(dot_product)
        p = specularity(angle, cf.side_wall_roughness, ph.wavelength)

        # Specular scattering:
        if random() < p:
            ph.theta = - ph.theta + 2 * normal_theta
            scattering_types.walls = Scattering.SPECULAR

        # Diffuse scattering:
        else:
            scattering_types.walls = Scattering.DIFFUSE
            for _ in range(10):
                # Lambertian distribution
                ph.theta = normal_theta + asin(2 * random() - 1) - pi / 2
                ph.phi = asin((asin(2 * random() - 1)) / (pi / 2))

                # Accept the angles only if they do not immediately cause new scattering:
                if no_new_scattering(ph, cf):
                    break

    def get_patch(self, color_holes, cf):
        """Create a patch in the shape of the hole to use in the plots"""
        eval_xs = linspace(-cf.width / 2, cf.width / 2, 50)
        parabola_ys = eval_xs**2 / (4 * self.focus) + self.tip
        parabola_points = column_stack((eval_xs, parabola_ys))
        polygon_point = vstack((parabola_points, [cf.width / 2, 0], [-cf.width / 2, 0]))
        return Polygon(polygon_point * 1e6, closed=True, facecolor=color_holes,)


class CircularPillar():
    """Shape of a circular pillar with inclined wall"""


    def __init__(self, x=0, y=0, diameter=200e-9, height=300e-9, wall_angle=pi / 2):
        self.x0 = x
        self.y0 = y
        self.diameter = diameter
        self.height = height
        self.wall_angle = wall_angle

    # did not add the is inside function because it requires the phonon speed

    def check_if_scattering(self, ph, scattering_types, x, y, z, cf):
        """Check if a phonon strikes a circular pillar and calculate new direction"""

        # Cone radius at a given z coordinate:
        radius = self.diameter / 2  # - (z - cf.thickness / 2) / tan(pillar.wall_angle)
        distance_from_pillar_center = sqrt((x - self.x0) ** 2 + (y - self.y0) ** 2)
        # distance_from_pillar_center_original = sqrt(
        #     (ph.x - self.x) ** 2 + (ph.y - self.y) ** 2
        # )
        step = 2 * ph.speed * cf.timestep

        # If phonon crosses the pillar boundary. Third condition is to exclude all other pillars:
        if (
            distance_from_pillar_center >= radius
            and z > cf.thickness / 2
            and distance_from_pillar_center < radius + step
        ):
            # Calculate angle to the surface and specular scattering probability:
            tangent_theta = atan((x - self.x0) / (y - self.y0))
            scattering_types.pillars = circle_inner_scattering(
                ph, tangent_theta, y, self.y0, cf.pillar_roughness
            )

    def get_patch(self, color_holes, cf):
        """Create a patch in the shape of the hole to use in the plots"""
        return Circle((1e6 * self.x0, 1e6 * self.y0), 1e6 * self.diameter / 2, facecolor=color_holes,)


class Interface:
    def is_crossed(self, ph, x, y, z) -> bool:
        """
        Check if phonon with given coordinates traverses the plane.
        It returns True or False depending whether x, y, z are on the other side.
        """
        pass

    def scatter(self, ph, scattering_types, x, y, z, cf):
        """
        Calculate the new direction after scattering on the interface wall.
        It returns ScatteringTypes object with the scattering type that occured.
        """
        pass

    def is_transmitted(self) -> bool:
        """
        Check if phonon transmitted without scattering, taking into account the probability.
        It return True if transmission occured or False if it must scatter.
        """
        pass

    def get_patch(self, color_holes, cf):
        """
        Create a patch in the shape of a thin line to use in the plots.
        It returns the matplotlib.patches objects like Rectangle etc.
        """
        pass


class VerticalPlane(Interface):
    """Vertical plane that represents an interface"""

    def __init__(self, position_x=0, transmission=0):
        self.position_x = position_x
        self.transmission = transmission

    def is_crossed(self, ph, x, y, z):
        """Check if phonon with traverses the vertical plane at given coordinate"""
        return (ph.x < self.position_x < x) or (ph.x > self.position_x > x)

    def is_transmitted(self):
        """Check if phonon traverses the plane given the transmission probability"""
        return random() < self.transmission

    def scatter(self, ph, scattering_types, x, y, z, cf):
        """Calculate the new direction after scattering on the interface wall"""

        # Scattering on the left wall:
        if x < self.position_x:
            scattering_types.interfaces = vertical_surface_left_scattering(ph, cf.interface_roughness, cf)

        # Scattering on the right wall:
        else:
            scattering_types.interfaces = vertical_surface_right_scattering(ph, cf.interface_roughness, cf)

    def get_patch(self, color_holes, cf):
        """Create a patch in the shape of a thin line to use in the plots"""
        return Rectangle(
            (1e6 * self.position_x, 0),
            1e6 * 0.005*cf.width, 1e6 * cf.length,
            facecolor=color_holes,
        )



class HorizontalPlane(Interface):
    """Horizontal plane that represents an interface"""

    def __init__(self, position_z=0, transmission=0):
        self.position_z = position_z
        self.transmission = transmission

    def is_crossed(self, ph, x, y, z):
        """Check if phonon with traverses the vertical plane at given coordinate"""
        return (ph.z < self.position_z < z) or (ph.z > self.position_z > z)

    def is_transmitted(self):
        """Check if phonon traverses the plane given the transmission probability"""
        return random() < self.transmission

    def scatter(self, ph, scattering_types, x, y, z, cf):
        """Calculate the new direction after scattering on the interface wall"""

        # Scattering on the left wall:
        if x < self.position_z:
            scattering_types.interfaces = horizontal_surface_up_scattering(ph, cf.interface_roughness, cf)

        # Scattering on the right wall:
        else:
            scattering_types.interfaces = horizontal_surface_down_scattering(ph, cf.interface_roughness, cf)

    def get_patch(self, color_holes, cf):
        """Create a patch in the shape of a thin line to use in the plots"""
        return Rectangle(
            (0, 1e6 * self.position_z),
            1e6 * cf.length, 1e6 * 0.005*cf.thickness,
            facecolor=color_holes,
        )
