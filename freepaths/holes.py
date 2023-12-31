"""
Module that contains the Hole classes. 
These classes contain all methods associated with the hole.

The functions from the scattering_parabolic, scattering_primitives, etc... are not in here yet.
"""


from math import atan
from numpy import pi, array, linspace, column_stack, vstack
from matplotlib.patches import Rectangle, Circle, Polygon


from freepaths.scattering_primitives import *


class CircularHole:
    """Shape of a circular hole"""

    def __init__(self, x=0, y=0, diameter=100e-9):
        self.x0 = x
        self.y0 = y
        self.diameter = diameter

    def check_if_scattering(self, ph, scattering_types, x, y, z, cf):
        """Check if a phonon strikes a circular hole and calculate the new direction"""

        # If phonon is inside the circle with given radius:
        radius = self.diameter / 2
        if (x - self.x0) ** 2 + (y - self.y0) ** 2 <= radius**2:
            if y == self.y0:
                y += 1e-9  # Prevent division by zero
            tangent_theta = atan((x - self.x0) / (y - self.y0))
            scattering_types.holes = circle_outer_scattering(
                ph, tangent_theta, y, self.y0, cf.hole_roughness, cf
            )
            return True
        return False

    def get_patch(self, color_holes, cf):
        return Circle(
            (1e6 * self.x0, 1e6 * self.y0),
            1e6 * self.diameter / 2,
            facecolor=color_holes,
        )


class RectangularHole:
    """Shape of a rectangular hole"""

    def __init__(self, x=0, y=0, size_x=100e-9, size_y=100e-9):
        self.x0 = x
        self.y0 = y
        self.size_x = size_x
        self.size_y = size_y

    def check_if_scattering(self, ph, scattering_types, x, y, z, cf):
        """Check if the phonon strikes a rectangular hole and calculate new direction"""

        # If the phonon is inside the rectangle:
        if (abs(x - self.x0) <= self.size_x / 2) and (
            abs(y - self.y0) <= self.size_y / 2
        ):
            # Coordinate y of the intersection with the hole side:
            y1 = (self.y0 - y) + cos(ph.theta) * (
                self.size_x / 2 - abs(self.x0 - x)
            ) / abs(sin(ph.theta))

            # Scattering on the left wall:
            if abs(y1) <= self.size_y / 2 and x < self.x0:
                scattering_types.holes = vertical_surface_left_scattering(
                    ph, cf.hole_roughness, cf
                )

            # Scattering on the right wall:
            elif abs(y1) <= self.size_y / 2 and x > self.x0:
                scattering_types.holes = vertical_surface_right_scattering(
                    ph, cf.hole_roughness, cf
                )

            # Scattering on the top wall:
            elif y > self.y0:
                scattering_types.holes = horizontal_surface_up_scattering(
                    ph, cf.side_wall_roughness
                )

            # Scattering on the bottom wall:
            else:
                scattering_types.holes = horizontal_surface_down_scattering(
                    ph, cf.side_wall_roughness
                )

            return True
        return False

    def get_patch(self, color_holes, cf):
        return Rectangle(
            (1e6 * (self.x0 - self.size_x / 2), 1e6 * (self.y0 - self.size_y / 2)),
            1e6 * self.size_x,
            1e6 * self.size_y,
            facecolor=color_holes,
        )


class TriangularUpHole:
    """Shape of a triangular hole facing up"""

    def __init__(self, x=0, y=0, size_x=100e-9, size_y=100e-9):
        self.x0 = x
        self.y0 = y
        self.size_x = size_x
        self.size_y = size_y

    def check_if_scattering(self, ph, scattering_types, x, y, z, cf):
        """Check if the phonon strikes a reverse triangular hole and calculate new direction"""

        # Angle of the triangle:
        beta = atan(0.5 * self.size_x / self.size_y)

        # If phonon is inside the triangle:
        if (
            self.size_y / 2 + (y - self.y0)
            <= (self.size_x / 2 - abs(x - self.x0)) / tan(beta)
        ) and (abs(y - self.y0) < self.size_y / 2):
            # Scattering on the bottom wall of the triangle:
            if (ph.y < self.y0 - self.size_y / 2) and (abs(ph.theta) < pi / 2):
                scattering_types.holes = horizontal_surface_down_scattering(
                    ph, cf.hole_roughness
                )

            # Scattering on the sidewalls of the triangle:
            else:
                scattering_types.holes = inclined_surfaces_up_scattering(
                    ph, beta, x, self.x0, cf.hole_roughness
                )

            return True
        return False

    def get_patch(self, color_holes, cf):
        return Polygon(
            [
                [1e6 * (self.x0 - self.size_x / 2), 1e6 * (self.y0 - self.size_y / 2)],
                [1e6 * (self.x0 + self.size_x / 2), 1e6 * (self.y0 - self.size_y / 2)],
                [1e6 * self.x0, 1e6 * (self.y0 + self.size_y / 2)],
            ],
            closed=True,
            facecolor=color_holes,
        )


class TriangularDownHole:
    """Shape of a triangular hole facing down"""

    def __init__(self, x=0, y=0, size_x=100e-9, size_y=100e-9):
        self.x0 = x
        self.y0 = y
        self.size_x = size_x
        self.size_y = size_y

    def check_if_scattering(self, ph, scattering_types, x, y, z, cf):
        """Check if the phonon strikes a reverse triangular hole and calculate new direction after the scattering"""

        # Angle of the triangle:
        beta = atan(0.5 * self.size_x / self.size_y)

        # If phonon is inside the triangle:
        if (
            self.size_y / 2 - (y - self.y0)
            <= (self.size_x / 2 - abs(x - self.x0)) / tan(beta)
        ) and (abs(y - self.y0) < self.size_y / 2):
            # Scattering on the top wall of the triangle:
            if (ph.y > self.y0 + self.size_y / 2) and (abs(ph.theta) > pi / 2):
                scattering_types.holes = horizontal_surface_up_scattering(
                    ph, cf.hole_roughness
                )

            # Scattering on the sidewalls of the triangle:
            else:
                scattering_types.holes = inclined_surfaces_down_scattering(
                    ph, beta, x, self.x0, cf.hole_roughness
                )

            return True
        return False

    def get_patch(self, color_holes, cf):
        return Polygon(
            [
                [1e6 * (self.x0 - self.size_x / 2), 1e6 * (self.y0 + self.size_y / 2)],
                [1e6 * (self.x0 + self.size_x / 2), 1e6 * (self.y0 + self.size_y / 2)],
                [1e6 * self.x0, 1e6 * (self.y0 - self.size_y / 2)],
            ],
            closed=True,
            facecolor=color_holes,
        )


class TriangularDownHalfHole:
    """Shape of a half triangular hole facing down"""

    def __init__(self, x=0, y=0, size_x=100e-9, size_y=100e-9, is_right_half=True):
        self.x0 = x
        self.y0 = y
        self.size_x = size_x
        self.size_y = size_y
        self.is_right_half = is_right_half

    def check_if_scattering(self, ph, scattering_types, x, y, z, cf):
        """Check if the phonon strikes a reverse triangular hole and calculate new direction"""

        # Angle of the triangle:
        beta = atan(0.5 * self.size_x / self.size_y)

        # If phonon is inside the right side of the triangle:
        if (
            self.is_right_half
            and (
                self.size_y / 2 - (y - self.y0)
                <= (self.size_x / 2 - abs(x - self.x0)) / tan(beta)
            )
            and (abs(y - self.y0) < self.size_y / 2)
            and (x > self.x0)
        ):
            # Scattering on the top wall of the triangle:
            if (ph.y > self.y0 + self.size_y / 2) and (abs(ph.theta) > pi / 2):
                scattering_types.holes = horizontal_surface_up_scattering(
                    ph, cf.hole_roughness
                )

            # Scattering on the vertical sidewall of the triangle:
            elif ph.x < self.x0:
                scattering_types.holes = vertical_surface_left_scattering(
                    ph, cf.hole_roughness, cf
                )

            # Scattering on the inclined sidewall of the triangle:
            else:
                scattering_types.holes = inclined_surfaces_down_scattering(
                    ph, beta, x, self.x0, cf.hole_roughness
                )

            return True

        # If phonon is inside the left side of the triangle:
        if (
            not self.is_right_half
            and (
                self.size_y / 2 - (y - self.y0)
                <= (self.size_x / 2 - abs(x - self.x0)) / tan(beta)
            )
            and (abs(y - self.y0) < self.size_y / 2)
            and (x < self.x0)
        ):
            # Scattering on the top wall of the triangle:
            if (ph.y > self.y0 + self.size_y / 2) and (abs(ph.theta) > pi / 2):
                scattering_types.holes = horizontal_surface_up_scattering(
                    ph, cf.hole_roughness
                )

            # Scattering on the vertical sidewall of the triangle:
            elif ph.x > self.x0:
                scattering_types.holes = vertical_surface_right_scattering(
                    ph, cf.hole_roughness, cf
                )

            # Scattering on the inclined sidewall of the triangle:
            else:
                scattering_types.holes = inclined_surfaces_down_scattering(
                    ph, beta, x, self.x0, cf.hole_roughness
                )

            return True
        return False

    def get_patch(self, color_holes, cf):
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


class TriangularUpHalfHole:
    """Shape of a half triangular hole facing up"""

    def __init__(self, x=0, y=0, size_x=100e-9, size_y=100e-9, is_right_half=True):
        self.x0 = x
        self.y0 = y
        self.size_x = size_x
        self.size_y = size_y
        self.is_right_half = is_right_half

    def check_if_scattering(self, ph, scattering_types, x, y, z, cf):
        """Check if the phonon strikes a reverse triangular hole and calculate new direction"""

        # Angle of the triangle:
        beta = atan(0.5 * self.size_x / self.size_y)

        # If phonon is inside the right side of the triangle:
        if (
            self.is_right_half
            and (
                self.size_y / 2 + (y - self.y0)
                <= (self.size_x / 2 - abs(x - self.x0)) / tan(beta)
            )
            and (abs(y - self.y0) < self.size_y / 2)
            and (x > self.x0)
        ):
            # Scattering on the bottom wall of the triangle:
            if (ph.y < self.y0 - self.size_y / 2) and (abs(ph.theta) < pi / 2):
                scattering_types.holes = horizontal_surface_down_scattering(
                    ph, cf.hole_roughness
                )

            # Scattering on the vertical sidewall of the triangle:
            elif ph.x < self.x0:
                scattering_types.holes = vertical_surface_left_scattering(
                    ph, cf.hole_roughness, cf
                )

            # Scattering on the inclined sidewall of the triangle:
            else:
                scattering_types.holes = inclined_surfaces_up_scattering(
                    ph, beta, x, self.x0, cf.hole_roughness
                )

            return True

        # If phonon is inside the left side of the triangle:
        if (
            not self.is_right_half
            and (
                self.size_y / 2 + (y - self.y0)
                <= (self.size_x / 2 - abs(x - self.x0)) / tan(beta)
            )
            and (abs(y - self.y0) < self.size_y / 2)
            and (x < self.x0)
        ):
            # Scattering on the bottom wall of the triangle:
            if (ph.y < self.y0 - self.size_y / 2) and (abs(ph.theta) < pi / 2):
                scattering_types.holes = horizontal_surface_down_scattering(
                    ph, cf.hole_roughness
                )

            # Scattering on the vertical sidewall of the triangle:
            elif ph.x > self.x0:
                scattering_types.holes = vertical_surface_right_scattering(
                    ph, cf.hole_roughness, cf
                )

            # Scattering on the inclined sidewall of the triangle:
            else:
                scattering_types.holes = inclined_surfaces_up_scattering(
                    ph, beta, x, self.x0, cf.hole_roughness
                )

            return True
        return False

    def get_patch(self, color_holes, cf):
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


class ParabolaTop:
    """Shape of a parabolic wall"""

    def __init__(self, tip=0, focus=0):
        self.tip = tip
        self.focus = focus

    def check_if_scattering(self, ph, scattering_types, x, y, z, cf):
        """Scattering on top parabolic boundary"""

        # If phonon is beyond the parabola:
        y_cept = -((cf.width / 2) ** 2) / (4 * self.focus) + self.tip
        if y > y_cept and (x**2 + 4 * self.focus * (y - self.tip)) >= 0:
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

            else:
                scattering_types.walls = Scattering.DIFFUSE
                for _ in range(10):
                    # Lambert distribution
                    ph.theta = normal_theta + asin(2 * random() - 1) - pi / 2
                    ph.phi = asin((asin(2 * random() - 1)) / (pi / 2))

                    # Accept the angles only if they do not immediately cause new scattering:
                    if no_new_scattering(ph, cf):
                        break

            return True
        return False

    def get_patch(self, color_holes, cf):
        eval_xs = linspace(-cf.width / 2, cf.width / 2, 50)
        parabola_ys = -(eval_xs**2) / (4 * self.focus) + self.tip
        parabola_points = column_stack((eval_xs, parabola_ys))
        polygon_point = vstack(
            (parabola_points, [cf.width / 2, cf.length], [-cf.width / 2, cf.length])
        )
        return Polygon(
            polygon_point * 1e6,
            closed=True,
            facecolor=color_holes,
        )


class ParabolaBottom:
    """Shape of a parabolic wall"""

    def __init__(self, tip=0, focus=0):
        self.tip = tip
        self.focus = focus

    def check_if_scattering(self, ph, scattering_types, x, y, z, cf):
        """Scattering on bottom parabolic boundary"""

        # If phonon is below the parabola:
        y_cept = (cf.width / 2) ** 2 / (4 * self.focus) + self.tip
        if y < y_cept and (x**2 - 4 * self.focus * (y - self.tip)) >= 0:
            # Calculate angle to the surface and specular scattering probability:
            normal_theta = pi * (x < 0) - atan(-2 * self.focus / x)
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

            else:
                scattering_types.walls = Scattering.DIFFUSE
                for _ in range(10):
                    # Lambertian distribution
                    ph.theta = normal_theta + asin(2 * random() - 1) - pi / 2
                    ph.phi = asin((asin(2 * random() - 1)) / (pi / 2))

                    # Accept the angles only if they do not immediately cause new scattering:
                    if no_new_scattering(ph, cf):
                        break

            return True
        return False

    def get_patch(self, color_holes, cf):
        eval_xs = linspace(-cf.width / 2, cf.width / 2, 50)
        parabola_ys = eval_xs**2 / (4 * self.focus) + self.tip
        parabola_points = column_stack((eval_xs, parabola_ys))
        polygon_point = vstack((parabola_points, [cf.width / 2, 0], [-cf.width / 2, 0]))
        return Polygon(
            polygon_point * 1e6,
            closed=True,
            facecolor=color_holes,
        )


class CircularPillar:
    """Shape of a circular pillar with inclined wall"""

    def __init__(self, x=0, y=0, diameter=200e-9, height=300e-9, wall_angle=pi / 2):
        self.x0 = x
        self.y0 = y
        self.diameter = diameter
        self.height = height
        self.wall_angle = wall_angle

    def check_if_scattering(self, ph, scattering_types, x, y, z, cf):
        """Check if a phonon strikes a circular pillar and calculate new direction"""

        # Cone radius at a given z coordinate:
        radius = self.diameter / 2  # - (z - cf.thickness / 2) / tan(pillar.wall_angle)
        distance_from_pillar_center = sqrt((x - self.x) ** 2 + (y - self.y) ** 2)
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
            tangent_theta = atan((x - self.x) / (y - self.y))
            scattering_types.pillars = circle_inner_scattering(
                ph, tangent_theta, y, self.y, cf.pillar_roughness
            )

            return True
        return False

    def get_patch(self, color_holes, cf):
        return Circle(
            (1e6 * self.x0, 1e6 * self.y0),
            1e6 * self.diameter / 2,
            facecolor=color_holes,
        )
