from numpy import pi, cos, sin, tan, exp, arctan, arcsin, sign, log, sqrt, arccos
from random import random
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as axes3d

fig = plt.figure()
ax = fig.add_subplot(1,1,1, projection='3d')

for n in range(1000):

    ##############################################
    # INTERNAL SCATTERING:

    # Random distribution:
    # theta = -pi + random()*2*pi
    # phi = arcsin(2*random()-1)

    ##############################################
    # INITIALIZATION:

    # Random distribution:
    # theta = -pi/2 + pi*random()
    # phi = arcsin(2*random()-1)

    # Lambert cosine distribution:
    # theta = arcsin(2*random()-1)
    # phi=arcsin((arcsin(2*random()-1))/(pi/2))

    ##############################################
    # BOTTOM SCATTERING:

    # Random distribution:
    # theta = -pi + random()*2*pi
    # phi = arcsin(random())

    # Lambert cosine distribution (correct):
    theta = -pi + random()*2*pi
    phi = random()*pi/2

    ##############################################
    # TOP SCATTERING:

    # Random distribution:
    # theta = -pi + random()*2*pi
    # phi = -arcsin(random())

    # Lambert cosine distribution (correct):
    # theta = -pi + random()*2*pi
    # phi = -random()*pi/2

    ##############################################

    # Random distribution:
    # x = 1
    # theta=-sign(x)*pi+sign(x)*pi*random()
    # phi=arcsin(2*random()-1)

    # Lambert cosine distribution (correct):
    # x = 1
    # theta = -sign(x)*pi/2 + arcsin(2*random()-1)
    # phi = arcsin((arcsin(2*random() - 1))/(pi/2))

    ############### RECTANGULAR HOLES ###############

    # LOWER WALL OF SQUARE HOLES:
    # Lambert cosine distribution:
    # rand_sign = sign((2*random() - 1))
    # theta = rand_sign*pi/2 + rand_sign*arccos(random())
    # phi = arcsin((arcsin(2*random() - 1))/(pi/2))

    # UPPER WALL OF SQUARE HOLES:
    # Lambert cosine distribution:
    # theta = arcsin(2*random() - 1)
    # phi = arcsin((arcsin(2*random() - 1))/(pi/2))

    # SIDEWALLS OF SQUARE HOLES:
    # original_theta = -pi/4 # Just to know which wall
    # theta = -sign(sin(original_theta))*pi/2 + arcsin(2*random()-1)
    # phi = arcsin((arcsin(2*random() - 1))/(pi/2))


    ############### CIRCULAR HOLES ##################

    # Parameters of the hole and phonon:
    y = 1
    y0 = 0
    x = 1
    x0 = 0
    if y == y0: y += 1e-9 # Prevent division by zero
    tangent_theta = arctan((x - x0)/(y - y0))

    # WALL OF CIRCULAR HOLES:

    # Random distribution:
    # theta = tangent_theta - (-pi/2 + pi*random()) - pi*(y < y0)
    # phi = arcsin(2*random() - 1)

    # Lambert cosine distribution:
    # theta = tangent_theta - (arcsin(2*random() - 1)) - pi*(y < y0)
    # phi = arcsin((arcsin(2*random() - 1))/(pi/2))

    ############## TRIANGULAR HOLES ################

    # Parameters of the hole and phonon:
    Lx = 1
    Ly = 1
    beta = arctan(0.5*Lx/Ly)
    x = 1
    x0 = 0

    # BOTTOM OF THE TRIANGLE FACING UP:

    # Lambert cosine distribution:
    # rand_sign = sign((2*random() - 1))
    # theta = rand_sign*pi/2 + rand_sign*arccos(random())
    # phi = arcsin((arcsin(2*random() - 1))/(pi/2))

    # SIDEWALLS OF THE TRIANGLE FACING UP:

    # Lambert cosine distribution:
    # theta=arcsin(2*random() - 1) + sign(x - x0)*(pi/2 - beta)
    # phi = arcsin((arcsin(2*random() - 1))/(pi/2))

    # TOP OF THE TRIANGLE FACING DOWN:
    # Lambert cosine distribution:
    # theta = arcsin(2*random() - 1)
    # phi = arcsin((arcsin(2*random() - 1))/(pi/2))

    # SIDEWALLS OF THE TRIANGLE FACING DOWN:

    # Lambert cosine distribution:
    # rand_sign = sign((2*random() - 1))
    # theta = rand_sign*pi - rand_sign*arcsin(random()) - sign(x-x0)*(pi/2 - beta)
    # phi = arcsin((arcsin(2*random() - 1))/(pi/2))

    ##############################################
    # PLOTTING THE DISTRIBUTION:
    X = sin(theta)*abs(cos(phi))
    Y = cos(theta)*abs(cos(phi))
    Z = sin(phi)
    plt.plot(X, Y, Z, '.', c='k', markersize='1')
    ax.set_xlabel('X')
    ax.set_xlim3d(-1, 1)
    ax.set_ylabel('Y')
    ax.set_ylim3d(-1, 1)
    ax.set_zlabel('Z')
    ax.set_zlim3d(-1, 1)
plt.show()
