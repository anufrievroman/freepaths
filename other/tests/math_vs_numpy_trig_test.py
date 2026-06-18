"""
Benchmark: math.sin/cos vs numpy.sin/cos for scalar arguments.
In FreePATHS, move() is called on every timestep with plain Python floats,
so the scalar performance is what matters.
"""

import timeit
import random
import math
import numpy as np

N = 1_000_000
theta = random.uniform(-math.pi, math.pi)
phi   = random.uniform(-math.pi / 2, math.pi / 2)

def using_math():
    math.cos(phi)
    math.sin(theta)
    math.cos(theta)

def using_numpy():
    np.cos(phi)
    np.sin(theta)
    np.cos(theta)

t_math  = timeit.timeit(using_math,  number=N)
t_numpy = timeit.timeit(using_numpy, number=N)

print(f"math  (scalar): {t_math:.3f} s  ({N/t_math/1e6:.1f} M calls/s)")
print(f"numpy (scalar): {t_numpy:.3f} s  ({N/t_numpy/1e6:.1f} M calls/s)")
print(f"numpy overhead: {t_numpy/t_math:.1f}x slower")
