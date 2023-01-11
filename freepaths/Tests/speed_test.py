import numpy as np
import math
from math import cos
import timeit
import random
from functools import lru_cache

phi = random.random()

def math_lib():
    a = []
    for i in range(1000):
        a.append(phi)
    return

def numpy_lib():
    a = np.zeros((1000))
    for i in range(1000):
        a[i] = phi
    return

n = 10000

t1 = timeit.timeit(math_lib, number = n)
print ('Math:', t1)
t2 = timeit.timeit(numpy_lib, number = n)
print ('Numpy:', t2)
