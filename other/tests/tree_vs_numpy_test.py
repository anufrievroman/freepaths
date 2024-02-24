import timeit
import numpy as np
from scipy.spatial import cKDTree

def time_search_with_ckdtree(points, query_points):
    tree = cKDTree(points)
    return tree.query(query_points)

def time_search_with_numpy(points, query_points):
    distances = np.linalg.norm(points[:, np.newaxis, :] - query_points, axis=2)
    nearest_indices = np.argmin(distances, axis=0)
    nearest_distances = distances[nearest_indices, range(len(query_points))]
    return nearest_distances, nearest_indices

def generate_random_points(num_points, dimension):
    return np.random.rand(num_points, dimension)

# Number of points to test
num_points_list = [10, 100, 1000, 10000]
dimension = 2
num_query_points = 100

for num_points in num_points_list:
    print(f"Number of points: {num_points}")
    points = generate_random_points(num_points, dimension)
    query_points = generate_random_points(num_query_points, dimension)

    time_ckdtree = timeit.timeit(lambda: time_search_with_ckdtree(points, query_points), number=10)
    time_numpy = timeit.timeit(lambda: time_search_with_numpy(points, query_points), number=10)

    print(f"Time taken with cKDTree: {time_ckdtree:.6f} seconds")
    print(f"Time taken with NumPy array: {time_numpy:.6f} seconds")
    print()
