from barnes_hut_v1 import Body, BarnesHutTreeV1, brute_force, create_points
from barnes_hut_v2 import BarnesHutTreeV2
from barnes_hut_v3 import BarnesHutTreeV3, create_points, AABB
import numpy as np
import time

x_min = 0
y_min = 0
x_max = 1000
y_max = 1000

theta = 0.5

n_points = 10000

# A = Body(3, 2, 1)
# B = Body(14, 2, 2)
# C = Body(13,3, 3)
# D = Body(18, 8, 4)
# E = Body(7,13, 5)
# F = Body(2,16, 6)
# G = Body(6, 18, 7)
# H = Body(12, 18, 8)

# points = [A, B, C, D, E, F, G, H]
#points = create_points(x_max, y_max, n_points)
positions, masses = create_points(x_max, y_max, n_points)

boundary = AABB(x_min, y_min, x_max, y_max)

start = time.perf_counter()

qt = BarnesHutTreeV3(boundary)

# for point in points:
#     qt.insert(point)

# qt.compute_mass()

for i in range(positions.shape[0]):
    qt.insert(i, positions)

qt.compute_mass(positions, masses)

tree_time = time.perf_counter() - start

bh_accelerations = []

start = time.perf_counter()

# for point in points:
#     bh_accelerations.append(qt.compute_force(point, theta))

for i in range(positions.shape[0]):
    bh_accelerations.append(qt.compute_force(i, positions, masses, theta))

force_time = time.perf_counter() - start

print(f"N={n_points}, x_max={x_max}, y_max={y_max}")
print(f"Theta={theta}")
print(f"Tree time: {tree_time}")
print(f"Force time: {force_time}")

# #print(f"BH: {bh_accelerations}")
# #print(f"Exact: {exact}")

# errors = []

# for i, body in enumerate(points):
#     ax_exact, ay_exact = exact[i]
#     ax_bh, ay_bh = bh_accelerations[i]

#     error = ((ax_bh - ax_exact) ** 2 + (ay_bh - ay_exact) ** 2) ** 0.5 #/ (ax_exact ** 2 + ay_exact ** 2) ** 0.5
#     errors.append(error)
#     #print(f"Error: {error}")

# print(f"Mean: {np.mean(errors)}")
# print(f"Max: {np.max(errors)}")
