#Third attempt of the Barnes - Hut algorithm
#Using Numpy arrays
import numpy as np

G = 1
EPS = 0.01

class Body:
    def __init__(self, x, y, mass):
        self.position = (x, y)
        self.mass = mass
        self.velocity = (0, 0)
        self.acceleration = (0, 0)

class BarnesHutTreeV3:
    def __init__(self, boundary):
        self.node_capacity = 4
        self.points = []
        self.boundary = boundary

        self.north_west = None
        self.north_east = None
        self.south_west = None
        self.south_east = None

        self.node_mass = 0.0
        self.center_of_mass = np.array([0.0, 0.0])

    def subdivide(self, positions):
        self.north_west = BarnesHutTreeV3(AABB(self.boundary.x_min, self.boundary.y_min, self.boundary.cx, self.boundary.cy))
        self.north_east = BarnesHutTreeV3(AABB(self.boundary.cx, self.boundary.y_min, self.boundary.x_max, self.boundary.cy))
        self.south_west = BarnesHutTreeV3(AABB(self.boundary.x_min, self.boundary.cy, self.boundary.cx, self.boundary.y_max))
        self.south_east = BarnesHutTreeV3(AABB(self.boundary.cx, self.boundary.cy, self.boundary.x_max, self.boundary.y_max))

        for index in self.points:
            if self.north_west.insert(index, positions): continue
            if self.north_east.insert(index, positions): continue
            if self.south_west.insert(index, positions): continue
            if self.south_east.insert(index, positions): continue

        self.points = []

    def insert(self, index, positions):
        if not self.boundary.contains_point(index, positions): return False

        if len(self.points) < self.node_capacity and self.north_west is None:
            self.points.append(index)
            return True

        if self.north_west is None:
            self.subdivide(positions)

        if self.north_west.insert(index, positions): return True
        if self.north_east.insert(index, positions): return True
        if self.south_west.insert(index, positions): return True
        if self.south_east.insert(index, positions): return True

        return False   

    def query_range(self, aabb):
        points_in_range = []

        if not self.boundary.intersect(aabb):
            return points_in_range

        for point in self.points:
            if aabb.contains_point(point):
                points_in_range.append(point)

        if self.north_west is None:
            return points_in_range

        points_in_range.extend(self.north_west.query_range(aabb))
        points_in_range.extend(self.north_east.query_range(aabb))
        points_in_range.extend(self.south_west.query_range(aabb))
        points_in_range.extend(self.south_east.query_range(aabb))

        return points_in_range

    def compute_mass(self, positions, masses):
        indices = np.asarray(self.points, dtype=np.int64)

        if self.north_west is None:
            self.node_mass = masses[indices].sum()

            if self.node_mass == 0.0:
                self.center_of_mass = np.array([0.0, 0.0])
                return

            self.center_of_mass = (masses[indices, None] * positions[indices]).sum(axis=0) / self.node_mass

        else:
            children = (self.north_west, self.north_east, self.south_west, self.south_east)

            for child in children:
                child.compute_mass(positions, masses)

            self.node_mass = sum(child.node_mass for child in children)

            if self.node_mass == 0:
                self.center_of_mass = np.array([0.0, 0.0])
                return

            self.center_of_mass = sum(child.node_mass * child.center_of_mass for child in children) / self.node_mass

    def compute_force(self, index, positions, masses, theta):
        x = positions[index, 0]
        y = positions[index, 1]

        acceleration = np.array([0.0, 0.0])

        stack = [self]

        visited = 0
        leaf_nodes = 0
        internal_nodes = 0

        while stack:
            node = stack.pop()
            visited +=1

            #Leaf
            if node.north_west is None:
                leaf_nodes += 1
                indices = np.asarray(node.points, dtype=np.int64)

                #Ignoring the particle itself
                mask = indices != index
                indices = indices[mask]

                dx = positions[indices, 0] - x
                dy = positions[indices, 1] - y

                d = np.column_stack((dx, dy))

                r2 = np.power(dx, 2) + np.power(dy, 2)

                acceleration += G * np.sum(masses[indices, None] * d / np.power(r2[:, None] + np.power(EPS, 2), 1.5), axis=0)

            #Internal node
            else:
                internal_nodes += 1
                #s is the width of the region represented by the internal node
                s = node.boundary.x_max - node.boundary.x_min 

                dx = node.center_of_mass[0] - x
                dy = node.center_of_mass[1] - y

                dd = np.array([dx, dy])

                #d is the distance between the body and the node’s center-of-mas
                d2 = np.power(dx, 2) + np.power(dy, 2)
                d = np.sqrt(d2)

                #contains_body = node.boundary.contains_point(index, positions)
                contains_body = (node.boundary.x_min <= x <= node.boundary.x_max and node.boundary.y_min <= y <= node.boundary.y_max)

                #the internal node is sufficiently far away
                if not contains_body and d > 0 and s / d < theta:
                    acceleration += G * node.node_mass * dd / np.power(d2 + np.power(EPS, 2), 1.5)  

                else:
                    stack.append(node.north_west)
                    stack.append(node.north_east)
                    stack.append(node.south_west)
                    stack.append(node.south_east)

        return acceleration, visited, leaf_nodes, internal_nodes

class AABB:
    def __init__(self, x_min, y_min, x_max, y_max):
        self.x_min = x_min
        self.y_min = y_min
        self.x_max = x_max
        self.y_max = y_max

        self.cx = (self.x_min + self.x_max) / 2
        self.cy = (self.y_min + self.y_max) / 2

    def contains_point(self, index, positions):
        x = positions[index, 0]
        y = positions[index, 1]

        if self.x_min <= x <= self.x_max and self.y_min <= y <= self.y_max: return True

        else: return False

    def intersect(self, aabb):
        if self.x_max < aabb.x_min or self.x_min > aabb.x_max:
            return False

        if self.y_max < aabb.y_min or self.y_min > aabb.y_max:
            return False

        return True

def brute_force(points):
    accelerations = []
    for body in points:
        ax, ay = 0, 0

        for point in points:
            if body is not point:
                dx = point.position[0] - body.position[0]
                dy = point.position[1] - body.position[1]

                r2 = dx**2 + dy**2

                ax += G * point.mass * dx / (r2 + EPS**2)**1.5
                ay += G * point.mass * dy / (r2 + EPS**2)**1.5

        accelerations.append((ax, ay))

    return accelerations

def create_points(x_max, y_max, n_points):
    positions = np.random.uniform(1, [x_max, y_max], (n_points, 2))
    masses = np.random.randint(1, 10, n_points)

    return positions, masses

import pygame
def draw(points, screen):
    for x, y in points.astype(np.int32):
        screen.set_at((x, y), (255, 255, 255))