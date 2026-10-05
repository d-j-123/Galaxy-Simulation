#First attempt of the Barnes - Hut algorithm
#Naive implementation of the algorithm
import numpy as np

G = 1
EPS = 0.01

class Body:
    def __init__(self, x, y, mass):
        self.position = (x, y)
        self.mass = mass
        self.velocity = (0, 0)
        self.acceleration = (0, 0)

class BarnesHutTreeV1:
    def __init__(self, boundary):
        self.node_capacity = 4
        self.points = []
        self.boundary = boundary

        self.north_west = None
        self.north_east = None
        self.south_west = None
        self.south_east = None

        self.node_mass = 0
        self.center_of_mass = (0, 0)

    def subdivide(self):
        self.north_west = BarnesHutTreeV1(AABB(self.boundary.x_min, self.boundary.y_min, self.boundary.cx, self.boundary.cy))
        self.north_east = BarnesHutTreeV1(AABB(self.boundary.cx, self.boundary.y_min, self.boundary.x_max, self.boundary.cy))
        self.south_west = BarnesHutTreeV1(AABB(self.boundary.x_min, self.boundary.cy, self.boundary.cx, self.boundary.y_max))
        self.south_east = BarnesHutTreeV1(AABB(self.boundary.cx, self.boundary.cy, self.boundary.x_max, self.boundary.y_max))

        for point in self.points:
            # self.north_west.insert(point)
            # self.north_east.insert(point)
            # self.south_west.insert(point)
            # self.south_east.insert(point)

            if self.north_west.insert(point): continue
            if self.north_east.insert(point): continue
            if self.south_west.insert(point): continue
            if self.south_east.insert(point): continue

        self.points = []

    def insert(self, point):
        if not self.boundary.contains_point(point):
            return False

        if len(self.points) < self.node_capacity and self.north_west is None:
            self.points.append(point)
            return True

        if self.north_west is None:
            self.subdivide()

        if self.north_west.insert(point):
            return True

        if self.north_east.insert(point):
            return True

        if self.south_west.insert(point):
            return True

        if self.south_east.insert(point):
            return True

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

    def compute_mass(self):
        if self.north_west is None:
            self.node_mass = sum(point.mass for point in self.points)

            if self.node_mass == 0:
                self.center_of_mass = (0, 0)
                return

            cx = sum(point.mass * point.position[0] for point in self.points) / self.node_mass
            cy = sum(point.mass * point.position[1] for point in self.points) / self.node_mass

            self.center_of_mass = (cx, cy)

        else:
            self.north_west.compute_mass()
            self.north_east.compute_mass()
            self.south_west.compute_mass()
            self.south_east.compute_mass()

            self.node_mass = self.north_west.node_mass + self.north_east.node_mass + self.south_west.node_mass + self.south_east.node_mass

            if self.node_mass == 0:
                self.center_of_mass = (0, 0)
                return

            cx = (self.north_west.node_mass * self.north_west.center_of_mass[0] + self.north_east.node_mass * self.north_east.center_of_mass[0] + 
                    self.south_west.node_mass * self.south_west.center_of_mass[0] + self.south_east.node_mass * self.south_east.center_of_mass[0]) / self.node_mass

            cy = (self.north_west.node_mass * self.north_west.center_of_mass[1] + self.north_east.node_mass * self.north_east.center_of_mass[1] + 
                    self.south_west.node_mass * self.south_west.center_of_mass[1] + self.south_east.node_mass * self.south_east.center_of_mass[1]) / self.node_mass

            self.center_of_mass = (cx, cy)

    def compute_force(self, body, theta):
        ax, ay = 0, 0
        #Leaf
        if self.north_west is None:
            for point in self.points:
                if point is not body:
                    dx = point.position[0] - body.position[0]
                    dy = point.position[1] - body.position[1]

                    r2 = dx ** 2 + dy ** 2
                    ax += G * point.mass * dx / (r2 + EPS ** 2) ** 1.5
                    ay += G * point.mass * dy / (r2 + EPS ** 2) ** 1.5

        else:
            #s is the width of the region represented by the internal node
            s = self.boundary.x_max - self.boundary.x_min 
            
            dx = self.center_of_mass[0] - body.position[0]
            dy = self.center_of_mass[1] - body.position[1]

            #d is the distance between the body and the node’s center-of-mass
            d2 = dx ** 2 + dy ** 2 
            d = d2 ** 0.5

            contains_body = self.boundary.contains_point(body)
            #the internal node is sufficiently far away
            if not contains_body and d > 0 and s / d < theta:
                ax += G * self.node_mass * dx /(d2 + EPS ** 2) ** 1.5
                ay += G * self.node_mass * dy /(d2 + EPS ** 2) ** 1.5

            else:
                ax1, ay1 = self.north_west.compute_force(body, theta)
                ax2, ay2 = self.north_east.compute_force(body, theta)
                ax3, ay3 = self.south_west.compute_force(body, theta)
                ax4, ay4 = self.south_east.compute_force(body, theta)

                ax = ax1 + ax2 + ax3 + ax4
                ay = ay1 + ay2 + ay3 + ay4

        return ax, ay

class AABB:
    def __init__(self, x_min, y_min, x_max, y_max):
        self.x_min = x_min
        self.y_min = y_min
        self.x_max = x_max
        self.y_max = y_max

        self.cx = (self.x_min + self.x_max) / 2
        self.cy = (self.y_min + self.y_max) / 2

    def contains_point(self, point):
        x, y = point.position

        if self.x_min <= x <= self.x_max and self.y_min <= y <= self.y_max:
            return True
        
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
    points = []

    for _ in range(n_points):
        x = np.random.uniform(0, x_max)
        y = np.random.uniform(0, y_max)
        mass = np.random.randint(1, 10)

        point = Body(x, y, mass)
        points.append(point)

    return points