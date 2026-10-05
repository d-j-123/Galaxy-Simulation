#Wiki pseudocode
class QuadTree:
    def __init__(self, boundary):
        self.node_capacity = 4
        self.points = []
        self.boundary = boundary

        self.north_west = None
        self.north_east = None
        self.south_west = None
        self.south_east = None

    def subdivide(self):
        self.north_west = QuadTree(AABB(self.boundary.x_min, self.boundary.y_min, self.boundary.cx, self.boundary.cy))
        self.north_east = QuadTree(AABB(self.boundary.cx, self.boundary.y_min, self.boundary.x_max, self.boundary.cy))
        self.south_west = QuadTree(AABB(self.boundary.x_min, self.boundary.cy, self.boundary.cx, self.boundary.y_max))
        self.south_east = QuadTree(AABB(self.boundary.cx, self.boundary.cy, self.boundary.x_max, self.boundary.y_max))

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

    def print_values(self):
        print(f"x_min: {self.x_min}")
        print(f"y_min: {self.y_min}")
        print(f"x_max: {self.x_max}")
        print(f"y_max: {self.y_max}")

class Point:
    def __init__(self, x, y):
        self.position = (x, y)
