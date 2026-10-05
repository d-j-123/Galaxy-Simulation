from quadtree import Point, QuadTree, AABB

x_min = 0
y_min = 0
x_max = 20
y_max = 20

A = Point(3, 2)
B = Point(14, 2)
C = Point(13,3)
D = Point(18, 8)
E = Point(7,13)
F = Point(2,16)
G = Point(6, 18)
H = Point(12, 18)

points = [A, B, C, D, E, F]

boundary = AABB(x_min, y_min, x_max, y_max)

qt = QuadTree(boundary)

for point in points:
    qt.insert(point)

range = AABB(0, 0, 10, 10)

result = qt.query_range(range)
print(result)
