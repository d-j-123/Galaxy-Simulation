#Fourth attempt of Barnes - Hut algorithm
#Using flat arrays insted of Python overhead
import numpy as np

G = 1
EPS = 0.01

class BarnesHutTreeV4:
    def __init__(self, boundary, max_nodes):
        self.max_nodes = max_nodes
        self.node_capacity = 4
        self.node_count = 0

        self.x_min = np.empty(max_nodes)
        self.y_min = np.empty(max_nodes)
        self.x_max = np.empty(max_nodes)
        self.y_max = np.empty(max_nodes)

        self.mass = np.zeros(max_nodes)
        self.com_x = np.zeros(max_nodes)
        self.com_y = np.zeros(max_nodes)

        #Children: NW, NE, SW, SE
        self.children = np.full((max_nodes, 4), -1, dtype=np.int32)

        #Particles in leaves
        self.points = [[] for _ in range(max_nodes)]

        #Create root
        self._create_node(boundary.x_min, boundary.y_min, boundary.x_max, boundary.y_max)

    def _create_node(self, x_min, y_min, x_max, y_max):
        if self.node_count >= self.max_nodes:
            raise RuntimeError("Maximum number of nodes exceeded")

        node = self.node_count
        self.node_count += 1

        self.x_min[node] = x_min
        self.y_min[node] = y_min
        self.x_max[node] = x_max
        self.y_max[node] = y_max

        return node

    def subdivide(self, node, positions):
        x_min = self.x_min[node]
        y_min = self.y_min[node]
        x_max = self.x_max[node]
        y_max = self.y_max[node]

        cx = (x_min + x_max) / 2
        cy = (y_min + y_max) / 2

        nw = self._create_node(x_min, y_min, cx, cy)
        ne = self._create_node(cx, y_min, x_max, cy)
        sw = self._create_node(x_min, cy, cx, y_max)
        se = self._create_node(cx, cy, x_max, y_max)

        self.children[node, 0] = nw 
        self.children[node, 1] = ne 
        self.children[node, 2] = sw 
        self.children[node, 3] = se 

        old_points = self.points[node]
        self.points[node] = []

        for index in old_points:
            if self.insert(nw, index, positions): continue
            if self.insert(ne, index, positions): continue
            if self.insert(sw, index, positions): continue
            if self.insert(se, index, positions): continue  

    def insert(self, node, index, positions):
        x_min = self.x_min[node]
        y_min = self.y_min[node]
        x_max = self.x_max[node]
        y_max = self.y_max[node]
        
        contains_point = (x_min <= positions[index, 0] <= x_max and y_min <= positions[index, 1] <= y_max)

        if not contains_point: return False

        if len(self.points[node]) < self.node_capacity and self.children[node, 0] == -1:
            self.points[node].append(index)
            return True

        if self.children[node, 0] == -1:
            self.subdivide(node, positions)

        for child in self.children[node]:
            if self.insert(child, index, positions):
                return True

        return False

    def compute_mass(self, positions, masses):
        def compute(node):
            #Leaf
            if self.children[node, 0] == -1:
                indices = np.asarray(self.points[node], dtype=np.int64)

                if len(indices) == 0:
                    self.mass[node] = 0.0
                    self.com_x[node] = 0.0
                    self.com_y[node] = 0.0
                    
                    return

                self.mass[node] = masses[indices].sum()

                if self.mass[node] == 0.0:
                    self.com_x[node] = 0.0
                    self.com_y[node] = 0.0

                    return

                com = (masses[indices, None] * positions[indices]).sum(axis=0) / self.mass[node]
                self.com_x = com[0]
                self.com_y = com[1]

            #Internal node
            else:
                for child in self.children[node]:
                    compute(child)

                children = self.childre[node]

                self.mass[node] = self.mass[children].sum()

                if self.mass[node] == 0.0:
                    self.com_x[node] = 0.0
                    self.com_y[node] = 0.0

                    return

                self.com_x = (self.mass[children] * self.com_x[children]).sum() / self.mass[node]
                self.com_y = (self.mass[children] * self.com_y[children]).sum() / self.mass[node]

    # def compute_force(self, index, positions, masses, theta):
    #         x = positions[index, 0]
    #         y = positions[index, 1]
    
    #         acceleration = np.array([0.0, 0.0])
    
    #         stack = [self]
    
    #         visited = 0
    #         leaf_nodes = 0
    #         internal_nodes = 0
    
    #         while stack:
    #             node = stack.pop()
    #             visited +=1
    
    #             #Leaf
    #             if node.north_west is None:
    #                 leaf_nodes += 1
    #                 indices = np.asarray(node.points, dtype=np.int64)
    
    #                 #Ignoring the particle itself
    #                 mask = indices != index
    #                 indices = indices[mask]
    
    #                 dx = positions[indices, 0] - x
    #                 dy = positions[indices, 1] - y
    
    #                 d = np.column_stack((dx, dy))
    
    #                 r2 = np.power(dx, 2) + np.power(dy, 2)
    
    #                 acceleration += G * np.sum(masses[indices, None] * d / np.power(r2[:, None] + np.power(EPS, 2), 1.5), axis=0)
    
    #             #Internal node
    #             else:
    #                 internal_nodes += 1
    #                 #s is the width of the region represented by the internal node
    #                 s = node.boundary.x_max - node.boundary.x_min 
    
    #                 dx = node.center_of_mass[0] - x
    #                 dy = node.center_of_mass[1] - y
    
    #                 dd = np.array([dx, dy])
    
    #                 #d is the distance between the body and the node’s center-of-mas
    #                 d2 = np.power(dx, 2) + np.power(dy, 2)
    #                 d = np.sqrt(d2)
    
    #                 #contains_body = node.boundary.contains_point(index, positions)
    #                 contains_body = (node.boundary.x_min <= x <= node.boundary.x_max and node.boundary.y_min <= y <= node.boundary.y_max)
    
    #                 #the internal node is sufficiently far away
    #                 if not contains_body and d > 0 and s / d < theta:
    #                     acceleration += G * node.node_mass * dd / np.power(d2 + np.power(EPS, 2), 1.5)  
    
    #                 else:
    #                     stack.append(node.north_west)
    #                     stack.append(node.north_east)
    #                     stack.append(node.south_west)
    #                     stack.append(node.south_east)
    
    #         return acceleration, visited, leaf_nodes, internal_nodes

    def compute_acceleration(self, index, positions, masses, theta):
        x = positions[index, 0]
        y = positions[index, 1]

        acceleration = np.array([0.0, 0.0])

        stack = [0]

        #visited = 0
        #leaf_nodes = 0
        #internal_nodes = 0

        while stack:
            node = stack.pop()

            #visited += 1

            if self.children[node, 0] == -1:
                #leaf_nodes += 1
                indices = np.asarray(self.points[node], dtype=np.int64)

                mask = indices != index
                indices = indices[mask]

                dx = positions[indices, 0] - x
                dy = positions[indices, 1] - y

                d = np.column_stack((dx, dy))

                r2 = np.power(dx, 2) + np.power(dy, 2)

                acceleration += G * np.sum(masses[indices, None] * d / np.power(r2[:, None] + np.power(EPS, 2), 1.5), axis=0)

            else:
                #internal_nodes += 1
                s = self.x_max[node] - self.x_min[node]

                dx = self.com_x[node] - x
                dy = self.com_y[node] - y

                dd = np.array([dx, dy])

                d2 = dx ** 2 + dy ** 2
                d = np.sqrt(d2)

                contains_body = (self.x_min[node] <= x <= self.x_max[node] and self.y_min[node] <= y <= self.y_max[node])

                if not contains_body and d > 0 and s / d < theta:
                    acceleration += G * self.mass[node] * dd / np.power(d2 + np.power(EPS, 2), 1.5)
                
                else:
                    stack.append(self.children[node, 0])
                    stack.append(self.children[node, 1])
                    stack.append(self.children[node, 2])
                    stack.append(self.children[node, 3])

        return acceleration
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

def create_points(x_max, y_max, n_points):
    positions = np.random.uniform(1, [x_max, y_max], (n_points, 2))
    masses = np.random.randint(1, 10, n_points)

    return positions, masses

def draw(points, screen):
    for x, y in points.astype(np.int32):
        screen.set_at((x, y), (255, 255, 255))

        