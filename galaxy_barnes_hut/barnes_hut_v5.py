#Fourth attempt of Barnes - Hut algorithm
#Using flat arrays insted of Python overhead
import numpy as np

G = 0.1
EPS = 1

GALAXY_MASS = 100000
HALO_MASS = 0.75 * GALAXY_MASS
MASS_DISK = 0.2 * GALAXY_MASS
MASS_BULGE = 0.05 * GALAXY_MASS
BH_MASS = 0.0001 * GALAXY_MASS

DISK_A = 150
DISK_B = 20

BULGE_A = 20
BULGE_B = 5
BULGE_OMEGA = 0.1

HALO_RHO = 0.01
HALO_R = 300
A = 100

from numba import njit, prange

@njit
def compute_acceleration_numba(index, positions, masses, theta, x_min, y_min, x_max, y_max, node_mass, com_x, com_y, children, points, point_count):
    x = positions[index, 0]
    y = positions[index, 1]

    acceleration_x = 0.0
    acceleration_y = 0.0

    stack = [0]

    while stack:
        node = stack.pop()

        if children[node, 0] == -1:
        #leaf_nodes += 1
            count = point_count[node]
            for i in range(count):
                point = points[node, i]

                if point == index: continue
            
                dx = positions[point, 0] - x
                dy = positions[point, 1] - y
        
        
                r2 = dx * dx + dy * dy
        
                acceleration = G * masses[point] / (r2 + EPS * EPS) ** 1.5

                acceleration_x += acceleration * dx
                acceleration_y += acceleration * dy

        else:
            s = x_max[node] - x_min[node]

            dx = com_x[node] - x
            dy = com_y[node] - y

            d2 = dx * dx + dy * dy
            d = np.sqrt(d2)

            contains_body = (x_min[node] <= x <= x_max[node] and y_min[node] <= y <= y_max[node])

            if not contains_body and d > 0.0 and s / d < theta:
                acceleration = G * node_mass[node] / (d2 + EPS * EPS) ** 1.5

                acceleration_x += acceleration * dx
                acceleration_y += acceleration * dy

            else:
                stack.append(children[node, 0])
                stack.append(children[node, 1])
                stack.append(children[node, 2])
                stack.append(children[node, 3])

    return acceleration_x, acceleration_y

@njit
def insert_numba(node, index, positions, x_min, y_min, x_max, y_max, children, points, point_count, node_count):
    x = positions[index, 0]
    y = positions[index, 1]

    contains_point = (x_min[node] <= x <= x_max[node] and y_min[node] <= y <= y_max[node])
    
    if not contains_point: return node_count, False

    if point_count[node] < 4 and children[node, 0] == -1:
        count = point_count[node]
        points[node, count] = index
        point_count[node] += 1
    
        return node_count, True
    
    if children[node, 0] == -1:
        #Subdivide
        new_index = index

        cx = (x_min[node] + x_max[node]) / 2
        cy = (y_min[node] + y_max[node]) / 2

        #NW
        nw = node_count
        node_count += 1

        x_min[nw] = x_min[node]
        y_min[nw] = y_min[node]
        x_max[nw] = cx
        y_max[nw] = cy

        #NE
        ne = node_count
        node_count += 1

        x_min[ne] = cx
        y_min[ne] = y_min[node]
        x_max[ne] = x_max[node]
        y_max[ne] = cy

        #SW
        sw = node_count
        node_count += 1

        x_min[sw] = x_min[node]
        y_min[sw] = cy
        x_max[sw] = cx
        y_max[sw] = y_max[node]

        #SE
        se = node_count 
        node_count += 1

        x_min[se] = cx
        y_min[se] = cy
        x_max[se] = x_max[node]
        y_max[se] = y_max[node]
        
        children[node, 0] = nw 
        children[node, 1] = ne 
        children[node, 2] = sw 
        children[node, 3] = se 
        
        old_count = point_count[node]
        
        for i in range(old_count):
            old_index = points[node, i]

            node_count, inserted = insert_numba(nw, old_index, positions, x_min, y_min, x_max, y_max, children,points, point_count, node_count)
            if inserted: continue

            node_count, inserted = insert_numba(ne, old_index, positions, x_min, y_min, x_max, y_max, children,points, point_count, node_count)
            if inserted: continue

            node_count, inserted = insert_numba(sw, old_index, positions, x_min, y_min, x_max, y_max, children,points, point_count, node_count)
            if inserted: continue

            node_count, inserted = insert_numba(se, old_index, positions, x_min, y_min, x_max, y_max, children,points, point_count, node_count)
            if inserted: continue
        
        point_count[node] = 0

        index = new_index

    for i in range(4):
        child = children[node, i]

        node_count, inserted = insert_numba(child, index, positions, x_min, y_min, x_max, y_max, children,points, point_count, node_count)

        if inserted:
            return node_count, True

    return node_count, False

@njit
def compute_mass_numba(node, positions, masses, children, points, point_count, node_mass, com_x, com_y):
    #Leaf
    if children[node, 0] == -1:
        count = point_count[node]
    
        if count == 0:
            node_mass[node] = 0.0
            com_x[node] = 0.0
            com_y[node] = 0.0
                        
            return
    
        total_mass = 0.0
        mass_x = 0.0
        mass_y = 0.0

        for i in range(count):
            index = points[node, i]
            mass = masses[index]

            total_mass += mass
            mass_x += mass * positions[index, 0]
            mass_y += mass * positions[index, 1]

        node_mass[node] = total_mass
    
        if total_mass == 0.0:
            com_x[node] = 0.0
            com_y[node] = 0.0
    
            return
    
        com_x[node] = mass_x / total_mass
        com_y[node] = mass_y / total_mass

        return
    
    #Internal node
    else:
        total_mass = 0.0
        mass_x = 0.0
        mass_y = 0.0

        for i in range(4):
            child = children[node, i]
    
            compute_mass_numba(child, positions, masses, children, points, point_count, node_mass, com_x, com_y)

            child_mass = node_mass[child]

            total_mass += child_mass
            mass_x += child_mass * com_x[child]
            mass_y += child_mass * com_y[child]

        node_mass[node] = total_mass

        if total_mass == 0.0:
            com_x[node] = 0.0
            com_y[node] = 0.0
    
            return
    
        com_x[node] = mass_x / total_mass
        com_y[node] = mass_y / total_mass

@njit(parallel=True)
def compute_accelerations_numba(n, positions, masses, theta, x_min, y_min, x_max, y_max, node_mass, com_x, com_y, children, points, point_count, cx, cy):
    accelerations = np.zeros_like(positions)

    for index in prange(n):
        x = positions[index, 0]
        y = positions[index, 1]
        
        acceleration_x = 0.0
        acceleration_y = 0.0
        
        stack = [0]
        
        while stack:
            node = stack.pop()
        
            if children[node, 0] == -1:
            #leaf_nodes += 1
                count = point_count[node]
                for i in prange(count):
                    point = points[node, i]
        
                    if point == index: continue
                    
                    dx = positions[point, 0] - x
                    dy = positions[point, 1] - y
                
                
                    r2 = dx * dx + dy * dy
                
                    acceleration = G * masses[point] / (r2 + EPS * EPS) ** 1.5
        
                    acceleration_x += acceleration * dx
                    acceleration_y += acceleration * dy
        
            else:
                s = x_max[node] - x_min[node]
        
                dx = com_x[node] - x
                dy = com_y[node] - y
        
                d2 = dx * dx + dy * dy
                d = np.sqrt(d2)
        
                contains_body = (x_min[node] <= x <= x_max[node] and y_min[node] <= y <= y_max[node])
        
                if not contains_body and d > 0.0 and s / d < theta:
                    acceleration = G * node_mass[node] / (d2 + EPS * EPS) ** 1.5
        
                    acceleration_x += acceleration * dx
                    acceleration_y += acceleration * dy
        
                else:
                    stack.append(children[node, 0])
                    stack.append(children[node, 1])
                    stack.append(children[node, 2])
                    stack.append(children[node, 3])

        accelerations[index, 0] = acceleration_x
        accelerations[index, 1] = acceleration_y

    ########
    #Halo acceleration
    halo_dx = cx - positions[:, 0]
    halo_dy = cy - positions[:, 1]

    r2 = halo_dx * halo_dx + halo_dy * halo_dy
    halo_acc = G * HALO_MASS / (r2 + A * A) ** 1.5

    halo_x = halo_acc * halo_dx
    halo_y = halo_acc * halo_dy

    accelerations[:, 0] += halo_x
    accelerations[:, 1] += halo_y

    return accelerations[:-1, :]

@njit
def build_tree_numba(positions, x_min, y_min, x_max, y_max, children, points, point_count, node_count, scree_w, screen_h):
    x_min[0] = 0
    y_min[0] = 0
    x_max[0] = scree_w
    y_max[0] = screen_h
    
    
    node_count = 1
    
    #children[:, :] = -1
    #point_count[:] = 0

    children[0, 0] = -1
    children[0, 1] = -1
    children[0, 2] = -1
    children[0, 3] = -1

    point_count[0] = 0

    n = positions.shape[0]

    for i in range(n):
        node_count, inserted = insert_numba(0, i, positions, x_min, y_min, x_max, y_max, children, points, point_count, node_count)

    return node_count

def create_points(x_max, y_max, n_points, R_max = 550):
    masses = np.random.randint(1, 10, n_points)

    R = np.random.uniform(20, R_max, n_points) + np.random.normal(0, 0.2, n_points)
    theta = np.random.uniform(0, 2 * np.pi, n_points)
    theta += np.random.normal(0, 0.1, n_points)

    cx = np.full(n_points, x_max / 2)
    cy = np.full(n_points, y_max / 2)

    x = cx + R * np.cos(theta)
    y = cy + R * np.sin(theta)

    positions = np.column_stack((x, y))

    return positions, masses

def create_disk(n_points, R_max, x_max, y_max):
    #masses = np.random.uniform(1, 10, n_points)
    masses = np.random.lognormal(0, 0.5,n_points)
    masses = np.clip(masses, 0.08, 15)
    masses *= MASS_DISK / masses.sum()

    R = np.random.uniform(20, R_max, n_points)
    theta = np.random.uniform(0, 2 * np.pi, n_points)
    theta += np.random.normal(0, 0.1, n_points)

    cx = np.full(n_points, x_max / 2)
    cy = np.full(n_points, y_max / 2)
    
    x = cx + R * np.cos(theta)
    y = cy + R * np.sin(theta)
    
    positions = np.column_stack((x, y))

    return positions, masses

def create_bulge(n_points, R_max, x_max, y_max):
    #masses = np.random.uniform(20, 40, n_points)
    masses = np.random.lognormal(0, 0.5,n_points)
    masses = np.clip(masses, 0.08, 15)
    masses *= MASS_BULGE / masses.sum()
    
    R = np.random.uniform(1, R_max, n_points)
    theta = np.random.uniform(0, 2 * np.pi, n_points)
    theta += np.random.normal(0, 0.1, n_points)
    
    cx = np.full(n_points, x_max / 2)
    cy = np.full(n_points, y_max / 2)
        
    x = cx + R * np.cos(theta)
    y = cy + R * np.sin(theta)
        
    positions = np.column_stack((x, y))
    
    return positions, masses

def create_black_hole(cx, cy):
    position = np.array([cx, cy])
    mass = np.array([BH_MASS])

    return position, mass

def draw(points, types, screen):
    for i in range(len(points)):
        x = int(points[i, 0])
        y = int(points[i, 1])

        if types[i] == 0:
            color = (255, 255, 255)
        elif types[i] == 1:
            color = (255, 180, 80)
        else:
            color = (0, 0, 255)
        screen.set_at((x, y), color)
        