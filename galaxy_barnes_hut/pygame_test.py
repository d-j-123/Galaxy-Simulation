import pygame
import numpy as np

from barnes_hut_v5 import create_points, draw, compute_mass_numba, compute_accelerations_numba, build_tree_numba, create_bulge, create_disk, create_black_hole

SCREEN_WIDHT = 1900
SCREEN_HEIGHT = 1200
FPS = 60

CX = SCREEN_WIDHT / 2
CY = SCREEN_HEIGHT / 2

N = 50000
N_DISK = int(0.9 * N)
N_BULGE = int(0.1 * N)
THETA = 0.5

pygame.init()

screen = pygame.display.set_mode((SCREEN_WIDHT, SCREEN_HEIGHT))
font = pygame.font.Font(None, 28)
running = True

pygame.display.set_caption("Galaxy")

disk_positions, disk_masses = create_disk(N_DISK, 550, x_max=SCREEN_WIDHT, y_max=SCREEN_HEIGHT)
bulge_positions, bulge_masses = create_bulge(N_BULGE, 20, x_max=SCREEN_WIDHT, y_max=SCREEN_HEIGHT)
black_hole_position, black_hole_mass = create_black_hole(CX, CY)

positions = np.vstack((disk_positions, bulge_positions, black_hole_position))
masses = np.concatenate((disk_masses, bulge_masses, black_hole_mass))
velocities = np.zeros((positions.shape[0] - 1, 2))

types = np.concatenate((np.zeros(N_DISK, dtype=np.int8), np.ones(N_BULGE, dtype=np.int8), np.full(1, 2, dtype=np.int8)))

#############################
max_nodes = 4 * N + 1

x_min = np.empty(max_nodes)
y_min = np.empty(max_nodes)
x_max = np.empty(max_nodes)
y_max = np.empty(max_nodes)

node_mass = np.zeros(max_nodes)
com_x = np.zeros(max_nodes)
com_y = np.zeros(max_nodes)

children = np.full((max_nodes, 4), -1, dtype=np.int32)

points = np.full((max_nodes, 4), -1, dtype=np.int32)
point_count = np.zeros(max_nodes, dtype=np.int32)

#Warm-up
node_count = 1

x_min[0] = 0
y_min[0] = 0
x_max[0] = SCREEN_WIDHT
y_max[0] = SCREEN_HEIGHT

children[:, :] = -1
point_count[:] = 0

dx = positions[:-1, 0] - SCREEN_WIDHT / 2
dy = positions[:-1, 1] - SCREEN_HEIGHT / 2

node_count = build_tree_numba(positions, x_min, y_min, x_max, y_max, children, points, point_count, node_count, SCREEN_WIDHT, SCREEN_HEIGHT)
compute_mass_numba(0, positions, masses, children, points, point_count, node_mass, com_x, com_y)
init_accelerations = compute_accelerations_numba(N, positions, masses, THETA, x_min, y_min, x_max, y_max, node_mass, com_x, com_y, children, points, point_count, CX, CY)

R = np.sqrt(dx * dx + dy * dy)

rx = dx / R
ry = dy / R

a_r = init_accelerations[:, 0] * rx + init_accelerations[:, 1] * ry

vc = np.sqrt(np.maximum(0.0, -R * a_r))
vr = np.random.normal(0, 0.1, N)

velocities[:, 0] = -ry * vc + vr * dx / R 
velocities[:, 1] = rx * vc + vr * dy / R 

clock = pygame.time.Clock()
t = 0
#########################
while running:
    t_text = font.render(f"t: {t:.0f}", True, "white")

    dt = clock.tick(FPS) / 1000

    fps = clock.get_fps()
    fps_text = font.render(f"FPS: {fps:.0f}", True, "white")

    for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False

    screen.fill("black")
###########################
    draw(positions, types, screen)

    #Reset
    x_min[0] = 0
    y_min[0] = 0
    x_max[0] = SCREEN_WIDHT
    y_max[0] = SCREEN_HEIGHT

    node_count = 1

    children[:, :] = -1
    point_count[:] = 0

    node_count = build_tree_numba(positions, x_min, y_min, x_max, y_max, children, points, point_count, node_count, SCREEN_WIDHT, SCREEN_HEIGHT)

    compute_mass_numba(0, positions, masses, children, points, point_count, node_mass, com_x, com_y)

    accelerations = compute_accelerations_numba(N, positions, masses, THETA, x_min, y_min, x_max, y_max, node_mass, com_x, com_y, children, points, point_count, CX, CY)

    velocities += accelerations * dt
    positions[:-1, :] += velocities * dt

    screen.blit(fps_text, (10, 10))
    screen.blit(t_text, (10, 30))

    pygame.display.flip()

    t += dt

pygame.quit()