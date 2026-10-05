import pygame
from physics import compute_acceleration, create_star, create_bulge_star, create_gas_cloud, create_dust_cloud
from objects import Bulge
import numpy as np

pygame.init()

SCREEN_WIDHT = 1900
SCREEN_HEIGHT = 1200
FPS = 60

STAR_COUNT = 3000
BULGE_STARS = 500
GAS_CLOUDS = 1000

MASS_BULGE = 1000
BULGE_OMEGA = 0.3

T = 5772

screen = pygame.display.set_mode((SCREEN_WIDHT, SCREEN_HEIGHT))
clock = pygame.time.Clock()
font = pygame.font.Font(None, 28)
running = True

pygame.display.set_caption("Galaxy")

stars = []
bulge_stars = []
gas_clouds = []
dust_clouds = []

for _ in range(STAR_COUNT):
    star = create_star()
    stars.append(star)

for _ in range(BULGE_STARS):
    bulge_star = create_bulge_star()
    bulge_stars.append(bulge_star)

bulge = Bulge(bulge_stars, BULGE_OMEGA)

for _ in range(GAS_CLOUDS):
    gas_cloud = create_gas_cloud()
    gas_clouds.append(gas_cloud)

    if np.random.random() < 0.3:
         dust_cloud = create_dust_cloud(gas_cloud)

gas_layer = pygame.Surface((SCREEN_WIDHT, SCREEN_HEIGHT), pygame.SRCALPHA)
dust_layer = pygame.Surface((SCREEN_WIDHT, SCREEN_HEIGHT), pygame.SRCALPHA)

t = 0
while running:
    t_text = font.render(f"t: {t:.0f}", True, "white")

    dt = clock.tick(FPS) / 1000

    fps = clock.get_fps()
    fps_text = font.render(f"FPS: {fps:.0f}", True, "white")

    for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False

    screen.fill("black")
    gas_layer.fill((0, 0, 0, 0))
    dust_layer.fill((0, 0, 0, 0))
########################################################
    for star in stars:
        star.draw(screen)
        compute_acceleration(star, t)
        star.update(dt)

    for gas_cloud in gas_clouds:
        gas_cloud.draw(gas_layer)
        compute_acceleration(gas_cloud, t)
        gas_cloud.update(dt)

    for dust_cloud in dust_clouds:
        dust_cloud.draw(dust_layer)

    bulge.draw(screen)
    bulge.update(dt)

    screen.blit(fps_text, (10, 10))
    screen.blit(t_text, (10, 30))
    screen.blit(gas_layer, (0, 0))
    screen.blit(dust_layer, (0, 0))

    pygame.display.flip()

    t += dt
        
pygame.quit()