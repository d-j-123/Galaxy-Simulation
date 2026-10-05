import pygame
import numpy as np

from objects import Star, BulgeStar, GasCloud, DustCloud

####################
SCREEN_WIDHT = 1900
SCREEN_HEIGHT = 1200

CX = SCREEN_WIDHT // 2
CY = SCREEN_HEIGHT // 2

EPS = 0.1
####################

####################
G = 0.1
T = 5772
AGE = 10e9

MASS_DISK = 5000
MASS_BULGE = 1000

DISK_A = 150
DISK_B = 20

BULGE_A = 20
BULGE_B = 5
BULGE_OMEGA = 0.1

HALO_RHO = 0.01
HALO_R = 300
####################

####################
N_ARMS = 4
GLOBAL_PATTERN_SPEED = 0.1
THETA = 0

SPIRAL_A = 1
SPIRAL_R0 = 1000
SPIRAL_B = -0.5
SPIRAL_SCALE = 1
####################
def create_star(Rd=100, sigma=0.2):
    R = np.random.gamma(2, Rd)

    arm = np.random.randint(0, N_ARMS)

    theta = THETA + np.log(R / SPIRAL_R0) / SPIRAL_B + 2 * np.pi * arm / N_ARMS 
    theta += np.random.normal(0, sigma)

    x = CX + R * np.cos(theta)
    y = CY + R * np.sin(theta)

    vx, vy = orbital_velocity(x, y, radial=True)

    mass = np.random.lognormal(-1.5, 1.3)
    mass = np.clip(mass, 0.08, 15.0)
    luminosity = np.power(mass, 3.5)
    radius = np.power(mass, 0.8)
    temperature = np.power(mass, 0.475) * T
    age = 1 / np.power(mass, 2.5) * AGE

    star = Star(x, y, vx, vy, age, luminosity, mass, temperature, radius)

    return star

def create_bulge_star(bulge_a=BULGE_A, bulge_b=BULGE_B, omega=BULGE_OMEGA):
    R = np.random.gamma(2, 0.3)
    R = min(R, 1)

    theta = np.random.uniform(0, 2 * np.pi)

    x = CX + bulge_a * R * np.cos(theta)
    y = CY + bulge_b * R * np.sin(theta)

    vx, vy = orbital_velocity(x, y, radial=True)

    mass = np.random.lognormal(-1.5, 0.8)
    mass = np.clip(mass, 0.08, 15.0)
    luminosity = np.power(mass, 3.5)
    radius = np.power(mass, 0.8)
    temperature = np.power(mass, 0.475) * T + 2000
    age = AGE / np.power(mass, 2.5)

    bulge_star = BulgeStar(x, y, vx, vy, age, luminosity, mass, temperature, radius, omega, CX, CY)

    return bulge_star

def create_gas_cloud():
    while True:
        R = np.random.gamma(2, 150)

        if R >= 60: break

    arm = np.random.randint(0, N_ARMS)

    theta = THETA + np.log(R / SPIRAL_R0) / SPIRAL_B + 2 * np.pi * arm / N_ARMS 
    theta += np.random.normal(0, 0.1)

    x = CX + R * np.cos(theta)
    y = CY + R * np.sin(theta)
    
    vx, vy = orbital_velocity(x, y, radial=True)

    # size = np.random.lognormal(np.log(10), 0.5)
    # size = np.clip(size, 4, 40)

    size = np.random.lognormal(np.log(10), 0.4)
    size *= 1 + 0.3 * R / 500
    size = np.clip(size, 4, 40)

    density = np.exp(-R / 500)
    density *= np.random.lognormal(0, 0.7)
    density = np.clip(density, 0.05, 3)

    gas_cloud = GasCloud(x, y, vx, vy, size, density)

    return gas_cloud

def create_dust_cloud(gas_cloud):
    dust_cloud = DustCloud(gas_cloud)

    return dust_cloud

def compute_acceleration(star, time):
    # for other in stars:
    #     if star != other:
    #         r2 = (star.position.x - other.position.x) ** 2 + (star.position.y - other.position.y) ** 2
    #         ax += other.mass * (other.position.x - star.position.x) / np.power(r2 + EPS ** 2, 1.5)
    #         ay += other.mass * (other.position.y - star.position.y) / np.power(r2 + EPS ** 2, 1.5)

    #Black hole interaction
    # r2 = (star.position.x - black_hole.position.x) ** 2 + (star.position.y - black_hole.position.y) ** 2
    # ax += black_hole.mass * (black_hole.position.x - star.position.x) / np.power(r2 + EPS ** 2, 1.5)
    # ay += black_hole.mass * (black_hole.position.y - star.position.y) / np.power(r2 + EPS ** 2, 1.5)

    # ax *= G
    # ay *= G

    ax, ay = disk_acceleration(star) + bulge_acceleration(star) + halo_acceleration(star) + spiral_acceleration(star, time)

    star.acceleration = pygame.Vector2(ax, ay) 

def orbital_velocity(x, y, radial=False, sigma=0.1):
    dx = x - CX
    dy = y - CY
    R2 = dx ** 2 + dy ** 2
    R = np.sqrt(R2)

    if R == 0: return 0, 0

    z = 0
    dphi_disk = G * MASS_DISK * R * 1 / np.power(R2 + (DISK_A + np.sqrt(z ** 2 + DISK_B ** 2)) ** 2, 1.5)
    dphi_bulge = G * MASS_BULGE / (R + BULGE_A) ** 2
    dphi_halo = 4 * np.pi * G * HALO_RHO * np.power(HALO_R, 3) / R * (np.log(1 + R / HALO_R) / R - 1 / (R + HALO_R))

    vc = np.sqrt(R * (dphi_disk + dphi_bulge + dphi_halo))
    vc += np.random.normal(0, sigma)

    vr = 0

    if radial:
        vr = np.random.normal(0, sigma)

    return -dy / R * vc + vr * dx / R, dx / R * vc + vr * dy / R

# def disk_potential(x, y, z=0, a=100, b=10):
#     dx = x - CX
#     dy = y - CY
#     R2 = dx ** 2 + dy ** 2
#     R = np.sqrt(R2)
    
#     phi_disk = - G * MASS_DISK / np.sqrt(R2 + (a + np.sqrt(z ** 2 + b ** 2)) ** 2)
#     return phi_disk

def disk_acceleration(star, z=0):
    dx = star.position.x - CX
    dy = star.position.y - CY
    R2 = dx ** 2 + dy ** 2
    R = np.sqrt(R2)

    #Disk acceleration
    #aR = d phi_disk / d R 
    ar = - G * MASS_DISK * R * 1 / np.power(R2 + (DISK_A + np.sqrt(z ** 2 + DISK_B ** 2)) ** 2, 1.5)
    a_phi = 0

    ax = ar * dx / R - a_phi * dy / R
    ay = ar * dy / R + a_phi * dx / R

    #return ax, ay
    return np.array([ax, ay])

# def bulge_potential(x, y, a=0):
#     dx = x - CX
#     dy = y - CY
#     R2 = dx ** 2 + dy ** 2
#     R = np.sqrt(R2)
    
#     #Bulge potential
#     phi_bulge = - G * MASS_BULGE / (R + a)
#     return phi_bulge

def bulge_acceleration(star):
    dx = star.position.x - CX
    dy = star.position.y - CY
    R2 = dx ** 2 + dy ** 2
    R = np.sqrt(R2)

    #Bulge acceleration
    #aR = d phi_bulge / d R
    ar = - G * MASS_BULGE / (R + BULGE_A) ** 2
    a_phi = 0
    
    ax = ar * dx / R - a_phi * dy / R
    ay = ar * dy / R + a_phi * dx / R
    
    #return ax, ay
    return np.array([ax, ay])

# def halo_potential(x, y, rho=1, r=1):
#     dx = x - CX
#     dy = y - CY
#     R2 = dx ** 2 + dy ** 2
#     R = np.sqrt(R2)
    
#     #Halo potential
#     phi_halo = - 4 * np.pi * G * rho * np.power(r, 3) * np.log(1 + R / r) / R
#     return phi_halo

def halo_acceleration(star):
    dx = star.position.x - CX
    dy = star.position.y - CY
    R2 = dx ** 2 + dy ** 2
    R = np.sqrt(R2)

    #Halo acceleration
    #aR = d phi_halo / d R
    ar = - 4 * np.pi * G * HALO_RHO * np.power(HALO_R, 3) / R * (np.log(1 + R / HALO_R) / R - 1 / (R + HALO_R))
    a_phi = 0
        
    ax = ar * dx / R - a_phi * dy / R
    ay = ar * dy / R + a_phi * dx / R
        
    #return ax, ay
    return np.array([ax, ay])

# def spiral_potential(x, y, theta, t=0):
#     dx = x - CX
#     dy = y - CY
#     R2 = dx ** 2 + dy ** 2
#     R = np.sqrt(R2)
    
#     R0 = 1000
#     b = 0.5
#     A = 1
    
#     phi = np.atan2(dy, dx)
#     psi = N_ARMS * (phi - GLOBAL_PATTERN_SPEED * t - theta - np.log(R / R0) / b)
    
#     #Spiral potential
#     phi_spiral = A * np.cos(psi)
#     return phi_spiral

def spiral_acceleration(star, t):
    dx = star.position.x - CX
    dy = star.position.y - CY
    R2 = dx ** 2 + dy ** 2
    R = np.sqrt(R2)

    phi = np.arctan2(dy, dx)
    psi = N_ARMS * (phi - GLOBAL_PATTERN_SPEED * t - THETA - np.log(R / SPIRAL_R0) / SPIRAL_B)

    #Spiral acceleration
    #aR = d phi_spiral / d R
    #aPHI = d phi_spiral / d PHI
    #ar = - SPIRAL_A * N_ARMS / (R * SPIRAL_B) * np.sin(psi)
    #a_phi = SPIRAL_A * N_ARMS * np.sin(psi) / R

    ar = - SPIRAL_A * np.exp(-R / SPIRAL_SCALE) * (- 1 / SPIRAL_SCALE) * np.cos(psi) - SPIRAL_A * np.exp(-R / SPIRAL_SCALE) * N_ARMS / (R * SPIRAL_B) * np.sin(psi)
    a_phi = SPIRAL_A * np.exp(-R / SPIRAL_SCALE) * N_ARMS * np.sin(psi) / R

    ax = ar * dx / R - a_phi * dy / R
    ay = ar * dy / R + a_phi * dx / R
            
    #return ax, ay
    return np.array([ax, ay])

