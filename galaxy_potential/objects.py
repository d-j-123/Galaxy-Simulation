import pygame
import numpy as np

class Star:
    def __init__(self, x, y, vx, vy, age, luminosity, mass, temperature, radius):
        self.position = pygame.Vector2(x, y)
        self.velocity = pygame.Vector2(vx, vy)
        self.acceleration = pygame.Vector2(0, 0)
        self.age = age
        self.luminosity = luminosity
        self.brightness = np.clip(140 + 35 * np.log10(self.luminosity), 70, 255)
        self.mass = mass
        self.temperature = temperature
        self.radius = radius
        self.color = temperature_to_color(self.temperature)

    def update(self, dt):
        self.velocity += self.acceleration * dt
        self.position += self.velocity * dt

    def draw(self, screen):
        pygame.draw.circle(screen, self.color, self.position, draw_radius(self.radius))

class Bulge:
    def __init__(self, bulge_stars, omega):
        self.theta = 0
        self.stars = bulge_stars
        self.omega = omega

    def update(self, dt):
        self.theta += self.omega * dt

        for star in self.stars:
            star.update(self.theta)

    def draw(self, screen):
        for star in self.stars:
            star.draw(screen)

class BulgeStar(Star):
    def __init__(self, x, y, vx, vy, age, luminosity, mass, temperature, radius, omega, CX, CY):
        super().__init__(x, y, vx, vy, age, luminosity, mass, temperature, radius)

        self.omega = omega
        self.bulge_center = pygame.Vector2(CX, CY)

        self.local_position = pygame.Vector2(x - CX, y - CY)

    def update(self, theta):
        #x = self.position.x - self.bulge_center.x
        #y = self.position.y - self.bulge_center.y

        self.position.x = self.local_position.x * np.cos(theta) - self.local_position.y * np.sin(theta) + self.bulge_center.x
        self.position.y = self.local_position.x * np.sin(theta) + self.local_position.y * np.cos(theta) + self.bulge_center.y
    
def temperature_to_color(T):
    #Red
    if T < 3500:
        return (255, 137, 18) 
    #Orange  
    elif T > 3500 and T < 5000:
        return (255, 190, 126) 
    #Yellow      
    elif T > 5000 and T < 6000:
        return (255, 230, 199)
    #Yellow-white     
    elif T < 6000 and T < 7500:
        return (255, 244, 234)    
    #White  
    elif T < 7500 and T < 10000:
        return (248, 247, 255)      
    #Blue-white
    elif T < 10000 and T < 30000:
        return (202, 215, 255)  
    #Blue    
    else:
        return (155, 176, 255)      
     

def apply_brightness(color, brightness):
    factor = brightness / 255

    return tuple(min(255, int(c * factor)) for c in color)

def draw_radius(radius):
    return int(np.clip(0.8 + 0.4 * radius, 1, 4))

class GasCloud:
    def __init__(self, x, y, vx, vy, size, density):
        self.position = pygame.Vector2(x, y)
        self.velocity = pygame.Vector2(vx, vy)
        self.acceleration = pygame.Vector2(0, 0)
        self.size = size
        self.density = density

        self.alpha = int(np.clip(self.density * 25, 3, 60))

    def update(self, dt):
        self.velocity += self.acceleration * dt
        self.position += self.velocity * dt

    def draw(self, screen):
        pygame.draw.circle(screen, (100, 150, 255, self.alpha), self.position, self.size)

# class DustCloud:
#     def __init__(self, gas_cloud):
#         self.gas_cloud = gas_cloud

#         angle = np.random.uniform(0, 2 * np.pi)
#         distance = np.random.uniform(0, self.gas_cloud.size)

#         self.offset = pygame.Vector2(np.cos(angle) * distance, np.sin(angle) * distance)
#         self.size = np.random.uniform(self.gas_cloud.size * 0.2, self.gas_cloud.size * 0.6)

#         self.density = np.random.uniform(0.3, 1.0)
#         self.alpha = int(np.clip(self.density * 50, 5, 70))

#     def get_position(self):
#         return self.gas_cloud.position + self.offset

#     def draw(self, screen):
#         pygame.draw.circle(screen, (20, 15, 15, self.alpha), self.get_position(), self.size)

class DustCloud:
    def __init__(self, gas_cloud):
        self.gas_cloud = gas_cloud
        self.particles = []

        n = np.random.randint(2, 6)

        for _ in range(n):
            angle = np.random.uniform(0, 2 * np.pi)
            distance = np.random.uniform(
                0,
                gas_cloud.size * 0.7
            )

            offset = pygame.Vector2(
                np.cos(angle) * distance,
                np.sin(angle) * distance
            )

            size = np.random.uniform(
                gas_cloud.size * 0.1,
                gas_cloud.size * 0.35
            )

            density = np.random.uniform(0.2, 1.0)
            alpha = int(np.clip(density * 40, 5, 40))

            self.particles.append(
                (offset, size, alpha)
            )

    def draw(self, screen):
        for offset, size, alpha in self.particles:
            position = self.gas_cloud.position + offset

            pygame.draw.circle(
                screen,
                (20, 15, 15, alpha),
                position,
                size
            )
