import numpy as np
from itertools import combinations

# Generates a matrix to perform an anticlockwise rotation
def generate_anti_clockwise_matrix(angle):
    M = np.array([
        [np.cos(angle), np.sin(angle)],
        [-np.sin(angle), np.cos(angle)]
    ])
    return M

# Generates a matrix to perform a clockwise rotation
def generate_clockwise_matrix(angle):
    M = np.array([
        [np.cos(angle), -np.sin(angle)],
        [np.sin(angle), np.cos(angle)]
    ])
    return M

# Performs a collision against a wall using the velocity and the rotation matrix
def wall_collision(v, M):
    v[1] *= -1
    vel = np.matmul(M, v)
    return vel

# Checks if two balls are overlapping or if 
def overlap(b1, b2, R):
    r1, r2 = b1.position, b2.position
    d = np.linalg.norm(r1)
    R1, R2 = b1.radius, b2.radius
    R = R1 + R2
    dr = r2 - r1
    D = np.linalg.norm(dr)
    if d + R1 > R:
        unit_r = (1 / d) * r1
        new_r1 = ((R-0.01) - R1) * unit_r
        new_r2 = r2
    else:
        new_r1 = r1
        new_r2 = r2
    if D < R:
        new_r1 = r1 - (3 * (R1 + R2 - D) / (5*D)) * dr
        new_r2 = r2 + (3*(R1 + R2 - D)/(5*D))*dr
    else:
        new_r1 = r1
        new_r2 = r2
    return new_r1, new_r2

# Performs collision between two balls
def collide_with_ball(b1, b2):
    R1, R2 = b1.radius, b2.radius
    R = R1**2 + R2**2
    r1, r2 = b1.position, b2.position
    dr = r2 - r1
    D = np.linalg.norm(dr)
    v1, v2 = b1.velocity, b2.velocity
    if D <= (R1 + R2):
        d = D**2
        u1 = v1 - 2*R2**2 / R * np.dot(v1-v2, r1-r2) / d * (r1 - r2)
        u2 = v2 - 2*R1**2 / R * np.dot(v2-v1, r2-r1) / d * (r2 - r1)
        print(b1, " and ", b2, " just collided")
        spacing = np.linalg.norm(r1 - r2)
        new_r1 = r1 - (3*(R1 + R2 - D)/(5*D))*dr
        new_r2 = r2 + (3*(R1 + R2 - D)/(5*D))*dr
    else:
        u1, u2 = v1, v2
        new_r1, new_r2 = r1, r2
    return u1, u2, new_r1, new_r2

# Performs collision against wall
def collide_with_wall(b,R):
    radius = b.radius
    pos = b.position
    vel = b.velocity
    d = np.linalg.norm(pos)
    if d + radius >= R:
        unit_r = (1/d) * pos
        r = R * unit_r
        m = -(r[0] / r[1])
        if r[0] < 0:
            theta = np.arctan(m) + np.pi
        else:
            theta = np.arctan(m)
        M = generate_anti_clockwise_matrix(theta)
        inverse_M = generate_clockwise_matrix(theta)
        v = np.matmul(M, vel)
        new_vel = wall_collision(v, inverse_M)
        new_pos = ((R-0.01) - radius) * unit_r
        print(b, " just collided with the wall")
    else:
        new_pos = pos
        new_vel = vel
    return new_pos, new_vel
        
