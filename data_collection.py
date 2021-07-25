import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations

# Function to produce a list of distances from the centre
def distance_from_centre(balls):
    d = []
    for ball in balls:
        x = np.linalg.norm(ball.position)
        d.append(x)
    return d

# Function to produce a list of inter-ball separations
def distance_between_balls(balls):
    d = []
    pairs = combinations(range(len(balls)), 2)
    for i,j in pairs:
        dx = balls[i].position - balls[j].position
        d.append(np.linalg.norm(dx))
    return d
