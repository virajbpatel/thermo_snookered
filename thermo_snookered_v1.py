import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
from matplotlib import gridspec
import collisions_in_circle as collisions
import data_collection as collector
import matplotlib.animation as animation
from itertools import combinations
import random

# Maxwell-Boltzmann constant
k_B = 1.38e-23

# Essential properties of simulation
BALLS = []
EVENTS = []
TIME = 0

# Lists to store data when animate is False
KE_time = []
MOMENTUM_time = []
DISTANCE_FROM_CENTRE = []
INTERBALL_SEPARATION = []
T = []
PRESSURE = []
v = []

# Array of data storage places to manually reset simulation so that it can run multiple times
storage_units = [BALLS, EVENTS, KE_time, MOMENTUM_time, DISTANCE_FROM_CENTRE, INTERBALL_SEPARATION, T, PRESSURE, v]

pause = False

# Procedure to reset simulation and set up a new initial state
def NEW_SIM():
    for unit in storage_units:
        unit.clear()
    TIME = 0

# Hexagonal numbers are used for circle packing to ensure that the user doesn't input too many balls in the container
def hexagonal(n):
    H = n*(2*n-1)
    if n == 1:
        return H
    else:
        S = (n-1)**2
        return H + S

# Procedure to pause animation
def onClick(event):
    global pause
    pause ^= True

# Function to generate a random initial position for a ball based on a circular distrbution
def GENERATE_RANDOM_POSITION(R):
    phi = 2 * np.pi * random.random()
    r = R * np.sqrt(random.random())
    x = r * np.cos(phi)
    y = r * np.sin(phi)
    pos = [x,y]
    return pos

# Function to generate a random initial velocity for a ball based on a normal distribution centred at 0
def GENERATE_RANDOM_VELOCITY(spread):
    v_x = np.random.normal(loc = 0, scale = spread)
    v_y = np.random.normal(loc = 0, scale = spread)
    vel = [v_x, v_y]
    return vel

class Event:
    '''
    An event is a collision event and has types wall and ball
    self (obj): Event itself
    time (float): Time event occurs
    '''
    numEvents = 0

    def __init__(self, time):
        # Initialises Event object
        self.time = time
        EVENTS.append(self)
        EVENTS.sort(key = lambda x: x.time)
        Event.numEvents += 1


class WallCollision(Event):
    '''
    A subclass of Event, specified for wall collision events
    ball (obj): ball to collide with wall
    '''
    def __init__(self, ball, time):
        # Initialises WallCollision object and inherits from Event
        super().__init__(time)
        self.ball = ball

    def __repr__(self):
        # Readable representation of the object
        label = "Event" + str(EVENTS.index(self) + 1) + ": " + self.ball.__repr__() + ", t=" + str(self.time)
        return label

class BallCollision(Event):
    '''
    A subclass of Event, specified for wall collision events
    b1 (obj): first ball in collision
    b2 (obj): second ball in collision
    '''
    def __init__(self, b1, b2, time):
        # Initialises BallCollision object and inherits from Event
        super().__init__(time)
        self.b1 = b1
        self.b2 = b2

    def __repr__(self):
        # Readable representation of the object
        label = "Event" + str(EVENTS.index(self) + 1) + ": " + self.b1.__repr__() + "," + self.b2.__repr__() + ", t=" + str(self.time)
        return label


class Ball:
    '''
    A ball object has attributes mass, radius, position and velocity. It also contains a patch attribute that will be used for animations.
    The position is a protected attribute. Position and velocity are numpy arrays.
    mass (float): mass of ball
    radius (float): radius of ball
    position (float list): position of ball
    velocity (float list): velocity of ball
    '''
    def __init__(self, mass, radius, position, velocity):
        # Initialises Ball object
        self.mass = mass
        self.radius = radius
        self.position = np.array(position, dtype = float)
        self.velocity = np.array(velocity, dtype = float)
        self.kinetic_energy = 0.5 * self.mass * np.linalg.norm(self.velocity)**2
        self.momentum = self.mass * self.velocity
        self.change_in_momentum = 0
        self.patch = pl.Circle(xy = self.position, radius = self.radius, fc = 'r')
        self.next_collision_time = 99999999999

    def __repr__(self):
        # Readable representation of the ball
        label = "Ball" + str(BALLS.index(self) + 1)
        return label

    # Method to move ball to new position after time dt
    def move(self, dt, R):
        self.position = self.position + dt * self.velocity
        dp = self.momentum
        self.position, self.velocity = collisions.collide_with_wall(self, R)
        self.updateMomentum()
        dp -= self.momentum
        self.change_in_momentum = np.linalg.norm(dp)

    '''
    For calculating collision times, the following error codes are used as times:
    Ends in:        Means:
    99              No other balls in simulation
    98              Ball collision time is complex
    97              Next collision time is negative
    96              Wall collision time is complex
    ''' 

    # Method to calculate the time until the ball collides with another ball or a wall
    def time_to_wall_collision(self):
        w_t = 999999996
        A = np.linalg.norm(self.velocity)**2
        B = 2 * np.dot(self.velocity, self.position)
        C = np.linalg.norm(self.position)**2 - (10 - self.radius)**2
        if B**2 - 4*A*C >= 0:
            w_times = np.roots([A,B,C])
            w_times.sort()
            if w_times[0] > 0 and w_times[1] > 0:
                w_t = min(w_times)
            elif w_times[0] > 0 or w_times[1] > 0:
                w_t = max(w_times)
        return w_t

    def time_to_ball_collision(self, other = None):
        b_t = 99999999999
        if other != None:
            v1 = self.velocity
            v2 = other.velocity
            r1 = self.position
            r2 = other.position
            R1 = self.radius
            R2 = other.radius
            dr = r1 - r2
            dv = v1 - v2
            a = np.linalg.norm(dv)**2
            b = 2*np.dot(dr,dv)
            c = np.linalg.norm(dr)**2 - (R1 + R2)**2
            b_times = np.roots([a, b, c])
            b_times.sort()
            if (b**2 - 4*a*c) >= 0:
                if b_times[0] > 0 and b_times[1] > 0:
                    b_t = min(b_times)
                elif b_times[0] > 0 or b_times[1] > 0:
                    b_t = max(b_times)
            else:
                b_t = 99999999998
        return b_t

    def time_to_collision(self, R):
        iter_list = [x for x in BALLS if x != self]
        t = 999999999999
        for ball in iter_list:
            self.position, ball.position = collisions.overlap(self, ball, R)
            b_t = self.time_to_ball_collision(other = ball)
            if b_t < t:
                t = b_t
                obj = ball
        w_t = self.time_to_wall_collision()
        min_t = min(t, w_t)
        if min_t < 0:
            self.next_collision_time = max(t, w_t)
        else:
            self.next_collision_time = min_t
        if self.next_collision_time == w_t:
            e = WallCollision(self, time = self.next_collision_time)
        else:
            e = BallCollision(self, obj, self.next_collision_time)

    # Procedure to redefine patch of the ball when position is changed
    def set_patch(self):
        self.patch = pl.Circle(xy = self.position, radius = self.radius, fc = 'r')

    # Getter for the patch attribute, used to output ball to the screen
    def get_patch(self):
        return self.patch

    # Method to update kinetic energy of the ball
    def updateKE(self):
        self.kinetic_energy = 0.5 * self.mass * np.linalg.norm(self.velocity)**2

    # Method to update momentum of the ball
    def updateMomentum(self):
        self.momentum = self.mass * self.velocity

    # Method to check if the ball overlaps with any other ball or the wall of the container
    def overlap(self, other, R):
        isOverlapping = (np.linalg.norm(self.position - other.position) < (self.radius + other.radius)) or (R - np.linalg.norm(self.position) < self.radius)
        return isOverlapping

    # Appends ball to list of balls once there is not overlap
    def finalise(self):
        BALLS.append(self)

class Container:
    '''
    A circular container object to hold the simulation. A physical limit of the positions of the balls.
    radius (float): radius of container
    '''
    def __init__(self, radius):
        # Initialises Container object
        self.radius = radius
        self.patch = pl.Circle((0,0), self.radius, fill = False)

    # Getter for patch attribute, used to output container to the screen
    def get_patch(self):
        return self.patch

    # More readable representation of object, used for testing
    def __repr__(self):
        return "Container(radius = " + str(self.radius) + ")"

class Simulation:
    '''
    A simulation object that brings all parts of the simulation together. Contains the balls, events and the container.
    container (obj): container for simulation
    n (int): number of balls
    balls (obj list): list of ball objects, can be empty upon initialisation to store new balls
    ball_radius (float): radius of the balls, same for each ball
    '''
    def __init__(self, container, n, balls, ball_radius):
        # Initialises Simulation object
        H = hexagonal(container.radius/ball_radius)
        if n > H:
            raise Exception("Too many particles, maximum number of particles is", H)
        # Creates n balls ensuring that they dont overlap
        for i in range(n):
            canRun = True
            while canRun:
                p = GENERATE_RANDOM_POSITION(container.radius)
                v = GENERATE_RANDOM_VELOCITY(3)
                b = Ball(1e-10, ball_radius, p, v)
                for ball in BALLS:
                    if b.overlap(ball, container.radius):
                        break
                else:
                    canRun = False
            b.finalise()
        self.balls = balls
        ke = 0
        p = np.array([0, 0], dtype = float)
        # Caalculate initial kinetic energy
        for ball in self.balls:
            ke += ball.kinetic_energy
            np.add(p, ball.momentum, out = p, casting = 'unsafe')
        self.kinetic_energy = ke
        self.momentum = np.linalg.norm(p)
        self.pressure = 0
        self.temperature = (2 * self.kinetic_energy) / (3 * len(self.balls) * k_B)
        self.container = container
        self.time = 0

    # Readable representation of object used for testing
    def __repr__(self):
        return "Simulation:\nContainer radius: %s \nNumber of balls: %s"%(str(self.container.radius), str(len(self.balls)))

    def update_system(self):
        # Updates the calculated attributes of the system
        ke = 0
        p = np.array([0,0], dtype = float)
        Dp = 0
        for ball in self.balls:
            ke += ball.kinetic_energy
            Dp += ball.change_in_momentum
            np.add(p, ball.momentum, out = p, casting = 'unsafe')
            v.append(np.linalg.norm(ball.velocity))
        self.kinetic_energy = ke
        self.pressure += Dp
        self.temperature = (2 * self.kinetic_energy) / (3 * len(self.balls) * k_B)
        self.momentum = np.linalg.norm(p)

    def do_collisions(self):
        # Performs any collisions that take place at that time. Only applies to animation
        pairs = combinations(range(len(self.balls)), 2)
        for i,j in pairs:
            self.balls[i].velocity, self.balls[j].velocity, self.balls[i].position, self.balls[j].position = collisions.collide_with_ball(self.balls[i], self.balls[j])

    # Method to initialise the system graphically. Also used to output patches to the screen at each iteration
    def init_balls(self):
        self.circles = []
        self.KE_points = []
        for ball in self.balls:
            self.circles.append(self.ax1.add_patch(ball.patch))
            b, = self.ax2.plot(0.5, ball.kinetic_energy, "bo")
            self.KE_points.append(b)
        r, = self.ax2.plot(0.5, self.kinetic_energy, "ro")
        P, = self.ax2.plot(0.75, self.momentum, "go")
        self.KE_points.append(r)
        self.KE_points.append(P)
        self.shapes = self.circles + self.KE_points
        return self.shapes

    # Method to animate the simulation. This is performed at each frame in the animation.
    def animate(self, i, animate):
        if not pause:
            self.time += 0.1
            for ball in self.balls:
                ball.move(0.1, self.container.radius)
            self.do_collisions()
            for b in self.balls:
                b.updateKE()
                b.updateMomentum()
                if animate:
                    b.set_patch()
            self.update_system()
        elif pause:
            pass
        if animate:
            circles = self.init_balls()
            return circles

    # Method to run the simulation with or without the animation.
    def run(self, save = False, animate = True, num_frames = 10, multiple = False):
        if animate: 
            # Animation is done using matplotlib.animation.FuncAnimation() using the init_balls() method as an initial state
            # This is for a smoother animation and to ensure the balls do not escape or form molecules between collisions
            f = plt.figure()
            spec = gridspec.GridSpec(ncols = 2, nrows = 1, width_ratios=[3,1])
            self.ax1 = f.add_subplot(spec[0])
            self.ax2 = f.add_subplot(spec[1])
            self.ax2.set_ylabel("KE")
            self.ax1.set_xlim(-10,10)
            self.ax1.set_ylim(-10,10)
            self.ax2.set_xlim(0,1)
            self.ax2.set_ylim(0, 2*self.kinetic_energy)
            self.ax1.grid()
            self.ax2.grid()
            self.ax1.add_patch(self.container.get_patch())
            f.canvas.mpl_connect('button_press_event', onClick)
            a = animation.FuncAnimation(f, self.animate, fargs = (animate,), init_func = self.init_balls, frames = 800, interval = 2, blit = True)
            if save:
                Writer = animation.writers['ffmpeg']
                writer = Writer(fps=60, bitrate=1800)
                a.save('collision.mp4', writer=writer)
            else:
                plt.show()
        else:
            # Runs without animation in a more efficient way
            # Outputs details of each frame at each iteration
            self.TIME = 0
            for frame in range(num_frames):
                print("------------------------------------------------------------------")
                print("FRAME: ", frame)
                print("------------------------------------------------------------------")
                for b in self.balls:
                    v.append(np.linalg.norm(b.velocity))
                    b.time_to_collision(self.container.radius)
                for E in EVENTS:
                    if isinstance(E, BallCollision):
                        r1 = E.b1.position + E.time * E.b1.velocity
                        r2 = E.b2.position + E.time * E.b2.velocity
                        dr = np.linalg.norm(r2 - r1)
                        print(E, dr)
                        if dr <= (E.b1.radius + E.b2.radius):
                            print(E)
                            EVENTS[0] = E
                            break
                    else:
                        if E.time <= 5:
                            print(E)
                            EVENTS[0] = E
                            break
                t = EVENTS[0].time
                if t >= 9999999999:
                    print(EVENTS)
                    for b in self.balls:
                        print(b.position, b.velocity)
                    EVENTS.clear()
                    raise Exception("Something went wrong. There is an error with the time. Error code: ", t)
                self.TIME += t
                print("Time elapsed: ", self.TIME)
                for ball in self.balls:
                    ball.move(t, self.container.radius)
                if isinstance(EVENTS[0], BallCollision):
                    print("Before:")
                    print("r1 = ", EVENTS[0].b1.position, "v1 = ", EVENTS[0].b1.velocity)
                    print("r2 = ", EVENTS[0].b2.position, "v2 = ", EVENTS[0].b2.velocity)
                    EVENTS[0].b1.velocity, EVENTS[0].b2.velocity, EVENTS[0].b1.position, EVENTS[0].b2.position = collisions.collide_with_ball(EVENTS[0].b1, EVENTS[0].b2)
                    print(EVENTS[0], " has been done")
                    print("After:")
                    print("r1 = ", EVENTS[0].b1.position, "v1 = ", EVENTS[0].b1.velocity)
                    print("r2 = ", EVENTS[0].b2.position, "v2 = ", EVENTS[0].b2.velocity)
                for each_ball in self.balls:
                    each_ball.updateMomentum()
                    each_ball.updateKE()
                self.update_system()
                if not multiple:
                    distances_from_centre = collector.distance_from_centre(self.balls)
                    interball_separation = collector.distance_between_balls(self.balls)
                    DISTANCE_FROM_CENTRE.extend(distances_from_centre)
                    INTERBALL_SEPARATION.extend(interball_separation)
                    KE_time.append([self.TIME, self.kinetic_energy])
                    MOMENTUM_time.append([self.TIME, self.momentum])
                del EVENTS[0]
                EVENTS.clear()

