import numpy as np
import matplotlib.pyplot as plt
import thermo_snookered_v1 as backend

m = 1e-10
k_B = 1.38e-23
N = 20

def run_animation(N, b_radius):
    '''
    Runs simulation with animation but without data collection.
    '''
    container = backend.Container(10)
    sim = backend.Simulation(container, N, backend.BALLS, b_radius)
    sim.run(animate = True)

def find_amp(hist, peak):
    # Find amplitude of the Maxwell-Boltzmann distribution fit
    M = max(hist)
    return M/peak

def maxwell_boltzmann(v, amp, T):
    # Formula for Maxwell-Boltzmann distribution
    ke = 0.5 * m * v**2
    return amp * v * np.exp(-ke/(k_B * T))

def vanDerWaals(T, N, container_radius, ball_radius, a):
    # Formula for van Der Waals force
    b = 2*ball_radius**2*(3*np.sqrt(3) - np.pi) + 2*np.pi
    V = np.pi * container_radius**2
    return ((N*k_B) / (V - N*b))*T - a*(N/V)**2

def distances_from_centre_hist(save = False):
    '''
    Plots a histogram of the distribution of distances from the centre for each ball over a period of time (100 frames).
    '''
    container = backend.Container(10)
    sim = backend.Simulation(container, 20, backend.BALLS, 1)
    sim.run(animate = False, num_frames = 100, multiple = False)
    plt.hist(backend.DISTANCE_FROM_CENTRE)
    plt.xlabel("Distance from centre")
    if save:
        plt.savefig("T-9-1.png")

def inter_ball_separation(save = False):
    '''
    Plots a histogram of the distribution of the inter-ball separation over a period of time (100  frames),
    '''
    container = backend.Container(10)
    sim = backend.Simulation(container, 20, backend.BALLS, 1)
    sim.run(animate = False, num_frames = 100, multiple = False)
    plt.hist(backend.INTERBALL_SEPARATION)
    plt.xlabel("Inter-ball separation")
    if save:
        plt.savefig("T-9-2.png")

def KE_time_graph(save = False):
    '''
    Plots a graph of kinetic energy of the system against time.
    '''
    container = backend.Container(10)
    sim = backend.Simulation(container, 20, backend.BALLS, 1)
    sim.run(animate = False, num_frames = 100, multiple = False)
    data = np.transpose(backend.KE_time)
    plt.plot(data[0], data[1])
    plt.xlabel("Time (s)")
    plt.ylabel("KE")
    if save:
        plt.savefig("T-11-1.png")

def momentum_time_graph(save = False):
    '''
    Plots a graph of momentum of the system against time.
    '''
    container = backend.Container(10)
    sim = backend.Simulation(container, 20, backend.BALLS, 1)
    sim.run(animate = False, num_frames = 200, multiple = False)
    data = np.transpose(backend.MOMENTUM_time)
    plt.plot(data[0], data[1])
    plt.xlabel("Time (s)")
    plt.ylabel("Momentum")
    if save:
        plt.savefig("T-11-2.png")

def pressure_temperature_graph(num_data_points = 10, graph = True, multiple = False, r = 1, save = False):
    '''
    Plots graph of pressure against temperature over a period of 200 frames. Can be run for different radii and for
    comparison to the van Der Waals force.
    '''
    T = []
    P = []
    container = backend.Container(10)
    V = np.pi * container.radius**2
    # Must take readings of different simulations to find a gradient of the relationship between P and T
    for i in range(num_data_points):
        sim = backend.Simulation(container, 20, backend.BALLS, r)
        sim.run(animate = False, num_frames = 200, multiple = True)
        T.append(sim.temperature)
        P.append(sim.pressure / (2 * np.pi * container.radius * sim.TIME))
        backend.NEW_SIM()
    if graph:
        plt.plot(T, P, "bo")
        m, c = np.polyfit(T, P, 1)
        P_fit = [m*t + c for t in T]
        plt.plot(T, P_fit, "r")
        g = 20*k_B/V
        P_model = [g * t for t in T]
        plt.plot(T, P_model, "g")
        plt.xlabel("Temperature (K)")
        plt.ylabel("Pressure (Pa)")
        if save:
            plt.savefig("T-11-3.png")
    elif multiple:
        # Used to see how the equation of state changes with the radius
        m, c = np.polyfit(T, P, 1)
        return m
    else:
        # Used for the vdW equation of state
        plt.plot(T, P, "bo", label = "Data")
        plt.xlabel("Temperature (K)")
        plt.ylabel("Pressure (Pa)")
        gradient, c = np.polyfit(T, P, 1)
        V = np.pi * container.radius**2
        '''
        b = (V/N) - (k_B/gradient)
        print(b)
        '''
        a = -1 * c * (V / N)**2
        return a, T

def radius_variation(save = False):
    '''
    Plots a graph of how the equation of state varies for different ball radii.
    '''
    radii = [0.25, 0.5, 0.75, 1, 1.25]
    m = []
    # Loops over different radii
    for radius in radii:
        gradient = pressure_temperature_graph(graph = False, r = radius)
        m.append(gradient)
    plt.plot(radii, m)
    plt.xlabel("Radius")
    plt.ylabel(r'$\frac{Nk_B}{V}$')
    if save:
        plt.savefig("T-12.png")

def ball_speed_hist(save = False):
    '''
    Plots histogram of distribution of ball velocities over a period of 200 frames and fits the shape to that of a 
    Maxwell-Boltzmann distribution.
    '''
    container = backend.Container(10)
    sim = backend.Simulation(container, 20, backend.BALLS, 1)
    sim.run(animate = False, num_frames = 200, multiple = False)
    T = sim.temperature
    peak = np.sqrt((k_B * T)/m) * np.exp(-0.5)
    hist, bin_edges = np.histogram(backend.v)
    amp = find_amp(hist, peak)
    x = np.arange(0, max(bin_edges), 0.2)
    v = [maxwell_boltzmann(i, amp, T) for i in x]
    plt.hist(backend.v)
    plt.plot(x, v)
    plt.xlabel("Speed")
    plt.ylabel("Number")
    if save:
        plt.savefig("T-13.png")

def equationOfState(save = False):
    '''
    Plots P-T relationship and compares to the vdW equation of state equation.
    '''
    a, T = pressure_temperature_graph(graph = False, multiple = False, r = 1) 
    P_vdw = [vanDerWaals(t, N, 10, 1, a) for t in T]
    plt.plot(T, P_vdw, label = "van Der Waals model")
    plt.legend()
    if save:
        plt.savefig("T-14.png")



'''
ANIMATION:
To run the simulation with the animation, uncomment the following line. Change the values in the brackets to adjust the ball radius
and the number of particles.
'''
#run_animation(20, 1)

'''
TASK 9:
To plot a distribution of the distance of the balls with the centre of the container, uncomment the following line.
'''
#distances_from_centre_hist()

'''
TASK 9:
To plot a distribution of the inter-ball separation, uncomment the following line.
'''
#inter_ball_separation()

'''
TASK 11:
To plot a graph of the kinetic energy against time, uncomment the following line.
'''
#KE_time_graph()

'''
TASK 11:
To plot a graph of the momentum against time, uncomment the following line.
'''
#momentum_time_graph()

'''
TASK 11:
To plot a graph of the pressure against temperature, uncomment the following line.
'''
#pressure_temperature_graph()

'''
TASK 12:
To plot a graph of how the equation of state varies with ball radius, uncomment the following line.
'''
#radius_variation()

'''
TASK 13:
To plot a graph of the distribution of the velocities, uncomment the following line.
'''
#ball_speed_hist()

'''
TASK 14:
To plot a graph comparing the pressure-temperature relationship to that of the van Der Waal's law, uncomment the following line.
'''
#equationOfState()

plt.show()