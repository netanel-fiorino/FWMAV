import numpy as np
######################################################################################
                #   Initial Conditions
######################################################################################
t_period = .05     #sec
freq = 1/t_period   #Hz
t_period = 1/freq   #sec (if you want to define the speed by frequency)
v_free_stream = 1*np.array([[1], [0], [0]])   #m/s
amplitude = np.radians(180)         #how many degrees the wing flaps
wing_length = 0.10      #m
wing_chord = 0.05       #m
beta = np.radians(50)   #stroke plane angle
phi = np.radians(-90)    #flapping angle
psi = 0                 #deviation angle
theta = 0
def phi_dot_f(t):
    return amplitude/2*2*np.pi*freq*np.cos(2*np.pi*freq*t + phi) #amplitude/2*2*np.pi*freq*np.cos(2*np.pi*freq*t + np.pi/2) 
def psi_dot_f(t):
    return 0
def theta_dot_f(t):
    return 0

######################################################################################
                #   Simulation Parameters
######################################################################################
num_elements = 5                    #how many elements along the wing
t_step = t_period/100

rho = 1.2682