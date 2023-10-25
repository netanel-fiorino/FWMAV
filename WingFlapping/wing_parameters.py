import numpy as np

######################################################################################
                #   Initial Conditions
######################################################################################
t_period = .05     #sec
freq = 1/t_period   #Hz
t_period = 1/freq   #sec (if you want to define the speed by frequency)
v_free_stream = 0*np.array([[1], [0], [0]])   #m/s
amplitude = np.radians(180)         #how many degrees the wing flaps
wing_length = 0.12      #m
num_elements = 20                    #how many elements along the wing
max_chord = 0.33*wing_length
wing_shape = 'rectangle' #'ellipse'                #ellipse or rectangle

beta = np.radians(90)   #stroke plane angle
phi = np.radians(-90)    # initial flapping angle/stroke angle
psi = 0                 #deviation angle/elevation nagle
theta = np.radians(0)               #pitch angle/rotation angle
def phi_dot_f(t):
    return amplitude/2*2*np.pi*freq*np.cos(2*np.pi*freq*t + phi) #amplitude/2*2*np.pi*freq*np.cos(2*np.pi*freq*t + np.pi/2) 
def psi_dot_f(t):
    return 0
def theta_dot_f(t):
    return theta*2*np.pi*freq*np.cos(2*np.pi*freq*t + phi) + 0
    return 0


######################################################################################
                #   Simulation Parameters
######################################################################################
t_step = t_period/100
rho = 1.2682
