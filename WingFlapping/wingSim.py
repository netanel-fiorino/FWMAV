import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('..')

from wingDynamics import WingDynamics
import wing_parameters as sim
plt.close('all')

freq = sim.freq #1/t_period   #Hz
amplitude = sim.amplitude #np.radians(180)
def run__update_once():
    wing_left = WingDynamics('left')
    wing_right = WingDynamics('right')
    time = 0
    #omega = amplitude/2*2*np.pi*freq*np.cos(2*np.pi*freq*time + np.pi/2)
    # wing_left.updateOmega(omega(time))
    # wing_right.updateOmega(omega(time))

    tmp_left = wing_left.force_calc()
    wing_left.update(time)
    tmp_right = wing_right.force_calc()
    wing_right.update(time)

def dual_flap():
    timeArr = []
    forceArr_left = []
    forceArr_right = []
    momentArr_left = []
    momentArr_right = []
    angleArr_left = []
    angleArr_right = []
    time = 0

    wing_left = WingDynamics('left')
    wing_right = WingDynamics('right')
    
    # omega = np.pi**2*freq*np.cos(2*np.pi*freq*time + np.pi/2)
    # wing_left.updateOmega(omega)
    # wing_right.updateOmega(omega)
    
    while time < sim.t_period/2:
        tmp_left = wing_left.force_calc()
        angleArr_left.append(wing_left._state.item(7))
        wing_left.update(time)
        tmp_right = wing_right.force_calc()
        wing_right.update(time)

        timeArr.append(time)

        forceArr_left.append(tmp_left[:3])
        momentArr_left.append(tmp_left[3:6])

        forceArr_right.append(tmp_right[:3])
        momentArr_right.append(tmp_right[3:6])

        
        time = time + sim.t_step
        # omega = np.pi**2*freq*np.cos(2*np.pi*freq*time + np.pi/2) #2*np.pi*freq *np.cos(2*np.pi*freq*time)
        # wing_left.updateOmega(omega)
        # wing_right.updateOmega(omega)

    while time < sim.t_period + sim.t_step:
        tmp_left = wing_left.force_calc()
        angleArr_left.append(wing_left._state.item(7))
        wing_left.update(time)
        tmp_right = wing_right.force_calc()
        wing_right.update(time)

        timeArr.append(time)

        forceArr_left.append(tmp_left[:3])
        momentArr_left.append(tmp_left[3:6])

        forceArr_right.append(tmp_right[:3])
        momentArr_right.append(tmp_right[3:6])

        
        
        time = time + sim.t_step
        # omega = amplitude/2*2*np.pi*freq*np.cos(2*np.pi*freq*time + np.pi/2)
        # wing_left.updateOmega(omega)
        # wing_right.updateOmega(omega)
        
        
    forceArr_left = np.array(forceArr_left) 
    forceArr_right = np.array(forceArr_right)
    momentArr_left = np.array(momentArr_left) 
    momentArr_right = np.array(momentArr_right)

    fig, ax = plt.subplots(1, 2, figsize=(12, 6))
    fig.suptitle("Forces and Moments from Dual Flapping")
    ##Force Plot
    #left wing
    # ax[0].plot(timeArr, forceArr_left[:, 0], label = 'fx from left')
    # ax[0].plot(timeArr, forceArr_left[:, 1], label = 'fy from left')
    # ax[0].plot(timeArr, forceArr_left[:, 2], label = 'fz from left')
    # #right wing
    # ax[0].plot(timeArr, forceArr_right[:, 0], label = 'fx from right')
    # ax[0].plot(timeArr, forceArr_right[:, 1], label = 'fy from right')
    # ax[0].plot(timeArr, forceArr_right[:, 2], label = 'fz from right')
    #total
    ax[0].plot(timeArr, forceArr_left[:, 0] + forceArr_right[:, 0], '-', linewidth=3, alpha = 0.75, label = 'fx')
    ax[0].plot(timeArr, forceArr_left[:, 1] + forceArr_right[:, 1], '--', linewidth=3, alpha = 0.75, label = 'fy')
    ax[0].plot(timeArr, forceArr_left[:, 2] + forceArr_right[:, 2], ':', linewidth=3, alpha = 0.75, label = 'fz')
    
    ax[0].set_title("Forces from dual flapping")
    ax[0].set_xlabel("Time [s]")
    ax[0].set_ylabel("Force [N]")
    ax[0].grid()
    ax[0].legend()

    ##Moment Plot
    # #left
    # ax[1].plot(timeArr, momentArr_left[:, 0], label = 'Mx left')
    # ax[1].plot(timeArr, momentArr_left[:, 1], label = 'My left')
    # ax[1].plot(timeArr, momentArr_left[:, 2], label = 'Mz left')
    # #right
    # ax[1].plot(timeArr, momentArr_right[:, 0], label = 'Mx right')
    # ax[1].plot(timeArr, momentArr_right[:, 1], label = 'My right')
    # ax[1].plot(timeArr, momentArr_right[:, 2], label = 'Mz right')
    #total
    ax[1].plot(timeArr, momentArr_left[:, 0] + momentArr_right[:, 0], '--', linewidth=3, alpha = 0.75, label = 'Mx')
    ax[1].plot(timeArr, momentArr_left[:, 1] + momentArr_right[:, 1], '-', linewidth=3, alpha = 0.75, label = 'My')
    ax[1].plot(timeArr, momentArr_left[:, 2] + momentArr_right[:, 2], '-', linewidth=3, alpha = 0.75, label = 'Mz')
    
    ax[1].set_title("Moments from dual flapping")
    ax[1].set_xlabel("Time [s]")
    ax[1].set_ylabel("Moment [Nm]")
    ax[1].grid()
    ax[1].legend()
    
    plt.tight_layout()
    plt.show()

    # plt.figure()
    # plt.plot(timeArr, np.degrees(angleArr_left), label = 'phi')
    # plt.legend()
    # plt.grid()
    # plt.show()
dual_flap()