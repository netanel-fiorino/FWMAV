import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('..')

from wingDynamics import WingDynamics
from bodyDynamics import FWMAVDynamics
from windSimulation import WindSimulation
import wing_parameters as sim
from message_types.msg_delta import MsgDelta
from tools.rotations import Quaternion2Rotation, Quaternion2Euler, Euler2Rotation
from graphicsViewer.viewer import FWMAVViewer
from graphicsViewer.dataViewer import DataViewer
import parameters.simulation_parameters as SIM

plt.close('all')

def run_update_once():
    wing_left = WingDynamics(side='left', wing_shape='ellipse')         #rectangle ellipse
    wing_right = WingDynamics(side='right', wing_shape='ellipse')
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
    phiArr_left = []
    phiArr_right = []
    thetaArr_left = []
    thetaArr_right = []
    heightArr = []
    forwardArr = []
    pitchArr = []
    time = 0
    VIDEO = False  # True==write video, False==don't write video
    fwmav_view = FWMAVViewer()  # initialize the mav viewer
    data_view = DataViewer()  # initialize view of data plots
    if VIDEO is True:
        from graphicsViewer.video import VideoWriter
        video = VideoWriter(video_name="controlledGlide.avi",
                            bounding_box=(0, 0, 1000, 1000),
                            output_rate=SIM.ts_video)
    wing_left = WingDynamics(side='left', wing_shape='rectangle')               #rectangle ellipse
    wing_right = WingDynamics(side='right', wing_shape='rectangle')
    wind = WindSimulation(sim.t_step)
    fwmav = FWMAVDynamics(sim.t_step, wing_left, wing_right)
    delta = MsgDelta()
    
    # omega = np.pi**2*freq*np.cos(2*np.pi*freq*time + np.pi/2)
    # wing_left.updateOmega(omega)
    # wing_right.updateOmega(omega)
    
    while time < 2*sim.t_period + 1*sim.t_step:
        # tmp_left = wing_left.force_calc()
        # phiArr_left.append(wing_left._state.item(7))
        # thetaArr_left.append(wing_left._state.item(9))
        # wing_left.update(time)
        
        # tmp_right = wing_right.force_calc()
        # phiArr_right.append(wing_right._state.item(7))
        # thetaArr_right.append(wing_right._state.item(9))
        # wing_right.update(time)

        # timeArr.append(time)

        # forceArr_left.append(tmp_left[:3])
        # momentArr_left.append(tmp_left[3:6])

        # forceArr_right.append(tmp_right[:3])
        # momentArr_right.append(tmp_right[3:6])
        delta.rudder = -0
        delta.throttle = 0
        delta.kappa = 0
        delta.tail_angle = np.radians(0)
        delta.flap = 0

        
        current_wind = wind.update()
        fwmav.update(delta, current_wind, time)
        heightArr.append(-fwmav._state.item(2))
        forwardArr.append(fwmav._state.item(0))
        pitchArr.append(Quaternion2Euler(fwmav._state[6:10])[1])
        time = time + sim.t_step
        timeArr.append(time)
        fwmav_view.update(fwmav.true_state)
        data_view.update(fwmav.true_state,  # true states
                     fwmav.true_state,  # estimated states
                     fwmav.true_state,  # commanded states
                     delta,  # inputs to aircraft
                     sim.t_step)
        if VIDEO is True:
            video.update(time)
    # forceArr_left = np.array(forceArr_left) 
    # forceArr_right = np.array(forceArr_right)
    # momentArr_left = np.array(momentArr_left) 
    # momentArr_right = np.array(momentArr_right)

    # fig, ax = plt.subplots(1, 2, figsize=(12, 6))
    # fig.suptitle(f"Forces and Moments from Dual Flapping, wing shape: {wing_left.wing_shape}")
    # ##Force Plot
    # #left wing
    # # ax[0].plot(timeArr, forceArr_left[:, 0], label = 'fx from left')
    # # ax[0].plot(timeArr, forceArr_left[:, 1], label = 'fy from left')
    # # ax[0].plot(timeArr, forceArr_left[:, 2], label = 'fz from left')
    # # #right wing
    # # ax[0].plot(timeArr, forceArr_right[:, 0], label = 'fx from right')
    # # ax[0].plot(timeArr, forceArr_right[:, 1], label = 'fy from right')
    # # ax[0].plot(timeArr, forceArr_right[:, 2], label = 'fz from right')
    # #total
    # ax[0].plot(timeArr, forceArr_left[:, 0] + forceArr_right[:, 0], '-', linewidth=3, alpha = 0.75, label = 'fx')
    # ax[0].plot(timeArr, forceArr_left[:, 1] + forceArr_right[:, 1], '--', linewidth=3, alpha = 0.75, label = 'fy')
    # ax[0].plot(timeArr, forceArr_left[:, 2] + forceArr_right[:, 2], '-', linewidth=3, alpha = 0.75, label = 'fz')
    
    # ax[0].set_title("Forces from dual flapping")
    # ax[0].set_xlabel("Time [s]")
    # ax[0].set_ylabel("Force [N]")
    # ax[0].grid()
    # ax[0].legend()

    # ##Moment Plot
    # # #left
    # # ax[1].plot(timeArr, momentArr_left[:, 0], label = 'Mx left')
    # # ax[1].plot(timeArr, momentArr_left[:, 1], label = 'My left')
    # # ax[1].plot(timeArr, momentArr_left[:, 2], label = 'Mz left')
    # # # # #right
    # # ax[1].plot(timeArr, momentArr_right[:, 0], label = 'Mx right')
    # # ax[1].plot(timeArr, momentArr_right[:, 1], label = 'My right')
    # # ax[1].plot(timeArr, momentArr_right[:, 2], label = 'Mz right')
    # #total
    # ax[1].plot(timeArr, momentArr_left[:, 0] + momentArr_right[:, 0], '--', linewidth=3, alpha = 0.75, label = 'Mx')
    # ax[1].plot(timeArr, momentArr_left[:, 1] + momentArr_right[:, 1], '-', linewidth=3, alpha = 0.75, label = 'My')
    # ax[1].plot(timeArr, momentArr_left[:, 2] + momentArr_right[:, 2], '--', linewidth=3, alpha = 0.75, label = 'Mz')
    
    # ax[1].set_title("Moments from dual flapping")
    # ax[1].set_xlabel("Time [s]")
    # ax[1].set_ylabel("Moment [Nm]")
    # ax[1].grid()
    # ax[1].legend()
    
    # plt.tight_layout()
    # #plt.show()

    # plt.figure()
    # plt.title("Angles")
    # plt.plot(timeArr, np.degrees(phiArr_left), label = 'phi left')
    # plt.plot(timeArr, np.degrees(thetaArr_left), label = 'theta left')
    # plt.plot(timeArr, np.degrees(phiArr_right), label = 'phi right')
    # plt.plot(timeArr, np.degrees(thetaArr_right), label = 'theta right')
    # plt.legend()
    # plt.grid()

    plt.figure()
    plt.title("height")
    plt.plot(timeArr, heightArr, label = 'height')
    plt.plot(timeArr, forwardArr, label = 'forward')
    plt.legend()
    plt.grid()

    plt.figure()
    plt.title("pitch")
    plt.plot(timeArr, np.degrees(pitchArr), label = 'pitch')
    plt.legend()
    plt.grid()
    plt.show()
    


    # avgFx = sum(forceArr_left[:, 0] + forceArr_right[:, 0])/len(forceArr_left[:, 0])
    # avgFy = sum(forceArr_left[:, 1] + forceArr_right[:, 1])/len(forceArr_left[:, 1])
    # avgFz = sum(forceArr_left[:, 2] + forceArr_right[:, 2])/len(forceArr_left[:, 2])
    # avgMx = sum(momentArr_left[:, 0] + momentArr_right[:, 0])/len(momentArr_left[:, 0])
    # avgMy = sum(momentArr_left[:, 1] + momentArr_right[:, 1])/len(momentArr_left[:, 1])
    # avgMz = sum(momentArr_left[:, 2] + momentArr_right[:, 2])/len(momentArr_left[:, 2])
    # print("Average Forces and Moments over one period:")
    # print(f"Fx: {avgFx}, Fy: {avgFy}, Fz: {avgFz}, Mx: {avgMx}, My: {avgMy}, Mz: {avgMz}")

    # plt.show()
#run_update_once()
dual_flap()