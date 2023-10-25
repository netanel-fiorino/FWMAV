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
    forceArr = []
    momentArr = []
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
    
    while time < 5*sim.t_period + 1*sim.t_step:
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
        moments = fwmav._moments
        momentArr.append([moments.item(0), moments.item(1), moments.item(2)])
        
        forces = fwmav._forces
        forceArr.append([forces.item(0), forces.item(1), forces.item(2)])
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

    # plt.figure()
    # plt.title("height")
    # plt.plot(timeArr, heightArr, label = 'height')
    # plt.plot(timeArr, forwardArr, label = 'forward')
    # plt.legend()
    # plt.grid()

    # plt.figure()
    # plt.title("pitch")
    # plt.plot(timeArr, np.degrees(pitchArr), label = 'pitch')
    # plt.legend()
    # plt.grid()
    # plt.show()

    forceArr = np.array(forceArr)
    momentArr = np.array(momentArr)
    
    plt.figure()
    plt.title("Moment")
    plt.plot(timeArr, momentArr[:, 0], label = 'Mx')
    plt.plot(timeArr, momentArr[:, 1], label = 'My')
    plt.plot(timeArr, momentArr[:, 2], label = 'Mz')
    plt.legend()
    plt.grid()

    
    plt.figure()
    plt.title("Force")
    plt.plot(timeArr, forceArr[:, 0], label = 'Fx')
    plt.plot(timeArr, forceArr[:, 1], label = 'Fy')
    plt.plot(timeArr, forceArr[:, 2], label = 'Fz')
    plt.legend()
    plt.grid()


    avgFx = sum(forceArr[:, 0])/len(forceArr[:, 0])
    avgFy = sum(forceArr[:, 1])/len(forceArr[:, 1])
    avgFz = sum(forceArr[:, 2])/len(forceArr[:, 2])
    avgMx = sum(momentArr[:, 0])/len(momentArr[:, 0])
    avgMy = sum(momentArr[:, 1])/len(momentArr[:, 1])
    avgMz = sum(momentArr[:, 2])/len(momentArr[:, 2])
    print("Average Forces and Moments over one period:")
    print(f"Fx: {avgFx}, Fy: {avgFy}, Fz: {avgFz}, Mx: {avgMx}, My: {avgMy}, Mz: {avgMz}")

    plt.show()
#run_update_once()
dual_flap()