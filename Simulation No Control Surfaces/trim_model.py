import sys
sys.path.append('..')
import numpy as np
import matplotlib.pyplot as plt
import parameters.simulation_parameters as SIM


from graphicsViewer.viewer import FWMAVViewer
from graphicsViewer.dataViewer import DataViewer
from dynamics import FWMAVDynamics
from windSimulation import WindSimulation
from message_types.msg_delta import MsgDelta
import parameters.FWMAV_parameters as FWMAV
from tools.rotations import Quaternion2Rotation, Quaternion2Euler, Euler2Rotation
from trim import compute_trim

# initialize the visualization
VIDEO = False  # True==write video, False==don't write video
fwmav_view = FWMAVViewer()  # initialize the mav viewer
data_view = DataViewer()  # initialize view of data plots
if VIDEO is True:
    from graphicsViewer.video import VideoWriter
    video = VideoWriter(video_name="controlledGlide.avi",
                        bounding_box=(0, 0, 1000, 1000),
                        output_rate=SIM.ts_video)

# initialize elements of the architecture
wind = WindSimulation(SIM.ts_simulation)
fwmav = FWMAVDynamics(SIM.ts_simulation)
# use compute_trim function to compute trim state and trim input
Va = 5.
gamma = 0.*np.pi/180.
trim_state, trim_input = compute_trim(fwmav, Va, gamma)
fwmav._state = trim_state  # set the initial state of the mav to the trim state
delta = trim_input  # set input to constant constant trim input

# initialize the simulation time
sim_time = SIM.start_time
plot_time = sim_time

sim_time = SIM.start_time

#main simulation loop
print("Press Command-Q to exit...")
while sim_time < SIM.end_time:

    # -------physical system-------------
    #current_wind = wind.update()  # get the new wind vector
    current_wind = np.zeros((6, 1))
    # this input excites the phugoid mode by adding an impulse at t=5.0
    # delta.elevator += input_signal.impulse(sim_time)
    # delta.rudder += input_signal.doublet(sim_time)
    fwmav.update(delta, current_wind, sim_time)

    # -------update viewer-------------
    if sim_time-plot_time > SIM.ts_plotting:
        fwmav_view.update(fwmav.true_state)  # plot body of MAV
        plot_time = sim_time
    data_view.update(fwmav.true_state,  # true states
                     fwmav.true_state,  # estimated states
                     fwmav.true_state,  # commanded states
                     delta,  # inputs to aircraft
                     SIM.ts_simulation)
    if VIDEO is True:
        video.update(sim_time)

    # -------increment time-------------
    sim_time += SIM.ts_simulation

if VIDEO is True:
    video.close()