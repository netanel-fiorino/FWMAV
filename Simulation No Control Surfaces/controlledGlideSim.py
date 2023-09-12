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

# initialize the visualization
VIDEO = False  # True==write video, False==don't write video
mav_view = FWMAVViewer()  # initialize the mav viewer
data_view = DataViewer()  # initialize view of data plots
if VIDEO is True:
    from graphicsViewer.video import VideoWriter
    video = VideoWriter(video_name="controlledGlide.avi",
                        bounding_box=(0, 0, 1000, 1000),
                        output_rate=SIM.ts_video)

# initialize elements of the architecture
wind = WindSimulation(SIM.ts_simulation)
fwmav = FWMAVDynamics(SIM.ts_simulation)
delta = MsgDelta()

# initialize the simulation time
sim_time = SIM.start_time
plot_time = sim_time

sim_time = SIM.start_time
tArr = []
uArr = []
vArr = []
wArr = []
northArr = []
eastArr = []
downArr = []
e0Arr = []
e1Arr = []
e2Arr = []
e3Arr = []
pArr = []
qArr = []
rArr = []
flapArr = []
# main simulation loop
print("Press Command-Q to exit...")
while sim_time < SIM.end_time:
    # -------set control surfaces-------------
    delta.rudder = -0#.0003026
    delta.throttle = 0#0.5
    delta.kappa = 0 # np.pi/2 #np.pi/4*np.sin(2*np.pi*10*sim_time)
    if ((sim_time % 10) < 5):
        delta.flap = 0
    elif ((sim_time % 10) >= 5):
        delta.flap = 1
    delta.flap = 0
    # -------physical system-------------
    current_wind = wind.update()  # get the new wind vector
    fwmav.update(delta, current_wind, sim_time)  # propagate the MAV dynamics

    # -------update viewer-------------
    if sim_time-plot_time > SIM.ts_plotting:
        mav_view.update(fwmav.true_state)  # plot body of MAV
        plot_time = sim_time
    data_view.update(fwmav.true_state,  # true states
                     fwmav.true_state,  # estimated states
                     fwmav.true_state,  # commanded states
                     delta,  # inputs to aircraft
                     SIM.ts_simulation)

    tArr.append(sim_time)
    northArr.append(fwmav._state.item(0))
    eastArr.append(fwmav._state.item(1))
    downArr.append(-fwmav._state.item(2))
    uArr.append(fwmav._state.item(3))
    vArr.append(fwmav._state.item(4))
    wArr.append(fwmav._state.item(5))
    e0Arr.append(fwmav._state.item(6))
    e1Arr.append(fwmav._state.item(7))
    e2Arr.append(fwmav._state.item(8))
    e3Arr.append(fwmav._state.item(9))
    pArr.append(fwmav._state.item(10))
    qArr.append(fwmav._state.item(11))
    rArr.append(fwmav._state.item(12))
    flapArr.append(delta.flap)

    if VIDEO is True:
        video.update(sim_time)

    # -------increment time-------------
    sim_time += SIM.ts_simulation

if VIDEO is True:
    video.close()


plt.figure()
# Loop through the data to create line segments with varying colors
for i in range(len(tArr) - 1):
    color = 'blue' if flapArr[i] == 1 else 'red'
    plt.plot(tArr[i:i+2], downArr[i:i+2], color=color)
# Add labels and title
plt.xlabel('Time')
plt.ylabel('Position')
# Add legend for colors
plt.plot([], [], color='blue', label='Flapping')
plt.plot([], [], color='red', label='Gliding')
plt.legend()
plt.grid()
plt.show()

print("FINISHED")
input("Finished")