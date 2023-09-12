import sys
sys.path.append('..')
import numpy as np
import matplotlib.pyplot as plt
import parameters.simulation_parameters as SIM


from graphicsViewer.viewer import FWMAVViewer
from graphicsViewer.dataViewer import DataViewer
from dynamics import FWMAVDynamics
from message_types.msg_delta import MsgDelta
import parameters.FWMAV_parameters as FWMAV
from tools.rotations import Quaternion2Rotation, Quaternion2Euler, Euler2Rotation

# initialize the visualization
VIDEO = False  # True==write video, False==don't write video
mav_view = FWMAVViewer()  # initialize the mav viewer
data_view = DataViewer()  # initialize view of data plots
if VIDEO is True:
    from graphicsViewer.video import VideoWriter
    video = VideoWriter(video_name="gliding.avi",
                        bounding_box=(0, 0, 1000, 1000),
                        output_rate=SIM.ts_video)
print("HERE")
# initialize elements of the architecture
fwmav = FWMAVDynamics(SIM.ts_simulation)
delta = MsgDelta()

# initialize the simulation time
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
# main simulation loop
print("Press Command-Q to exit...")
while sim_time < SIM.end_time:
    # -------vary forces and moments to check dynamics-------------
    f = Euler2Rotation(fwmav.true_state.theta, fwmav.true_state.phi, fwmav.true_state.psi).T @ np.array([[0], [0], [FWMAV.mass*FWMAV.gravity]])
    fx = f.item(0)
    fy = f.item(1)  # 10
    fz = f.item(2)  # 10
    Mx = 0.  # 0.1
    My = 0.  # 0.1
    Mz = 0.  # 0.1
    forces_moments = np.array([[fx, fy, fz, Mx, My, Mz]]).T

    # -------physical system-------------
    fwmav.updateForce(forces_moments)  # propagate the MAV dynamics
  
    # -------update viewer-------------
    mav_view.update(fwmav.true_state)  # plot body of MAV
    data_view.update(fwmav.true_state,  # true states
                     fwmav.true_state,  # estimated states
                     fwmav.true_state,  # commanded states
                     delta,  # inputs to the aircraft
                     SIM.ts_simulation)
    
    tArr.append(sim_time)
    northArr.append(fwmav._state.item(0))
    eastArr.append(fwmav._state.item(1))
    downArr.append(fwmav._state.item(2))
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


    if VIDEO is True:
        video.update(sim_time)

    # -------increment time-------------
    sim_time += SIM.ts_simulation

if VIDEO is True:
    video.close()
