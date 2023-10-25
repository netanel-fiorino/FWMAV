import sys
sys.path.append('..')
import numpy as np

######################################################################################
                #   sample times, etc
######################################################################################
ts_simulation = 0.0005  # smallest time step for simulation
start_time = 0.  # start time for simulation
end_time = .5  # end time for simulation

ts_plotting = 0.001  # refresh rate for plots

#ts_video = 60  # write rate for video

ts_control = ts_simulation  # sample rate for the controller