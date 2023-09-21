import sys
sys.path.append('..')
import numpy as np
from scipy.optimize import minimize
from tools.rotations import Euler2Quaternion
from message_types.msg_delta import MsgDelta
import parameters.FWMAV_parameters as FWMAV

def compute_trim(fwmav, Va, gamma):
    # define initial state and input
    e0 = Euler2Quaternion(0., gamma, 0.)
    state0 = np.array([
        [fwmav._state.item(0)],  # pn
        [fwmav._state.item(1)],  # pe
        [fwmav._state.item(2)],  # pd
        [fwmav._state.item(3)],  # u
        [fwmav._state.item(4)],  # v
        [fwmav._state.item(5)],  # w
        [e0.item(0)],  # e0
        [e0.item(1)],  # e1
        [e0.item(2)],  # e2
        [e0.item(3)],  # e3
        [fwmav._state.item(10)],  # p
        [fwmav._state.item(11)],  # q
        [fwmav._state.item(12)]   # r
    ])

    delta0 = MsgDelta()
    
    x0 = np.concatenate((state0, delta0.to_array()), axis=0)
    # define equality constraints
    # bnds = ((None, None), (None, None), (None, None), (None, None),
    #         (None, None), (None, None), (None, None), (None, None), 
    #         (None, None), (None, None), (None, None), (None, None), (None, None),
    #         (-1.0, 1.0), (-1.0, 1.0), (-1.0, 1.0), (0.0, 1.0))
    cons = ({'type': 'eq',
             'fun': lambda x: np.array([
                                x[3]**2 + x[4]**2 + x[5]**2 - Va**2,  # magnitude of velocity vector is Va
                                x[4],  # v=0, force side velocity to be zero
                                x[6]**2 + x[7]**2 + x[8]**2 + x[9]**2 - 1.,  # force quaternion to be unit length
                                x[7],   # e1=0  - forcing e1=e3=0 ensures zero roll and zero yaw in trim
                                x[9],   # e3=0
                                x[10],  # p=0  - angular rates should all be zero
                                x[11],  # q=0
                                x[12],  # r=0
                                ]),
             'jac': lambda x: np.array([
                                [0., 0., 0., 2*x[3], 2*x[4], 2*x[5], 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0., 0., 2*x[6], 2*x[7], 2*x[8], 2*x[9], 0., 0., 0., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0.],
                                ])
             })
    # solve the minimization problem to find the trim states and inputs
    res = minimize(trim_objective_fun, np.squeeze(x0), method='SLSQP', args=(fwmav, Va, gamma), #bounds=bnds,
                   constraints=cons, options={'ftol': 1e-10, 'disp': True})
    # extract trim state and input and return
    trim_state = np.array([res.x[0:13]]).T
    trim_input = MsgDelta(rudder=res.x.item(13),
                          throttle=res.x.item(14),
                          kappa=res.x.item(15),
                          flap=res.x.item(16),
                          tail_angle=res.x.item(17))
    trim_input.print()
    print('trim_state=', trim_state.T)

    return trim_state, trim_input

def trim_objective_fun(x, fwmav, Va, gamma):
    # objective function to be minimized

    #slide 19 chap5
    state = x[0:13]
    delta = MsgDelta(
        rudder=x.item(13),
        throttle=x.item(14),
        kappa=x.item(15),
        flap=x.item(16),
        tail_angle=x.item(17)
    )
    delta = MsgDelta(
        rudder=x.item(13),
        throttle=x.item(14),
        kappa=x.item(15),
        flap=0,
        tail_angle=x.item(17)
    )
    desired_trim_state_dot = np.array([
        [0., 0., -Va*np.sin(gamma), 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]
    ]).T
    fwmav._state = state
    fwmav._update_velocity_data(delta)
    forces_moments = fwmav._forces_moments(delta, time = 0)
    f = fwmav._derivatives(state, forces_moments)
    tmp = desired_trim_state_dot - f
    if (np.abs(delta.tail_angle) < 10**10):
        print(delta.tail_angle)
        print(tmp)
    J = np.linalg.norm(tmp[2:13])**2.0
    return J