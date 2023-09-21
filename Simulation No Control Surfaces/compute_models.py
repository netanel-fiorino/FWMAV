import sys
sys.path.append('..')
import numpy as np
from scipy.optimize import minimize
from tools.rotations import Euler2Quaternion, Quaternion2Euler
import parameters.FWMAV_parameters as FWMAV
from parameters.simulation_parameters import ts_simulation as Ts
from message_types.msg_delta import MsgDelta


def compute_model(FWMAV, trim_state, trim_input):
    A_lon, B_lon, A_lat, B_lat = compute_ss_model(FWMAV, trim_state, trim_input)
    Va_trim, alpha_trim, theta_trim, a_phi1, a_phi2, a_theta1, a_theta2, a_theta3, \
    a_V1, a_V2, a_V3 = compute_tf_model(FWMAV, trim_state, trim_input)

    # write transfer function gains to file
    file = open('model_coef.py', 'w')
    file.write('import numpy as np\n')
    file.write('x_trim = np.array([[%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f]]).T\n' %
               (trim_state.item(0), trim_state.item(1), trim_state.item(2), trim_state.item(3),
                trim_state.item(4), trim_state.item(5), trim_state.item(6), trim_state.item(7),
                trim_state.item(8), trim_state.item(9), trim_state.item(10), trim_state.item(11),
                trim_state.item(12)))
    file.write('u_trim = np.array([[%f, %f, %f, %f]]).T\n' %
               (trim_input.elevator, trim_input.aileron, trim_input.rudder, trim_input.throttle))
    file.write('Va_trim = %f\n' % Va_trim)
    file.write('alpha_trim = %f\n' % alpha_trim)
    file.write('theta_trim = %f\n' % theta_trim)
    file.write('a_phi1 = %f\n' % a_phi1)
    file.write('a_phi2 = %f\n' % a_phi2)
    file.write('a_theta1 = %f\n' % a_theta1)
    file.write('a_theta2 = %f\n' % a_theta2)
    file.write('a_theta3 = %f\n' % a_theta3)
    file.write('a_V1 = %f\n' % a_V1)
    file.write('a_V2 = %f\n' % a_V2)
    file.write('a_V3 = %f\n' % a_V3)
    file.write('A_lon = np.array([\n    [%f, %f, %f, %f, %f],\n    '
               '[%f, %f, %f, %f, %f],\n    '
               '[%f, %f, %f, %f, %f],\n    '
               '[%f, %f, %f, %f, %f],\n    '
               '[%f, %f, %f, %f, %f]])\n' %
    (A_lon[0][0], A_lon[0][1], A_lon[0][2], A_lon[0][3], A_lon[0][4],
     A_lon[1][0], A_lon[1][1], A_lon[1][2], A_lon[1][3], A_lon[1][4],
     A_lon[2][0], A_lon[2][1], A_lon[2][2], A_lon[2][3], A_lon[2][4],
     A_lon[3][0], A_lon[3][1], A_lon[3][2], A_lon[3][3], A_lon[3][4],
     A_lon[4][0], A_lon[4][1], A_lon[4][2], A_lon[4][3], A_lon[4][4]))
    file.write('B_lon = np.array([\n    [%f, %f],\n    '
               '[%f, %f],\n    '
               '[%f, %f],\n    '
               '[%f, %f],\n    '
               '[%f, %f]])\n' %
    (B_lon[0][0], B_lon[0][1],
     B_lon[1][0], B_lon[1][1],
     B_lon[2][0], B_lon[2][1],
     B_lon[3][0], B_lon[3][1],
     B_lon[4][0], B_lon[4][1],))
    file.write('A_lat = np.array([\n    [%f, %f, %f, %f, %f],\n    '
               '[%f, %f, %f, %f, %f],\n    '
               '[%f, %f, %f, %f, %f],\n    '
               '[%f, %f, %f, %f, %f],\n    '
               '[%f, %f, %f, %f, %f]])\n' %
    (A_lat[0][0], A_lat[0][1], A_lat[0][2], A_lat[0][3], A_lat[0][4],
     A_lat[1][0], A_lat[1][1], A_lat[1][2], A_lat[1][3], A_lat[1][4],
     A_lat[2][0], A_lat[2][1], A_lat[2][2], A_lat[2][3], A_lat[2][4],
     A_lat[3][0], A_lat[3][1], A_lat[3][2], A_lat[3][3], A_lat[3][4],
     A_lat[4][0], A_lat[4][1], A_lat[4][2], A_lat[4][3], A_lat[4][4]))
    file.write('B_lat = np.array([\n    [%f, %f],\n    '
               '[%f, %f],\n    '
               '[%f, %f],\n    '
               '[%f, %f],\n    '
               '[%f, %f]])\n' %
    (B_lat[0][0], B_lat[0][1],
     B_lat[1][0], B_lat[1][1],
     B_lat[2][0], B_lat[2][1],
     B_lat[3][0], B_lat[3][1],
     B_lat[4][0], B_lat[4][1],))
    file.write('Ts = %f\n' % Ts)
    file.close()


def compute_tf_model(mav, trim_state, trim_input):
    # trim values
    mav._state = trim_state
    mav._update_velocity_data()
    Va_trim = mav._Va
    alpha_trim = mav._alpha
    phi, theta_trim, psi = Quaternion2Euler(trim_state[6:10])

    # define transfer function constants
    #CHAP5 - slide 21
    a_phi1 = -(MAV.rho*Va_trim**2*MAV.S_wing*MAV.b/2) * MAV.C_p_p * MAV.b/(2*Va_trim)
    a_phi2 = (MAV.rho*Va_trim**2*MAV.S_wing*MAV.b/2) * MAV.C_p_delta_a
    #CHAP5 - slide 28
    a_theta1 = -(MAV.rho*Va_trim**2*MAV.c*MAV.S_wing/(2*MAV.Jy)) * MAV.C_m_q * MAV.c / (2*Va_trim)
    a_theta2 = -(MAV.rho*Va_trim**2*MAV.c*MAV.S_wing/(2*MAV.Jy)) * MAV.C_m_alpha
    a_theta3 = (MAV.rho*Va_trim**2*MAV.c*MAV.S_wing/(2*MAV.Jy)) * MAV.C_m_delta_e

    # Compute transfer function coefficients using new propulsion model
    # CHAP5 - slide 34
    a_V1 = MAV.rho*Va_trim*MAV.S_wing/MAV.mass * (MAV.C_D_0 + MAV.C_D_alpha*alpha_trim + MAV.C_D_delta_e*trim_input.elevator) \
        - (1/MAV.mass) * dT_dVa(mav, Va_trim, trim_input.throttle)
    a_V2 = (1/MAV.mass) * dT_ddelta_t(mav, Va_trim, trim_input.throttle)
    a_V3 = MAV.gravity*np.cos(theta_trim-alpha_trim)

    return Va_trim, alpha_trim, theta_trim, a_phi1, a_phi2, a_theta1, a_theta2, a_theta3, a_V1, a_V2, a_V3


def compute_ss_model(mav, trim_state, trim_input):
    x_euler = euler_state(trim_state)
    Aq = df_dx(mav, x_euler, trim_input)
    Bq = df_du(mav, x_euler, trim_input)

    T = dt_dq(mav, trim_state, trim_input)
    Tinv = dtI_dq(mav, trim_state, trim_input)

    A = T @ Aq @ Tinv
    A[2, :] = -A[2, :]

    B = T @ Bq
    B[2, :] = - B[2, :]

    # extract longitudinal states (u, w, q, theta, pd) and change pd to h
    #CHAP5 - slide 51

    E1 = np.array([
        [0.0, 0.0,  0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0,  0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ])

    E2 = np.array([
        [0.0, 1.0, 0.0, 0.0],
        [1.0, 0.0, 0.0, 0.0],
    ])

    A_lon = E1 @ A @ E1.T
    B_lon = E1 @ B @ E2.T

    # extract lateral states (v, p, r, phi, psi)
    #CHAP5 - slide 52

    E3 = np.array([
        [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
    ])

    E4 = np.array([
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ])

    A_lat = E3 @ A @ E3.T
    B_lat = E3 @ B @ E4.T
    return A_lon, B_lon, A_lat, B_lat

def euler_state(x_quat):
    # convert state x with attitude represented by quaternion
    # to x_euler with attitude represented by Euler angles

    #CHAP5 - slide 48
    e = np.array([x_quat.item(6), x_quat.item(7), x_quat.item(8), x_quat.item(9)])
    [phi, theta, psi] = Quaternion2Euler(e)
    x_euler = np.array([[x_quat.item(0)],
                        [x_quat.item(1)],
                        [x_quat.item(2)],
                        [x_quat.item(3)],
                        [x_quat.item(4)],
                        [x_quat.item(5)],
                        [phi],
                        [theta],
                        [psi],
                        [x_quat.item(10)],
                        [x_quat.item(11)],
                        [x_quat.item(12)]
                        ])

    return x_euler

def quaternion_state(x_euler):
    # convert state x_euler with attitude represented by Euler angles
    # to x_quat with attitude represented by quaternions

    #CHAP5 - slide 48
    e = Euler2Quaternion(x_euler.item(6), x_euler.item(7), x_euler.item(8))
    x_quat = np.array([[x_euler.item(0)],
                       [x_euler.item(1)],
                       [x_euler.item(2)],
                       [x_euler.item(3)],
                       [x_euler.item(4)],
                       [x_euler.item(5)],
                       [e.item(0)],
                       [e.item(1)],
                       [e.item(2)],
                       [e.item(3)],
                       [x_euler.item(9)],
                       [x_euler.item(10)],
                       [x_euler.item(11)]
                       ])

    return x_quat

def f_euler(mav, x_euler, delta):

    forces_moments = mav._forces_moments(delta)
    x_quat = quaternion_state(x_euler)
    mav._state = x_quat
    mav._update_velocity_data(np.zeros((6,1)))
    f_quat = mav._derivatives(x_quat, forces_moments)
    f_euler_ = euler_state(f_quat)
    return f_euler_

def df_dx(mav, x_euler, delta):
    # CHAP5 - slide 47

    # take partial of f_euler with respect to x_euler

    # A = jacobian(f_euler, mav, x_euler, delta, 0)
    # print("HAVE")
    # print(A)
    # return A

    m = 12
    n = 12
    eps = 0.001
    A = np.zeros((m, n))
    f_at_x = f_euler(mav, x_euler, delta)
    for i in range(0, n):
        x_eps = np.copy(x_euler)
        x_eps[i][0] += eps
        f_at_x_eps = f_euler(mav, x_eps, delta)
        df_dxi = (f_at_x_eps-f_at_x)/eps
        A[:, i] = df_dxi[:, 0]
    print("HELLO")
    print(A)
    return A


def df_du(mav, x_euler, delta):
    # CHAP5 - slide 47 (modified for u)
    # take partial of f_euler with respect to input
    m = 12
    n = 4
    eps = 0.001
    B = np.zeros((m, n))
    f_at_u = f_euler(mav, x_euler, delta)
    for i in range(0, n):
        u_eps = delta.to_array()
        u_eps[i] += eps
        delta_eps = MsgDelta()
        delta_eps.from_array(u_eps)
        f_at_u_eps = f_euler(mav, x_euler, delta_eps)
        df_dui = (f_at_u_eps-f_at_u)/eps
        B[:, i] = df_dui[:, 0]
    return B


def dT_dVa(mav, Va, delta_t):
    # returns the derivative of motor thrust with respect to Va
    eps = 0.001
    T_eps, Q_eps = mav._motor_thrust_torque(Va + eps, delta_t)
    T, Q = mav._motor_thrust_torque(Va, delta_t)
    return (T_eps - T) / eps

def dT_ddelta_t(mav, Va, delta_t):
    # returns the derivative of motor thrust with respect to delta_t
    print("SOMETHING", delta_t)
    eps = 0.001
    T_eps, Q_eps = mav._motor_thrust_torque(Va, delta_t + eps)
    T, Q = mav._motor_thrust_torque(Va, delta_t)
    print((T_eps - T)/eps/11)
    return (T_eps - T) / eps

def dt_dq(mav, x_euler, delta):
    dt_dq = np.zeros([12, 13])
    dt_dq[0:6, 0:6] = np.eye(6)
    dt_dq[9:, 10:] = np.eye(3)
    dt_dq[6:9, 6:10] = jacobian(Quaternion2Euler, mav, x_euler[6:10], delta, 2)
    return dt_dq

def dtI_dq(mav, x_euler, delta):
    [phi, theta, psi] = Quaternion2Euler(x_euler[6:10])
    dtI_dq = np.zeros([13, 12])
    dtI_dq[0:6, 0:6] = np.eye(6)
    dtI_dq[10:, 9:] = np.eye(3)
    dtI_dq[6:10, 6:9] = jacobian(Euler2Quaternion, mav, (phi, theta, psi), delta, 3)
    return dtI_dq

def jacobian(fun, mav, x, input, spot):
    # compute jacobian of fun with respect to x
    if spot == 0:
        f = fun(mav, x, input)
        m = f.shape[0]
        n = x.shape[0]
        eps = 0.01  # deviation
        A = np.zeros((m, n))
        for i in range(0, n):
            x_eps = np.copy(x)
            x_eps[i][0] += eps
            f_eps = fun(mav, x_eps, input)
            df = (f_eps - f) / eps
            A[:, i] = df[:, 0]
        return A
    elif spot == 1:
        f = fun(mav, x, input)
        m = f.shape[0]
        n = input.shape[0]
        eps = 0.01  # deviation
        B = np.zeros((m, n))
        for i in range(0, n):
            input_eps = np.copy(input)
            input_eps[i][0] += eps
            f_eps = fun(mav, x, input_eps)
            df = (f_eps - f) / eps
            B[:, i] = df[:, 0]
        return B
    elif spot == 2:
        f = np.asarray(fun(x))
        f = f.reshape([len(f), 1])
        m = f.shape[0]
        n = x.shape[0]
        eps = 0.01  # deviation
        dthde = np.zeros((m, n))
        for i in range(0, n):
            x_eps = np.copy(x)
            x_eps[i][0] += eps

            f_eps = np.asarray(fun(x_eps))
            f_eps = f_eps.reshape([len(f), 1])
            df = np.array([[(f_eps.item(0) - f.item(0)) / eps, (f_eps.item(1) - f.item(1)) / eps, (f_eps.item(2) - f.item(2)) / eps]]).T
            # df = df.reshape([3,1])
            dthde[:, i] = df[:, 0]
        return dthde
    elif spot == 3:
        phi = x[0]
        theta = x[1]
        psi = x[2]
        f = np.asarray(fun(phi, theta, psi))
        f = f.reshape([len(f), 1])
        m = f.shape[0]
        n = len(x)
        eps = 0.01  # deviation
        dthde = np.zeros((m, n))
        for i in range(0, n):
            x_eps = np.copy(x).reshape([1, 3])
            x_eps[0][i] += eps
            phi = x_eps.item(0)
            theta = x_eps.item(1)
            psi = x_eps.item(2)
            f_eps = np.asarray(fun(phi, theta, psi))
            f_eps = f_eps.reshape([len(f), 1])
            df = (f_eps - f) / eps
            # df = df.reshape([3,1])
            dthde[:, i] = df[:, 0]
        return dthde