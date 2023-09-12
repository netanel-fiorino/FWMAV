import sys
sys.path.append('..')
import numpy as np

# load message types
from message_types.msg_state import MsgState
import parameters.FWMAV_parameters as FWMAV
from tools.rotations import Quaternion2Euler, Quaternion2Rotation


class FWMAVDynamics:
    def __init__(self, Ts):
        self.ts_simulation = Ts
        # set initial states based on parameter file
        # _state is the 13x1 internal state of the aircraft that is being propagated:
        # _state = [pn, pe, pd, u, v, w, e0, e1, e2, e3, p, q, r]
        self._state = np.array([[FWMAV.north0],  # (0)
                               [FWMAV.east0],   # (1)
                               [FWMAV.down0],   # (2)
                               [FWMAV.u0],    # (3)
                               [FWMAV.v0],    # (4)
                               [FWMAV.w0],    # (5)
                               [FWMAV.e0],    # (6)
                               [FWMAV.e1],    # (7)
                               [FWMAV.e2],    # (8)
                               [FWMAV.e3],    # (9)
                               [FWMAV.p0],    # (10)
                               [FWMAV.q0],    # (11)
                               [FWMAV.r0]])   # (12)
        self.true_state = MsgState()
        self._wind = np.array([[0.], [0.], [0.]])  # wind in NED frame in meters/sec
        self._update_velocity_data()
        # store forces to avoid recalculation in the sensors function
        self._forces = np.array([[0.], [0.], [0.]])
        self._Va = FWMAV.Va0
        self._alpha = 0
        self._beta = 0
    ###################################
    # public functions
    
    #update states based on explecit forces
    def updateForce(self, forces_moments):
        '''
            Integrate the differential equations defining dynamics. 
            Inputs are the forces and moments on the aircraft.
            Ts is the time step between function calls.
        '''

        # Integrate ODE using Runge-Kutta RK4 algorithm
        time_step = self.ts_simulation
        k1 = self._derivatives(self._state, forces_moments)
        k2 = self._derivatives(self._state + time_step/2.*k1, forces_moments)
        k3 = self._derivatives(self._state + time_step/2.*k2, forces_moments)
        k4 = self._derivatives(self._state + time_step*k3, forces_moments)
        self._state += time_step/6 * (k1 + 2*k2 + 2*k3 + k4)

        # normalize the quaternion
        e0 = self._state.item(6)
        e1 = self._state.item(7)
        e2 = self._state.item(8)
        e3 = self._state.item(9)
        normE = np.sqrt(e0**2+e1**2+e2**2+e3**2)
        self._state[6][0] = self._state.item(6)/normE
        self._state[7][0] = self._state.item(7)/normE
        self._state[8][0] = self._state.item(8)/normE
        self._state[9][0] = self._state.item(9)/normE

        # update the message class for the true state
        self._update_true_state()

    #update states based on wind and control inputs
    def update(self, delta, wind):
        '''
            Integrate the differential equations defining dynamics. 
            Inputs are the forces and moments on the aircraft.
            Ts is the time step between function calls.
        '''
        #get forces and moments acting on the FWMAV
        forces_moments = self._forces_moments(delta)
        # Integrate ODE using Runge-Kutta RK4 algorithm
        time_step = self.ts_simulation
        k1 = self._derivatives(self._state, forces_moments)
        k2 = self._derivatives(self._state + time_step/2.*k1, forces_moments)
        k3 = self._derivatives(self._state + time_step/2.*k2, forces_moments)
        k4 = self._derivatives(self._state + time_step*k3, forces_moments)
        self._state += time_step/6 * (k1 + 2*k2 + 2*k3 + k4)

        # normalize the quaternion
        e0 = self._state.item(6)
        e1 = self._state.item(7)
        e2 = self._state.item(8)
        e3 = self._state.item(9)
        normE = np.sqrt(e0**2+e1**2+e2**2+e3**2)
        self._state[6][0] = self._state.item(6)/normE
        self._state[7][0] = self._state.item(7)/normE
        self._state[8][0] = self._state.item(8)/normE
        self._state[9][0] = self._state.item(9)/normE

        # update the airspeed, angle of attack, and side slip angles using new state
        self._update_velocity_data(wind)
        # update the message class for the true state
        self._update_true_state()

    def external_set_state(self, new_state):
        self._state = new_state

    ###################################
    # private functions
    def _derivatives(self, state, forces_moments):
        """
        for the dynamics xdot = f(x, u), returns f(x, u)
        """
        # extract the states
        north = state.item(0)
        east = state.item(1)
        down = state.item(2)
        u = state.item(3)
        v = state.item(4)
        w = state.item(5)
        e0 = state.item(6)
        e1 = state.item(7)
        e2 = state.item(8)
        e3 = state.item(9)
        p = state.item(10)
        q = state.item(11)
        r = state.item(12)
        #   extract forces/moments
        fx = forces_moments.item(0)
        fy = forces_moments.item(1)
        fz = forces_moments.item(2)
        l = forces_moments.item(3)
        m = forces_moments.item(4)
        n = forces_moments.item(5)

        # position kinematics
        pos_dot = Quaternion2Rotation(np.array([e0, e1, e2, e3])) @ np.array([u, v, w])     #THIS MIGHT NEED TO BE Quaternion2Rotation(np.array([e0, e1, e2, e3])) @ np.array([u, v, w]).T
        north_dot = pos_dot[0]
        east_dot = pos_dot[1]
        down_dot = pos_dot[2]

        # position dynamics
        u_dot = r * v - q * w + fx / FWMAV.mass
        v_dot = p * w - r * u + fy / FWMAV.mass
        w_dot = q * u - p * v + fz / FWMAV.mass

        # rotational kinematics
        # skew = np.array([
        #     [0, -p, -q, -r], 
        #     [p, 0, r, -q],
        #     [q, -r, 0, p],
        #     [r, q, -p, 0]])
        # e_dot = 0.5 * skew @ np.array([e0, e1, e2, e3]).T
        # e0_dot = e_dot.item(0) 
        # e1_dot = e_dot.item(1) 
        # e2_dot = e_dot.item(2) 
        # e3_dot = e_dot.item(3) 
        e0_dot = 0.5 * (-p*e1 - q*e2 - r*e3)
        e1_dot = 0.5 * ( p*e0 + r*e2 - q*e3)
        e2_dot = 0.5 * ( q*e0 - r*e1 + p*e3)
        e3_dot = 0.5 * ( r*e0 + q*e1 - p*e2)

        # rotatonal dynamics
        p_dot = FWMAV.gamma1*p*q - FWMAV.gamma2*q*r + FWMAV.gamma3*l + FWMAV.gamma4*n
        q_dot = FWMAV.gamma5*p*r - FWMAV.gamma6*(p**2 - r**2) + m / FWMAV.Jy
        r_dot = FWMAV.gamma7*p*q - FWMAV.gamma1*q*r + FWMAV.gamma4*l + FWMAV.gamma8*n

        # collect the derivative of the states
        x_dot = np.array([[north_dot, east_dot, down_dot, u_dot, v_dot, w_dot,
                           e0_dot, e1_dot, e2_dot, e3_dot, p_dot, q_dot, r_dot]]).T
        return x_dot

    def _update_velocity_data(self, wind=np.zeros((6,1))):
        steady_state = wind[0:3]
        gust = wind[3:6]
        Ri_b = Quaternion2Rotation(self._state[6:10]).T
        wind_body_frame = Ri_b @ steady_state + gust
        # velocity vector relative to the airmass
        v_air = self._state[3:6] - wind_body_frame
        ur = v_air.item(0)
        vr = v_air.item(1)
        wr = v_air.item(2)
        # compute airspeed
        self._Va = np.linalg.norm(v_air)
        # compute angle of attack
        if ur == 0:
            self._alpha = 0
        else:
            self._alpha = np.arctan(wr/ur)
        # compute sideslip angle
        if self._Va == 0:
            self._beta = 0
        else:
            self._beta = np.arcsin(vr/self._Va)

    def _forces_moments(self, delta):
        """
        return the forces on the UAV based on the state, wind, and control surfaces
        :param delta: np.matrix(delta_a, delta_e, delta_r, delta_t)
        :return: Forces and Moments on the UAV np.matrix(Fx, Fy, Fz, Ml, Mn, Mm)
        """
        # phi, theta, psi = Quaternion2Euler(self._state[6:10])
        p = self._state.item(10)
        q = self._state.item(11)
        r = self._state.item(12)

        # Sigma --> Blending Function
        sigma_num = (1 + np.exp(-FWMAV.M*(self._alpha - FWMAV.alpha0)) + np.exp(FWMAV.M*(self._alpha + FWMAV.alpha0))) 
        sigma_denom = (1 + np.exp(-FWMAV.M*(self._alpha - FWMAV.alpha0))) * (1 + np.exp(FWMAV.M*(self._alpha + FWMAV.alpha0)))
        sigma = sigma_num / sigma_denom

        # compute gravitaional forces
        Ri_b = Quaternion2Rotation(self._state[6:10]).T
        f_g = Ri_b @ np.array([
                            [0],
                            [0],
                            [FWMAV.mass*FWMAV.gravity],
                            ])

        # compute Lift and Drag coefficients
        CL = ((1 - sigma)*(FWMAV.C_L_0 + FWMAV.C_L_alpha*self._alpha)) + ((sigma)*2*np.sign(self._alpha)*(np.sin(self._alpha)**2)*np.cos(self._alpha))
        CD = FWMAV.C_D_p + ((FWMAV.C_L_0 + FWMAV.C_L_alpha*self._alpha)**2)/(np.pi*FWMAV.e*FWMAV.AR)

        # compute Lift and Drag Forces
        q_dynamic = 0.5*FWMAV.rho*(self._Va**2)*FWMAV.S_wing        # dynamic pressure --> force term (*S) for lift/drag
        F_lift = q_dynamic*(CL + (FWMAV.C_L_q*0.5*FWMAV.c/self._Va*q) + FWMAV.C_L_delta_e*delta.elevator)
        F_drag = q_dynamic*(CD + (FWMAV.C_D_q*0.5*FWMAV.c/self._Va*q) + FWMAV.C_D_delta_e*delta.elevator)

        #compute propeller thrust and torque
        thrust_prop, torque_prop = self._motor_thrust_torque(self._Va, delta.throttle)

        # compute longitudinal forces in body frame
        #fx = -F_drag*np.cos(self._alpha) + F_lift*np.sin(self._alpha) + f_g.item(0) + thrust_prop
        #fz = -F_drag*np.sin(self._alpha) - F_lift*np.cos(self._alpha) + f_g.item(2)

        fx = -F_drag*np.cos(self._alpha) + F_lift*np.cos(delta.kappa)*np.sin(self._alpha) + f_g.item(0) + thrust_prop
        fz = -F_drag*np.sin(self._alpha) - F_lift*np.cos(delta.kappa)*np.cos(self._alpha) + f_g.item(2)

        # compute lateral forces in body frame
        fy = q_dynamic * (FWMAV.C_Y_0 + FWMAV.C_Y_beta*self._beta + FWMAV.C_Y_p*0.5*FWMAV.b*p/self._Va + FWMAV.C_Y_r*0.5*FWMAV.b*r/self._Va + FWMAV.C_Y_delta_a*delta.aileron + FWMAV.C_Y_delta_r*delta.rudder) \
            + f_g.item(1)

        # compute logitudinal torque in body frame
        My = q_dynamic*FWMAV.c*(FWMAV.C_m_0 + (FWMAV.C_m_alpha*self._alpha) + (FWMAV.C_m_q*FWMAV.c*0.5*q/self._Va) + (FWMAV.C_m_delta_e*delta.elevator))
       
        # compute lateral torques in body frame
        Mx = q_dynamic*FWMAV.b*(FWMAV.C_ell_0 + FWMAV.C_ell_beta*self._beta + FWMAV.C_ell_p*0.5*FWMAV.b*p/self._Va + FWMAV.C_ell_r*0.5*FWMAV.b*r/self._Va + FWMAV.C_ell_delta_a*delta.aileron + FWMAV.C_ell_delta_r*delta.rudder) \
            + torque_prop
        Mz = q_dynamic*FWMAV.b*(FWMAV.C_n_0 + FWMAV.C_n_beta*self._beta + FWMAV.C_n_p*0.5*FWMAV.b*p/self._Va + FWMAV.C_n_r*0.5*FWMAV.b*r/self._Va + FWMAV.C_n_delta_a*delta.aileron + FWMAV.C_n_delta_r*delta.rudder)

        self._forces[0] = fx
        self._forces[1] = fy 
        self._forces[2] = fz
        #print(np.array([[fx, fy, fz, Mx, My, Mz]])) 
        return np.array([[fx, fy, fz, Mx, My, Mz]]).T

    def _motor_thrust_torque(self, Va, delta_t):
        # compute thrust and torque due to propeller  (See addendum by McLain)

        # map delta_t throttle command(0 to 1) into motor input voltage
        V_in = FWMAV.V_max * delta_t

        # quadratic formula to solve for motor speed
        a = FWMAV.C_Q0 * FWMAV.rho * np.power(FWMAV.D_prop, 5) / ((2*np.pi)**2)
        b = (FWMAV.C_Q1 * FWMAV.rho * np.power(FWMAV.D_prop, 4) * Va / (2*np.pi)) + (FWMAV.KQ**2)/FWMAV.R_motor
        c = (FWMAV.C_Q2 * FWMAV.rho * np.power(FWMAV.D_prop, 3) * Va**2) - (FWMAV.KQ/FWMAV.R_motor)*V_in + FWMAV.KQ*FWMAV.i0

        # operating propeller speed
        Omega_op = (-b + np.sqrt(b**2 - 4*a*c)) / (2.*a) # rad/sec
        # compute advance ratio
        J_op = 2 * np.pi * Va / (Omega_op * FWMAV.D_prop)

        # thrust and torque coefficients
        C_T = FWMAV.C_T2 * J_op**2 + FWMAV.C_T1 * J_op + FWMAV.C_T0
        C_Q = FWMAV.C_Q2 * J_op**2 + FWMAV.C_Q1 * J_op + FWMAV.C_Q0

        n = Omega_op / (2 * np.pi)  # rev/sec

        # thrust and torque due to propeller
        thrust_prop = FWMAV.rho * n**2 * np.power(FWMAV.D_prop, 4) * C_T
        torque_prop = - FWMAV.rho * n**2 * np.power(FWMAV.D_prop, 5) * C_Q
        return 0,0 #thrust_prop, torque_prop

    def _update_true_state(self):
        # update the true state message:
        phi, theta, psi = Quaternion2Euler(self._state[6:10])
        pdot = Quaternion2Rotation(self._state[6:10]) @ self._state[3:6]
        self.true_state.north = self._state.item(0)
        self.true_state.east = self._state.item(1)
        self.true_state.altitude = -self._state.item(2)
        self.true_state.Va = self._Va
        self.true_state.alpha = self._alpha
        self.true_state.beta = self._beta
        self.true_state.phi = phi
        self.true_state.theta = theta
        self.true_state.psi = psi
        self.true_state.Vg = np.linalg.norm(pdot)
        self.true_state.gamma = np.arcsin(pdot.item(2) / self.true_state.Vg)
        self.true_state.chi = np.arctan2(pdot.item(1), pdot.item(0))
        self.true_state.p = self._state.item(10)
        self.true_state.q = self._state.item(11)
        self.true_state.r = self._state.item(12)
        
