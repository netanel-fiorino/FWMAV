import sys
sys.path.append('..')
import numpy as np

# load message types
from message_types.msg_state import MsgState
import parameters.FWMAV_parameters as FWMAV
from tools.rotations import Quaternion2Euler, Quaternion2Rotation, skew
from wingFlapping.wingDynamics import WingDynamics, R_theta, R_psi, R_phi, R_beta

class FWMAVDynamics:
    def __init__(self, Ts, wingLeft, wingRight):
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
        self._update_velocity_data_init()
        # store forces to avoid recalculation in the sensors function
        self._forces = np.array([[0.], [0.], [0.]])
        self._Va = FWMAV.Va0
        self._Va_tail = FWMAV.Va0 + np.linalg.norm(skew(self._state[10:13]) @ FWMAV.r_tail)
        self.wing_left = wingLeft
        self.wing_right = wingRight
        #self._alpha = 0
        #self._alpha_tail = 0 - FWMAV.tail_angle
        #self._beta = 0
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
    def update(self, delta, wind, time):
        '''
            Integrate the differential equations defining dynamics. 
            Inputs are the forces and moments on the aircraft.
            Ts is the time step between function calls.
        '''
        #get forces and moments acting on the FWMAV
        forces_moments = self._forces_moments(delta, time)
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
        self._update_velocity_data(delta, wind)
        #self._update_velocity_data_init(wind)
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
        phi, theta, psi = Quaternion2Euler(np.array([e0,e1,e2,e3]))
        L_inv = np.array([
            [1, np.sin(phi)*np.tan(theta), np.cos(phi)*np.tan(theta)],
            [0, np.cos(phi), -np.sin(phi)],
            [0, np.sin(phi)*1/np.cos(theta), np.cos(phi)*1/np.cos(theta)]
            ])
        omega = L_inv * np.array([p,q,r])
        J = np.array([
            [FWMAV.Jx, 0, FWMAV.Jxz],
            [0, FWMAV.Jy, 0],
            [FWMAV.Jxz, 0, FWMAV.Jz]
        ])
        #omega_dot = np.linalg.inv(J) @ (-(np.array([p,q,r])) @ (J*np.array([p,q,r])).T + np.array([l,m,n]))
        omega_dot = np.linalg.inv(J).dot(np.cross(-(np.array([p,q,r])), J.dot(np.array([p,q,r]))) + np.array([l,m,n]))
        p_dot = omega_dot[0]
        q_dot = omega_dot[1]
        r_dot = omega_dot[2]

        #print(q_dot, FWMAV.gamma5*p*r - FWMAV.gamma6*(p**2 - r**2) + m / FWMAV.Jy)
        #print(p_dot, FWMAV.gamma1*p*q - FWMAV.gamma2*q*r + FWMAV.gamma3*l + FWMAV.gamma4*n)
        #print(r_dot, FWMAV.gamma7*p*q - FWMAV.gamma1*q*r + FWMAV.gamma4*l + FWMAV.gamma8*n)
        # p_dot = FWMAV.gamma1*p*q - FWMAV.gamma2*q*r + FWMAV.gamma3*l + FWMAV.gamma4*n
        # q_dot = FWMAV.gamma5*p*r - FWMAV.gamma6*(p**2 - r**2) + m / FWMAV.Jy
        # r_dot = FWMAV.gamma7*p*q - FWMAV.gamma1*q*r + FWMAV.gamma4*l + FWMAV.gamma8*n

        # collect the derivative of the states
        x_dot = np.array([[north_dot, east_dot, down_dot, u_dot, v_dot, w_dot,
                           e0_dot, e1_dot, e2_dot, e3_dot, p_dot, q_dot, r_dot]]).T
        return x_dot

    def _update_velocity_data_init(self, wind=np.zeros((6,1))):
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
        v_air_tail = v_air + skew(self._state[10:13]) @ FWMAV.r_tail
        #self._Va_tail = np.linalg.norm(v_air_tail)[0]
        self._Va_tail = self._Va
        #print(self._Va_tail)
        # compute angle of attack
        if ur == 0:
            self._alpha = 0
        else:
            self._alpha = np.arctan(wr/ur)
        if v_air_tail.item(0) == 0:
            self._alpha_tail = 0
        else:
            #self._alpha_tail = np.arctan(wr/ur) - FWMAV.tail_angle
            self._alpha_tail = np.arctan(v_air_tail[2,0].item(0)/v_air_tail[0,0].item(0)) - FWMAV.tail_angle
        # compute sideslip angle
        if self._Va == 0:
            self._beta = 0
        else:
            self._beta = np.arcsin(vr/self._Va)
        

    def _update_velocity_data(self, delta, wind=np.zeros((6,1))):
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
        v_air_tail = v_air.T + skew(self._state[10:13]) @ FWMAV.r_tail
        #print(np.linalg.norm(v_air_tail))
        #self._Va_tail = np.linalg.norm(v_air_tail)[0]
        self._Va_tail = np.linalg.norm(v_air_tail)
        # compute angle of attack
        if ur == 0:
            self._alpha = 0
        else:
            self._alpha = np.arctan(wr/ur)
        if v_air_tail.item(0) == 0:
            self._alpha_tail = 0
        else:
            self._alpha_tail = np.arctan(v_air_tail[0,2].item(0)/v_air_tail[0,0].item(0)) - delta.tail_angle
        # compute sideslip angle
        if self._Va == 0:
            self._beta = 0
        else:
            self._beta = np.arcsin(vr/self._Va)
        #print("update")

    def _forces_moments(self, delta, time):
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
        F_lift = q_dynamic*(CL + (FWMAV.C_L_q*0.5*FWMAV.c/self._Va*q))
        F_drag = q_dynamic*(CD + (FWMAV.C_D_q*0.5*FWMAV.c/self._Va*q))

        #sigma_num_tail = (1 + np.exp(-FWMAV.M*(self._alpha_tail - FWMAV.alpha0)) + np.exp(FWMAV.M*(self._alpha_tail + FWMAV.alpha0))) 
        #sigma_denom_tail = (1 + np.exp(-FWMAV.M*(self._alpha_tail - FWMAV.alpha0))) * (1 + np.exp(FWMAV.M*(self._alpha_tail + FWMAV.alpha0)))
        #sigma_tail = sigma_num_tail / sigma_denom_tail

        #CL_tail = (2*np.pi*self._alpha_tail)/(1+(2*np.pi*self._alpha_tail)/(np.pi*FWMAV.AR)) #((1 - sigma_tail)*(FWMAV.C_L_0 + FWMAV.C_L_alpha*self._alpha_tail)) + ((sigma)*2*np.sign(self._alpha_tail)*(np.sin(self._alpha_tail)**2)*np.cos(self._alpha_tail))
        #CD_tail = 1.28*np.abs(np.sin(self._alpha_tail)) + CL_tail**2/(np.pi*FWMAV.AR*FWMAV.e) #FWMAV.C_D_p + ((FWMAV.C_L_0 + FWMAV.C_L_alpha*self._alpha_tail)**2)/(np.pi*FWMAV.e*FWMAV.AR)

        a0=7.159
        a1=-3.741
        b1=-0.08485
        a2=-0.7499
        b2=0.02863
        w=1.99
        c_d=a0+a1*np.cos(self._alpha_tail*w)+b1*np.sin(self._alpha_tail*w)+a2*np.cos(2*self._alpha_tail*w)+b2*np.sin(2*self._alpha_tail*w)
        CD_tail = c_d/3.8

        a0_2=-0.2495
        a1_2=-0.02425
        b1_2=5.896
        a2_2=0.1975
        b2_2=1.429
        a3_2=0.01897
        b3_2=0.4558
        w_2=2.011
        c_l=a0_2+a1_2*np.cos(self._alpha_tail*w_2)+b1_2*np.sin(self._alpha_tail*w_2)+a2_2*np.cos(2*self._alpha_tail*w_2)+b2_2*np.sin(2*self._alpha_tail*w_2)+a3_2*np.cos(3*self._alpha_tail*w_2)+b3_2*np.sin(3*self._alpha_tail*w_2)
        CL_tail = c_l/5.02

        a0=6.411
        a1=0.5402
        b1=-0.8709
        a2=-4.047
        b2=-0.1551
        w=1.654
        c_d=a0+a1*np.cos(self._alpha*w)+b1*np.sin(self._alpha*w)+a2*np.cos(2*self._alpha*w)+b2*np.sin(2*self._alpha*w)
        c_l=-1.045-0.7912*np.cos(1.835*self._alpha)+4.43*np.sin(1.835*self._alpha)
        CD_body = c_d
        CL_body = c_l


        q_dynamic_tail = 0.5*FWMAV.rho*(self._Va_tail**2)*FWMAV.S_tail        # dynamic pressure --> force term (*S) for lift/drag
        F_lift_tail = q_dynamic_tail*(CL_tail) #q_dynamic_tail*(CL_tail + (FWMAV.C_L_q*0.5*FWMAV.c/self._Va_tail*q))
        F_drag_tail = q_dynamic_tail*(CD_tail) #q_dynamic_tail*(CD_tail + (FWMAV.C_D_q*0.5*FWMAV.c/self._Va_tail*q))

        q_dynamic_body = 0.5*FWMAV.rho*(self._Va**2)*FWMAV.S_body        # dynamic pressure --> force term (*S) for lift/drag
        F_lift_body = q_dynamic_body*(CL_body) #q_dynamic_tail*(CL_tail + (FWMAV.C_L_q*0.5*FWMAV.c/self._Va_tail*q))
        F_drag_body = q_dynamic_body*(CD_body) #q_dynamic_tail*(CD_tail + (FWMAV.C_D_q*0.5*FWMAV.c/self._Va_tail*q))
        
        q_dynamic_wing = 0.5*FWMAV.rho*(self._Va**2)*FWMAV.S_wing       # dynamic pressure --> force term (*S) for lift/drag
        F_lift_wing = q_dynamic_wing*(CL_tail) #q_dynamic_tail*(CL_tail + (FWMAV.C_L_q*0.5*FWMAV.c/self._Va_tail*q))
        F_drag_wing = q_dynamic_wing*(CD_tail) #q_dynamic_tail*(CD_tail + (FWMAV.C_D_q*0.5*FWMAV.c/self._Va_tail*q))

        # compute longitudinal forces in body frame
        #fx = -F_drag*np.cos(self._alpha) + F_lift*np.sin(self._alpha) + f_g.item(0) + thrust_prop
        #fz = -F_drag*np.sin(self._alpha) - F_lift*np.cos(self._alpha) + f_g.item(2)

        if delta.flap == 0:
            force_tail =  np.array([[np.cos(self._alpha_tail), -np.sin(self._alpha_tail)], [0,0], [np.sin(self._alpha_tail), np.cos(self._alpha_tail)]]) @ np.array([-F_drag_tail, -F_lift_tail])
            force_body = np.array([[np.cos(self._alpha), -np.sin(self._alpha)], [0,0], [np.sin(self._alpha), np.cos(self._alpha)]]) @ np.array([-F_drag_body, -F_lift_body])
            force_wing = np.array([[np.cos(self._alpha), -np.sin(self._alpha)], [0,0], [np.sin(self._alpha), np.cos(self._alpha)]]) @ np.array([-F_drag_wing, -F_lift_wing])
            fx = force_wing.item(0) + f_g.item(0) + force_tail.item(0) + force_body.item(0) + 0.1633
            fy = force_wing.item(1) + f_g.item(1) + force_tail.item(1) + force_body.item(1) 
            fz = force_wing.item(2) + f_g.item(2) + force_tail.item(2) + force_body.item(2)
            #print(np.linalg.norm(self._state[3:6] + skew(self._state[10:13]) @ FWMAV.r_tail), self._Va, self._state[3:6] + skew(self._state[10:13]) @ FWMAV.r_tail, self._state[3:6])
            #print(force_wing.item(0), f_g.item(0), force_tail.item(0), force_body.item(0))
            #print(self._Va, self._Va_tail)
            #print(50*"-")
            M = np.cross(FWMAV.r_tail, force_tail.T) + np.cross(FWMAV.r_body, force_body.T) + np.cross(FWMAV.r_wing, force_wing.T)
            #print(force_tail.item(2), force_body.item(2), (180/(2*np.pi)*self._alpha), force_wing.item(2))
            My = M.item(1)
            #print(My, np.cross(FWMAV.r_tail, force_tail.T), np.cross(FWMAV.r_body, force_body.T))
            Mx = 0
            Mz = 0
            
            #fx = -F_drag*np.cos(self._alpha) + F_lift*np.sin(self._alpha) + f_g.item(0) + -F_drag_tail*np.cos(self._alpha_tail) + F_lift_tail*np.sin(self._alpha_tail)
            #fz = -F_drag*np.sin(self._alpha) - F_lift*np.cos(self._alpha) + f_g.item(2) + -F_drag_tail*np.sin(self._alpha_tail) - F_lift_tail*np.cos(self._alpha_tail)
            #fy = q_dynamic * (FWMAV.C_Y_0 + FWMAV.C_Y_beta*self._beta + FWMAV.C_Y_p*0.5*FWMAV.b*p/self._Va + FWMAV.C_Y_r*0.5*FWMAV.b*r/self._Va + FWMAV.C_Y_delta_r*delta.rudder) \
            #+ f_g.item(1)
            # compute logitudinal torque in body frame
            #My = q_dynamic*FWMAV.c *(FWMAV.C_m_0 + (FWMAV.C_m_alpha*self._alpha) + (FWMAV.C_m_q*FWMAV.c*0.5*q/self._Va))
            #My = q_dynamic*FWMAV.c *(FWMAV.C_m_0 + (FWMAV.C_m_alpha*self._alpha) + (FWMAV.C_m_q*FWMAV.c*0.5*q/self._Va))
       
            # compute lateral torques in body frame
            #Mx = q_dynamic*FWMAV.b*(FWMAV.C_ell_0 + FWMAV.C_ell_beta*self._beta + FWMAV.C_ell_p*0.5*FWMAV.b*p/self._Va + FWMAV.C_ell_r*0.5*FWMAV.b*r/self._Va + FWMAV.C_ell_delta_r*delta.rudder) \
            #Mz = q_dynamic*FWMAV.b*(FWMAV.C_n_0 + FWMAV.C_n_beta*self._beta + FWMAV.C_n_p*0.5*FWMAV.b*p/self._Va + FWMAV.C_n_r*0.5*FWMAV.b*r/self._Va + FWMAV.C_n_delta_r*delta.rudder)
        else:
            thrust_angle = 0.122
            force_tail = np.array([[np.cos(self._alpha_tail), -np.sin(self._alpha_tail)], [0,0], [np.sin(self._alpha_tail), np.cos(self._alpha_tail)]]) @ np.array([-F_drag_tail, -F_lift_tail])
            force_body = np.array([[np.cos(self._alpha), -np.sin(self._alpha)], [0,0], [np.sin(self._alpha), np.cos(self._alpha)]]) @ np.array([-F_drag_body, -F_lift_body])
            force_wing = 1.1*(0.1*np.sin(2*np.pi*time)+0.055) * np.array([np.cos(thrust_angle), 0, np.sin(thrust_angle)])
            fx = force_wing.item(0) + f_g.item(0) + force_tail.item(0) + force_body.item(0)
            fy = force_wing.item(1) + f_g.item(1) + force_tail.item(1) + force_body.item(1)
            fz = force_wing.item(2) + f_g.item(2) + force_tail.item(2) + force_body.item(2)
            #fx = 1.1*(0.1*np.sin(2*np.pi*2*time)+0.055)*np.cos(0.61) + f_g.item(0) + -F_drag_tail*np.cos(self._alpha_tail) + F_lift_tail*np.sin(self._alpha_tail) + -F_drag_body*np.cos(self._alpha) + F_lift_body*np.sin(self._alpha)#+ aerodynamic forces from body and tail
            #fy = 0
            #fz = -1.1*(0.1*np.sin(2*np.pi*2*time)+0.055)*np.sin(0.61) + f_g.item(2) + -F_drag_tail*np.sin(self._alpha_tail) - F_lift_tail*np.cos(self._alpha_tail) + -F_drag_body*np.sin(self._alpha) - F_lift_body*np.cos(self._alpha)#+ aerodynamic forces from body and tail
            #My = FWMAV.wing_COP * -1*-1.1*(0.1*np.sin(2*np.pi*2*time)+0.055)*np.sin(thrust_angle) + FWMAV.tail_length * (-F_drag_tail*np.sin(self._alpha_tail) - F_lift_tail*np.cos(self._alpha_tail)) - 0.0069*(-F_drag_tail*np.cos(self._alpha_tail) + F_lift_tail*np.sin(self._alpha_tail)) + 0.0114*(-F_drag_body*np.sin(self._alpha) - F_lift_body*np.cos(self._alpha))#+ aerodynamic forces from body and tail 
            
            M = (np.cross(FWMAV.r_body, force_body.T) + np.cross(FWMAV.r_tail, force_tail.T) + np.cross(FWMAV.r_wing, force_wing.T))
            #print(M.item(1))
            #print(force_tail.item(2) + force_body.item(2))
            My = M.item(1)
            Mx = 0#q_dynamic*FWMAV.b * (FWMAV.C_ell_0 + FWMAV.C_ell_beta*self._beta + FWMAV.C_ell_p*0.5*FWMAV.b*p/self._Va + FWMAV.C_ell_r*0.5*FWMAV.b*r/self._Va)
            Mz = 0#q_dynamic*FWMAV.b * (FWMAV.C_n_0 + FWMAV.C_n_beta*self._beta + FWMAV.C_n_p*0.5*FWMAV.b*p/self._Va)        

        self._forces[0] = fx
        self._forces[1] = fy 
        self._forces[2] = fz
        #print(fx, fz, My, np.degrees(self._alpha), np.degrees(self._alpha_tail))
        #return np.array([[fx, fy, fz, 0, 0, 0]]).T
        return np.array([[fx, fy, fz, Mx, My, Mz]]).T

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
        
