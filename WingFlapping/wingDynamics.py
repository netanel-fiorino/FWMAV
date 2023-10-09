import numpy as np
import matplotlib.pyplot as plt
import wing_parameters as sim
from wing_parameters import phi_dot_f, psi_dot_f, theta_dot_f
plt.close('all')

class WingDynamics:
    def __init__(self, side, v0=sim.v_free_stream, wing_length=sim.wing_length, max_chord=sim.max_chord, freq=sim.freq, \
                amplitude=sim.amplitude, num_elements=sim.num_elements, wing_shape=sim.wing_shape, \
                beta=sim.beta, phi=sim.phi, psi=sim.psi, theta=sim.theta, phi_dot=phi_dot_f(0), psi_dot=psi_dot_f(0), theta_dot=theta_dot_f(0)):
        self._state = np.array([[0],    # (x-pos) [0]
                               [0],     # (y-pos) [1]
                               [0],     # (z-pos) [2]
                               [v0.item(0)],     # (u) [3]
                               [v0.item(1)],     # (v) [4]
                               [v0.item(2)],     # (w) [5]
                               [beta],     # (beta) angle to stroke plane [6]
                               [phi],     # (phi) causes flapping motion [7]
                               [psi],     # (psi) deviation angle [8]
                               [theta],      # (theta) [9]
                               [phi_dot],      # (phi_dot) [10]
                               [psi_dot],      # (psi_dot) [11]
                               [theta_dot]])     # (theta_dot) [12]
        self.num_elements = num_elements
        self.wing_length = wing_length
        self.wing_shape = wing_shape
        self.max_chord = max_chord
        self.avg_chord()
        self.freq = freq
        self.amplitude = amplitude
        self.u_inf = np.sqrt(self._state.item(3)**2 + self._state.item(4)**2 + self._state.item(5)**2)
        self.J = self.u_inf/(2*self.amplitude*self.wing_length*self.freq)
        self.side = side
        self.rhat2()
        self.rhatM()
        
        

    def chord_f(self, r):
        if self.wing_shape == 'rectangle':
            return self.max_chord                                            #rectangle wing
        elif self.wing_shape == 'ellipse':
            #print(1-r**2/(self.wing_length**2))
            return self.max_chord*np.sqrt(1-r**2/(self.wing_length**2))           #ellipse wing 
    
    def avg_chord(self):
        avg_chord = 0
        for i in range(self.num_elements):
            avg_chord = avg_chord + self.chord_f((i+0.5)*self.wing_length/self.num_elements)/self.num_elements
        self.avg_chord = avg_chord        #m      avg chord length
        #print(avg_chord)
    
    def rhat2(self):
        rhat2 = 0
        for i in range(1000):
            deltar = self.wing_length/1000/self.wing_length
            rhat2 = rhat2 + ((i+1)*deltar)**2*self.chord_f((i+1)*deltar*self.wing_length)/self.avg_chord * deltar
        self.rhat2 = np.sqrt(rhat2)
        #print(rhat2)

    def rhatM(self):
        rhatM = 0
        for i in range(1000):
            deltar = self.wing_length/1000/self.wing_length
            rhatM = rhatM + ((i+1)*deltar)**2*(self.chord_f((i+1)*deltar*self.wing_length)/self.avg_chord)**2 * deltar
        self.rhatM = np.sqrt(rhatM)
        #print(rhatM)


    def update(self, time):
        self._state[10][0] = phi_dot_f(time)
        self._state[11][0] = psi_dot_f(time)
        self._state[12][0] = theta_dot_f(time)

        self._state[7][0] = self._state.item(7) + self._state.item(10)*sim.t_step
        self._state[8][0] = self._state.item(8) + self._state.item(11)*sim.t_step
        self._state[9][0] = self._state.item(9) + self._state.item(12)*sim.t_step
        
    def updateOmega(self, omega):
        self._state[10][0] = omega
        
    
    def get_aero_coeff(self, alpha):
        self.u_inf = np.sqrt(self._state.item(3)**2 + self._state.item(4)**2 + self._state.item(5)**2)
        self.J = self.u_inf/(2*self.amplitude*self.wing_length*self.freq)
        xhat_0 = 0.25
        #rhat_2**2 = integral 1_0(rhat**2*chat*d(rhat))        rhat = r/R_wing, chat = c(r)/mean(c)
        #rhat_2 = np.sqrt(1/3)                   #non dimensional second moment of the wing
        rhat_2 = self.rhat2
        #rhat_M**2 = integral 1_0(rhat**2*chat**2*d(rhat))        rhat = r/R_wing, chat = c(r)/mean(c)
        #rhat_M = np.sqrt(1/3)
        rhat_M = self.rhatM
        K_PL = -2.109*((self.J + rhat_2)**-0.606) + 4.136
        K_VL = 2.659*((self.J + rhat_2)**-0.666) + -0.344
        K_PD = -0.182*((self.J + rhat_2)**-2.414) + 1.370
        K_VD = 0.765*((self.J + rhat_2)**-1.497) + 2.078
        K_PM = 0.803*((self.J + rhat_M)**-0.972) + -0.363
        K_VM = -0.242*((self.J + rhat_M)**-1.354) + -0.554


        CL = K_PL*np.sin(alpha)*(np.cos(alpha)**2) + K_VL*(np.sin(alpha)**2)*np.cos(alpha)
        CD = K_PD*(np.sin(alpha)**2)*np.cos(alpha) + K_VD*(np.sin(alpha)**3)
        CM = K_PM*(np.sin(alpha)**2)*np.cos(alpha) + K_VM*(np.sin(alpha)**2)
        CR = np.pi*(0.75-xhat_0)
        return CL, CD, CM, CR

    def force_calc(self):
        beta = self._state.item(6)
        phi = self._state.item(7)
        psi = self._state.item(8)
        theta = self._state.item(9)
        phi_dot = self._state.item(10)
        psi_dot = self._state.item(11)
        theta_dot = self._state.item(12)

        if self.side == 'right':              #side when looking at fwmav, (side along negative y axis)
            #Dont know if all these should be negatives
            phi = -phi
            phi_dot = -phi_dot
            theta = theta
            theta_dot = theta_dot
            psi = -psi
            psi_dot = -psi_dot
        
        R_BW = R_theta(theta) @ R_psi(psi) @ R_phi(phi) @ R_beta(beta)

        #convention for naming: V_W_body is velocity of the body expressed in the wing frame (W)
        V_W_body = R_BW @ self._state[3:6]      #rotation @ body velocity
        omega_W_wing = np.array([[0], [theta_dot], [0]]) + R_theta(theta) @ np.array([[psi_dot], [0], [0]]) + R_theta(theta) @ R_psi(psi) @ np.array([[0], [0], [phi_dot]])
        #print(self.side, omega_W_wing)
        F_W_trans = np.array([[0],[0],[0]])
        F_W_rot = np.array([[0],[0],[0]])
        M_W_trans = np.array([[0],[0],[0]])
        M_W_rot = np.array([[0],[0],[0]])

        for i in range(self.num_elements):
            #aero calculations
            delta_r = self.wing_length/self.num_elements
            if self.side == 'left':                                            #side when looking at fwmav, (side along positive y axis)
                r_W_i = np.array([[0], [(i+1)*delta_r], [0]])                 #distance from wing pivot to i-th blade element 
                #r_W_i = np.array([[0], [(i+1)*delta_r - delta_r/2], [0]])
            elif self.side == 'right':
                delta_r = self.wing_length/self.num_elements
                r_W_i = np.array([[0], [-(i+1)*delta_r], [0]])
                #r_W_i = np.array([[0], [-(i+1)*delta_r + delta_r/2], [-self.wing_chord/4]])
                #r_W_i = np.array([[0], [-(i+1)*delta_r + delta_r/2], [0]])
            chord_i = self.chord_f(r_W_i.item(1))
            
            V_W_inflow = V_W_body + np.cross(omega_W_wing.T, r_W_i.T).T
            #print(self.side, V_W_inflow)
            V_W_i = np.array([[1, 0, 0],[0, 0, 0],[0, 0, 1]]) @ V_W_inflow
            chat_W_i = np.array([[0], [0], [-1]])                   #c hat in Wing frame for element i   
            if (V_W_i.T @ chat_W_i) == 0:
                alpha = np.pi/2
            else:
                alpha = np.arctan((np.sqrt(np.cross(V_W_i.T, chat_W_i.T) @ np.cross(V_W_i.T, chat_W_i.T).T))/(V_W_i.T @ chat_W_i))
                #print(self.side, np.degrees(alpha))
            CL, CD, CM, CR = self.get_aero_coeff(alpha)
            #print(CM)
            #calculate frames
            i_W = np.array([[1], [0], [0]])
            j_W = np.array([[0], [1], [0]])
            s_W_i = np.array([[0], [np.dot(r_W_i.T, j_W).item(0)], [chord_i/4]])
            lhat_W_i = (np.dot(V_W_i.T, i_W))/(np.sqrt(np.dot(V_W_i.T, i_W) * np.dot(V_W_i.T, i_W))) * np.array([[0, 0, 1], [0, 1, 0], [-1, 0, 0]]) @ (V_W_i)/(np.sqrt(V_W_i.T @ V_W_i))    
            dhat_W_i = np.array([[-1, 0, 0], [0, 0, 0], [0, 0, -1]]) @ (V_W_i)/(np.sqrt(V_W_i.T @ V_W_i))
            if self.side == 'right':
                j_W = np.array([[0], [1], [0]])
                s_W_i = np.array([[0], [np.dot(r_W_i.T, j_W).item(0)], [chord_i/4]])
                lhat_W_i = (np.dot(V_W_i.T, i_W))/(np.sqrt(np.dot(V_W_i.T, i_W) * np.dot(V_W_i.T, i_W))) * np.array([[0, 0, 1], [0, -1, 0], [-1, 0, 0]]) @ (V_W_i)/(np.sqrt(V_W_i.T @ V_W_i))    
                dhat_W_i = np.array([[-1, 0, 0], [0, 0, 0], [0, 0, -1]]) @ (V_W_i)/(np.sqrt(V_W_i.T @ V_W_i))

            #print(self.side, V_W_i.T @ V_W_i)
            #if self._state.item(10) > 0:
                #print(self.side, V_W_i.T @ V_W_i)
            #print(lhat_W_i, dhat_W_i)
            #Force Calculations
            #print(chord_f(r_W_i.item(1)))
            F_W_trans_i = (CL* 0.5 * sim.rho * (V_W_i.T @ V_W_i) * chord_i * delta_r) * lhat_W_i + (CD * 0.5 * sim.rho * (V_W_i.T @ V_W_i) * chord_i * delta_r) * dhat_W_i
            F_W_rot_i = (CR * sim.rho * np.dot(omega_W_wing.T, j_W) * np.sqrt(V_W_i.T @ V_W_i) * chord_i**2 * delta_r) * i_W
            F_W_trans = F_W_trans + F_W_trans_i
            F_W_rot = F_W_rot + F_W_rot_i
            #print(F_W_trans)
            #print(np.cross(r_W_i.T, F_W_trans_i.T))
            #print(self.side, np.sign(self._state.item(10)))
            M_W_trans = M_W_trans + np.sign(self._state.item(10))*(CM * 0.5 * sim.rho * V_W_i.T @ V_W_i * chord_i**2 * delta_r) * j_W + np.cross(r_W_i.T, F_W_trans_i.T).T
            #print(np.cross(r_W_i.T, F_W_trans_i.T).T)
            #print(np.cross(s_W_i.T, F_W_rot_i.T).T)
            #M_W_trans = M_W_trans + (CM * 0.5 * sim.rho * V_W_i.T @ V_W_i * self.wing_chord**2 * delta_r) * j_W
            #M_W_trans = M_W_trans + 0*np.cross(r_W_i.T, F_W_trans_i.T).T
            M_W_rot = M_W_rot + np.cross(s_W_i.T, F_W_rot_i.T).T
            #print(self.side, (s_W_i.T, F_W_rot_i.T))

        F_W = F_W_trans + F_W_rot
        M_W = M_W_trans + M_W_rot
        #print(F_W)
        F_B = R_BW.T @ F_W
        M_B = R_BW.T @ M_W
        #print(R_BW.T)
        
        return np.array([[F_B.item(0), F_B.item(1), F_B.item(2), M_B.item(0), M_B.item(1), M_B.item(2)]]).T

            

def R_theta(theta):
    c_theta = np.cos(theta)
    s_theta = np.sin(theta)
    return np.array([[c_theta, 0, -s_theta],
                        [0, 1, 0],
                        [s_theta, 0, c_theta]])
def R_psi(psi):
    c_psi = np.cos(psi)
    s_psi = np.sin(psi)
    return np.array([[1, 0, 0],
                      [0, c_psi, -s_psi],
                      [0, s_psi, c_psi]])
def R_phi(phi):
    c_phi = np.cos(phi)
    s_phi = np.sin(phi)
    return np.array([[c_phi, s_phi, 0],
                      [-s_phi, c_phi, 0],
                      [0, 0, 1]])
def R_beta(beta):
    c_beta = np.cos(beta)
    s_beta = np.sin(beta)
    return np.array([[c_beta, 0, s_beta],
                        [0, 1, 0],
                        [-s_beta, 0, c_beta]])



