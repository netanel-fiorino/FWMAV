import numpy as np
import matplotlib.pyplot as plt
plt.close('all')

class WingDynamics:
    def __init__(self, v0, freq, side):
        self._state = np.array([[0],    # (x-pos) [0]
                               [0],     # (y-pos) [1]
                               [0],     # (z-pos) [2]
                               [v0.item(0)],     # (u) [3]
                               [v0.item(1)],     # (v) [4]
                               [v0.item(2)],     # (w) [5]
                               [np.radians(0)],     # (beta) angle to stroke plane [6]
                               [np.radians(-90)],     # (phi) causes flapping motion [7]
                               [0],     # (psi) deviation angle [8]
                               [0],      # (theta) [9]
                               [2*np.pi*freq],      # (phi_dot) [10]
                               [0],      # (psi_dot) [11]
                               [0]])     # (theta_dot) [12]
        self.wing_length = .05
        self.wing_chord = .01
        self.freq = freq
        self.u_inf = np.sqrt(self._state.item(3)**2 + self._state.item(4)**2 + self._state.item(5)**2)
        self.J = self.u_inf/(2*np.radians(180)*self.wing_length*self.freq)
        self.side = side
        self.num_elements = 5

    
    def update(self, ts):
        self._state[7][0] = self._state.item(7) + self._state.item(10)*ts
        
    def updateOmega(self, omega):
        self._state[10][0] = omega
        
    
    def get_aero_coeff(self, alpha):
        self.u_inf = np.sqrt(self._state.item(3)**2 + self._state.item(4)**2 + self._state.item(5)**2)
        self.J = self.u_inf/(2*np.radians(180)*self.wing_length*self.freq)
        xhat_0 = 0.75
        rhat_2 = np.sqrt(1/3)                   #non dimensional second moment of the wing
        rhat_M = np.sqrt(1/3)
        K_PL = -2.109*((self.J + rhat_2)**-0.606) + 4.136
        K_VL = 2.659*((self.J + rhat_2)**-0.666) + -0.344
        K_PD = -0.182*((self.J + rhat_2)**-2.414) + 1.370
        K_VD = 0.765*((self.J + rhat_2)**-1.497) + 2.078
        K_PM = 0.803*((self.J + rhat_M)**-0.972) + -0.363
        K_VM = -0.242*((self.J + rhat_M)**-1.354) + -0.554


        CL = K_PL*np.sin(alpha)*(np.cos(alpha)**2)+K_VL*(np.sin(alpha)**2)*np.cos(alpha)
        CD = K_PD*(np.sin(alpha)**2)*np.cos(alpha) + K_VD*(np.sin(alpha)**3)
        CM = K_PM*(np.sin(alpha)**2)*np.cos(alpha) + K_VM*(np.sin(alpha)**2)
        CR = np.pi*(0.75-xhat_0)
        #print(CL, CD, CM, CR)

        return CL, CD, CM, CR

    def force_calc(self):
        beta = self._state.item(6)
        phi = self._state.item(7)
        psi = self._state.item(8)
        theta = self._state.item(9)
        phi_dot =self._state.item(10)
        psi_dot = self._state.item(11)
        theta_dot = self._state.item(12)


        R_BW = R_theta(theta) @ R_psi(psi) @ R_phi(phi) @ R_beta(beta)

        #convention for naming: V_W_body is velocity of the body expressed in the wing frame (W)
        V_W_body = R_BW @ self._state[3:6]      #rotation @ body velocity
        omega_W_wing = np.array([[0], [theta_dot], [0]]) + R_theta(theta) @ np.array([[psi_dot], [0], [0]]) + R_theta(theta) @ R_psi(psi) @ np.array([[0], [0], [phi_dot]])
        
        F_W_trans = np.array([[0],[0],[0]])
        F_W_rot = np.array([[0],[0],[0]])
        M_W_trans = np.array([[0],[0],[0]])
        M_W_rot = np.array([[0],[0],[0]])

        for i in range(self.num_elements): #range(self.num_elements):
            #aero calculations
            delta_r = self.wing_length/self.num_elements
            if self.side == 'left':                                            #side when looking at fwmav, (side along positive y axis)
                r_W_i = np.array([[0], [(i+1)*delta_r], [0]])                 #distance from wing pivot to i-th blade element
            elif self.side == 'right':
                r_W_i = np.array([[0], [-(i+1)*delta_r], [0]])
            V_W_inflow = V_W_body + np.cross(omega_W_wing.T, r_W_i.T).T
            #print(V_W_inflow)
            V_W_i = np.array([[1, 0, 0],[0, 0, 0],[0, 0, 1]]) @ V_W_inflow
            chat_W_i = np.array([[0], [0], [-1]])                   #c hat in Wing frame for element i   
            #print(V_W_i.T @ chat_W_i)
            if (V_W_i.T @ chat_W_i) == 0:
                alpha = np.pi/2
            else:
                #print(np.cross(V_W_i.T, chat_W_i.T).shape, np.cross(V_W_i.T, chat_W_i.T).shape)
                alpha = np.arctan((np.sqrt(np.cross(V_W_i.T, chat_W_i.T) @ np.cross(V_W_i.T, chat_W_i.T).T))/(V_W_i.T @ chat_W_i))
                #print(np.degrees(alpha))
            CL, CD, CM, CR = self.get_aero_coeff(alpha)
            #calculate frames
            i_W = np.array([[1], [0], [0]])
            j_W = np.array([[0], [1], [0]])

            
            lhat_W_i = (np.dot(V_W_i.T, i_W))/(np.sqrt(np.dot(V_W_i.T, i_W) * np.dot(V_W_i.T, i_W))) * np.array([[0, 0, 1], [0, 1, 0], [-1, 0, 0]]) @ (V_W_i)/(np.sqrt(V_W_i.T @ V_W_i))    
            dhat_W_i = np.array([[-1, 0, 0], [0, 0, 0], [0, 0, -1]]) @ (V_W_i)/(np.sqrt(V_W_i.T @ V_W_i))
            #print(lhat_W_i, dhat_W_i)
            #Force Calculations
            F_W_trans_i = (CL* 0.5 * rho * (V_W_i.T @ V_W_i) * self.wing_chord * delta_r) * lhat_W_i + (CD * 0.5 * rho * (V_W_i.T @ V_W_i) * self.wing_chord * delta_r) * dhat_W_i
            F_W_rot_i = (CR * rho * np.dot( omega_W_wing.T, j_W) * np.sqrt(V_W_i.T @ V_W_i) * self.wing_chord**2 * delta_r) * i_W
            F_W_trans = F_W_trans + F_W_trans_i
            #F_W_rot = F_W_rot + F_W_rot_i
            #print(F_W_trans)
            #print(np.cross(r_W_i.T, F_W_trans_i.T))
            M_W_trans = M_W_trans + (CM * 0.5 * rho * V_W_i.T @ V_W_i * self.wing_chord**2 * delta_r) * j_W #+ np.cross(r_W_i.T, F_W_trans_i.T)
            #M_W_rot = M_W_rot + np.cross(s_W_i, F_W_rot_i)

        F_W = F_W_trans + F_W_rot
        M_W = M_W_trans + M_W_rot
        #print(F_W)
        F_B = R_BW.T @ F_W
        M_B = R_BW.T @ M_W
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



t_period = .05     #sec
freq = 1/t_period   #Hz
t_step = t_period/100
v_free_stream = 5*np.array([[1], [0], [0]])   #m/s
rho = 1.2682
amplitude = np.radians(180)

#print(R_theta(0) @ R_psi(0) @ R_phi(np.radians(-0)) @ R_beta(np.radians(90)))
wing_left = WingDynamics(v_free_stream, freq, 'left')
wing_right = WingDynamics(v_free_stream, freq, 'right')
time = 0
omega = 2*np.pi*freq*np.cos(2*np.pi*freq*time)
wing_left.updateOmega(omega)
wing_right.updateOmega(omega)
timeArr = []
forceArr_left = []
forceArr_right = []
momentArr_left = []
momentArr_right = []
angleArr_left = []
i = 0
while time < t_period/2:
    tmp_left = wing_left.force_calc()
    wing_left.update(t_step)
    tmp_right = wing_left.force_calc()
    wing_right.update(t_step)

    timeArr.append(time)

    forceArr_left.append(tmp_left[:3])
    momentArr_left.append(tmp_left[3:6])

    forceArr_right.append(tmp_left[:3])
    momentArr_right.append(tmp_left[3:6])

    angleArr_left.append(wing_left._state.item(7))
    time = time + t_step
    omega = 2*np.pi*freq *np.cos(2*np.pi*freq*time)
    wing_left.updateOmega(omega)
    wing_right.updateOmega(omega)
    i+=1
    


while time < t_period:
    tmp_left = wing_left.force_calc()
    wing_left.update(t_step)
    tmp_right = wing_left.force_calc()
    wing_right.update(t_step)

    timeArr.append(time)

    forceArr_left.append(tmp_left[:3])
    momentArr_left.append(tmp_left[3:6])

    forceArr_right.append(tmp_left[:3])
    momentArr_right.append(tmp_left[3:6])

    angleArr_left.append(wing_left._state.item(7))
    time = time + t_step
    omega = 2*np.pi*freq *np.cos(2*np.pi*freq*time)
    wing_left.updateOmega(omega)
    wing_right.updateOmega(omega)
    i+=1
    
    
forceArr_left = np.array(forceArr_left) 
forceArr_right = np.array(forceArr_right)
momentArr_left = np.array(momentArr_left) 
momentArr_right = np.array(momentArr_right)

plt.figure()
plt.plot(timeArr, forceArr_left[:, 0], label = 'fx from left')
plt.plot(timeArr, forceArr_left[:, 1], label = 'fy from left')
plt.plot(timeArr, forceArr_left[:, 2], label = 'fz from left')

plt.plot(timeArr, forceArr_right[:, 0], label = 'fx from right')
plt.plot(timeArr, forceArr_right[:, 1], label = 'fy from right')
plt.plot(timeArr, forceArr_right[:, 2], label = 'fz from right')
plt.legend()
plt.show()

plt.figure()
plt.plot(timeArr, momentArr_left[:, 0], label = 'Mx left')
plt.plot(timeArr, momentArr_left[:, 1], label = 'My left')
plt.plot(timeArr, momentArr_left[:, 2], label = 'Mz left')

plt.plot(timeArr, momentArr_right[:, 0], label = 'Mx right')
plt.plot(timeArr, momentArr_right[:, 1], label = 'My right')
plt.plot(timeArr, momentArr_right[:, 2], label = 'Mz right')
plt.legend()
plt.show()

plt.figure()
plt.plot(timeArr, np.degrees(angleArr_left), label = 'phi')
plt.legend()
plt.show()