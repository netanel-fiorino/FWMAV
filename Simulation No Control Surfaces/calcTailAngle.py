
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
CL_body = -c_l


q_dynamic_tail = 0.5*FWMAV.rho*(self._Va_tail**2)*FWMAV.S_tail        # dynamic pressure --> force term (*S) for lift/drag
F_lift_tail = q_dynamic_tail*(CL_tail) #q_dynamic_tail*(CL_tail + (FWMAV.C_L_q*0.5*FWMAV.c/self._Va_tail*q))
F_drag_tail = q_dynamic_tail*(CD_tail) #q_dynamic_tail*(CD_tail + (FWMAV.C_D_q*0.5*FWMAV.c/self._Va_tail*q))

q_dynamic_body = 0.5*FWMAV.rho*(self._Va**2)*FWMAV.S_body        # dynamic pressure --> force term (*S) for lift/drag
F_lift_body = q_dynamic_body*(CL_body) #q_dynamic_tail*(CL_tail + (FWMAV.C_L_q*0.5*FWMAV.c/self._Va_tail*q))
F_drag_body = q_dynamic_body*(CD_body) #q_dynamic_tail*(CD_tail + (FWMAV.C_D_q*0.5*FWMAV.c/self._Va_tail*q))

q_dynamic_wing = 0.5*FWMAV.rho*(self._Va**2)*FWMAV.S_wing       # dynamic pressure --> force term (*S) for lift/drag
F_lift_wing = q_dynamic_wing*(CL_tail) #q_dynamic_tail*(CL_tail + (FWMAV.C_L_q*0.5*FWMAV.c/self._Va_tail*q))
F_drag_wing = q_dynamic_wing*(CD_tail) #q_dynamic_tail*(CD_tail + (FWMAV.C_D_q*0.5*FWMAV.c/self._Va_tail*q))


force_tail =  np.array([[np.cos(self._alpha_tail), -np.sin(self._alpha_tail)], [0,0], [np.sin(self._alpha_tail), np.cos(self._alpha_tail)]]) @ np.array([-F_drag_tail, -F_lift_tail])
force_body = np.array([[np.cos(self._alpha), -np.sin(self._alpha)], [0,0], [np.sin(self._alpha), np.cos(self._alpha)]]) @ np.array([-F_drag_body, -F_lift_body])
force_wing = np.array([[np.cos(self._alpha), -np.sin(self._alpha)], [0,0], [np.sin(self._alpha), np.cos(self._alpha)]]) @ np.array([-F_drag_wing, -F_lift_wing])
fx = force_wing.item(0) + f_g.item(0) + force_tail.item(0) + force_body.item(0)
fy = force_wing.item(1) + f_g.item(1) + force_tail.item(1) + force_body.item(1) 
fz = force_wing.item(2) + f_g.item(2) + force_tail.item(2) + force_body.item(2)

M = np.cross(FWMAV.r_tail, force_tail.T) + np.cross(FWMAV.r_body, force_body.T) + np.cross(FWMAV.r_wing, force_wing.T)
My = 0#M.item(1)
Mx = 0
Mz = 0