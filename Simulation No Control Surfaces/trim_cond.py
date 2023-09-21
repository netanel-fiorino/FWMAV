import sys
sys.path.append('..')
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from message_types.msg_state import MsgState
import parameters.FWMAV_parameters as FWMAV
from tools.rotations import Quaternion2Euler, Quaternion2Rotation, Euler2Rotation, skew

def calculate_force(tail_angle, alpha):
    Ri_b = Euler2Rotation(0, alpha, 0)
    f_g = Ri_b @ np.array([
                            [0],
                            [0],
                            [FWMAV.mass*FWMAV.gravity],
                            ])

    alpha_tail = alpha - tail_angle
    #print(np.degrees(alpha_tail))
    Va = FWMAV.Va0
    Va_tail = FWMAV.Va0
    
    a0=7.159
    a1=-3.741
    b1=-0.08485
    a2=-0.7499
    b2=0.02863
    w=1.99
    c_d=a0+a1*np.cos(alpha_tail*w)+b1*np.sin(alpha_tail*w)+a2*np.cos(2*alpha_tail*w)+b2*np.sin(2*alpha_tail*w)
    CD_tail = c_d/3.8

    a0_2=-0.2495
    a1_2=-0.02425
    b1_2=5.896
    a2_2=0.1975
    b2_2=1.429
    a3_2=0.01897
    b3_2=0.4558
    w_2=2.011
    c_l=a0_2+a1_2*np.cos(alpha_tail*w_2)+b1_2*np.sin(alpha_tail*w_2)+a2_2*np.cos(2*alpha_tail*w_2)+b2_2*np.sin(2*alpha_tail*w_2)+a3_2*np.cos(3*alpha_tail*w_2)+b3_2*np.sin(3*alpha_tail*w_2)
    CL_tail = c_l/5.02

    a0=6.411
    a1=0.5402
    b1=-0.8709
    a2=-4.047
    b2=-0.1551
    w=1.654
    c_d=a0+a1*np.cos(alpha*w)+b1*np.sin(alpha*w)+a2*np.cos(2*alpha*w)+b2*np.sin(2*alpha*w)
    c_l=-1.045-0.7912*np.cos(1.835*alpha)+4.43*np.sin(1.835*alpha)
    CD_body = c_d
    CL_body = c_l


    q_dynamic_tail = 0.5*FWMAV.rho*(Va_tail**2)*FWMAV.S_tail        # dynamic pressure --> force term (*S) for lift/drag
    F_lift_tail = q_dynamic_tail*(CL_tail) #q_dynamic_tail*(CL_tail + (FWMAV.C_L_q*0.5*FWMAV.c/self._Va_tail*q))
    F_drag_tail = q_dynamic_tail*(CD_tail) #q_dynamic_tail*(CD_tail + (FWMAV.C_D_q*0.5*FWMAV.c/self._Va_tail*q))

    q_dynamic_body = 0.5*FWMAV.rho*(Va**2)*FWMAV.S_body        # dynamic pressure --> force term (*S) for lift/drag
    F_lift_body = q_dynamic_body*(CL_body) #q_dynamic_tail*(CL_tail + (FWMAV.C_L_q*0.5*FWMAV.c/self._Va_tail*q))
    F_drag_body = q_dynamic_body*(CD_body) #q_dynamic_tail*(CD_tail + (FWMAV.C_D_q*0.5*FWMAV.c/self._Va_tail*q))

    q_dynamic_wing = 0.5*FWMAV.rho*(Va**2)*FWMAV.S_wing       # dynamic pressure --> force term (*S) for lift/drag
    F_lift_wing = q_dynamic_wing*(CL_tail) #q_dynamic_tail*(CL_tail + (FWMAV.C_L_q*0.5*FWMAV.c/self._Va_tail*q))
    F_drag_wing = q_dynamic_wing*(CD_tail) #q_dynamic_tail*(CD_tail + (FWMAV.C_D_q*0.5*FWMAV.c/self._Va_tail*q))

    # compute longitudinal forces in body frame
    #fx = -F_drag*np.cos(self._alpha) + F_lift*np.sin(self._alpha) + f_g.item(0) + thrust_prop
    #fz = -F_drag*np.sin(self._alpha) - F_lift*np.cos(self._alpha) + f_g.item(2)

    force_tail =  np.array([[np.cos(alpha_tail), -np.sin(alpha_tail)], [0,0], [np.sin(alpha_tail), np.cos(alpha_tail)]]) @ np.array([-F_drag_tail, -F_lift_tail])
    force_body = np.array([[np.cos(alpha), -np.sin(alpha)], [0,0], [np.sin(alpha), np.cos(alpha)]]) @ np.array([-F_drag_body, -F_lift_body]) #np.array([0,0,0]) 
    force_wing = np.array([[np.cos(alpha), -np.sin(alpha)], [0,0], [np.sin(alpha), np.cos(alpha)]]) @ np.array([-F_drag_wing, -F_lift_wing])
    fx = force_wing.item(0) + f_g.item(0) + force_tail.item(0) + force_body.item(0)
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
    return fx, fz, My
print(calculate_force(np.radians(9.611), np.radians(11.1523)))
tail_a = np.linspace(0, np.pi/12, 1000)
alph = np.linspace(0, np.pi/12, 1000)
# fxArr = []
# fzArr = []
# myArr = []
# alphArr = []
# for i in range(len(alph)):
#     fxArr = []
#     fzArr = []
#     myArr = []
#     for j in range(len(tail_a)):
#         fx_temp, fz_temp, my_temp = calculate_force(tail_a[j], alph[i])
#         fxArr.append(fx_temp)
#         fzArr.append(fz_temp)
#         myArr.append(my_temp)
#     alphArr.append([fxArr, fzArr, myArr])


def align_yaxis(ax1, v1, ax2, v2):
    """adjust ax2 ylimit so that v2 in ax2 is aligned to v1 in ax1"""
    _, y1 = ax1.transData.transform((0, v1))
    _, y2 = ax2.transData.transform((0, v2))
    adjust_yaxis(ax2,(y1-y2)/2,v2)
    adjust_yaxis(ax1,(y2-y1)/2,v1)

def adjust_yaxis(ax,ydif,v):
    """shift axis ax by ydiff, maintaining point v at the same location"""
    inv = ax.transData.inverted()
    _, dy = inv.transform((0, 0)) - inv.transform((0, ydif))
    miny, maxy = ax.get_ylim()
    miny, maxy = miny - v, maxy - v
    if -miny>maxy or (-miny==maxy and dy > 0):
        nminy = miny
        nmaxy = miny*(maxy+dy)/(miny+dy)
    else:
        nmaxy = maxy
        nminy = maxy*(miny+dy)/(maxy+dy)
    ax.set_ylim(nminy+v, nmaxy+v)

# fig, ax1 = plt.subplots(projection='3d')
# ax1.set_xlabel('tail angle (rad)')
# ax1.set_ylabel('Force (n)')
# ax1.plot(alphArr[], label = "fx")
# ax1.plot(tail_a, fxArr, label = "fx")
# #ax1.plot(tail_a, fxArr, label = "fx")
# #ax1.plot(tail_a, fzArr, label = "fz")

# ax2 = ax1.twinx()
# ax2.set_ylabel("moment (nm)")
# ax2.plot(tail_a, myArr, 'green', label = "my")
# fig.tight_layout()
# fig.legend()
# align_yaxis(ax1, 0, ax2, 0)
# plt.grid()
# plt.show()

fxArr = np.zeros((len(alph), len(tail_a)))
fzArr = np.zeros((len(alph), len(tail_a)))
myArr = np.zeros((len(alph), len(tail_a)))

# Populate the arrays with calculated values
for i in range(len(alph)):
    for j in range(len(tail_a)):
        fxArr[i, j], fzArr[i, j], myArr[i, j] = calculate_force(tail_a[j], alph[i])

# Create a 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Create a meshgrid for the x and y values
X, Y = np.meshgrid(tail_a, alph)

# Normalize the data to fit within a similar range
max_force = max(np.max(fxArr), np.max(fzArr))
max_moment = np.max(myArr)
fxArr /= max_force
fzArr /= max_force
myArr /= max_moment

# Plot the normalized forces and moments using color coding
# cmap = plt.get_cmap('viridis')
# ax.plot_surface(X, Y, fxArr, cmap=cmap, label='Fx', alpha=0.7)
cmap = plt.get_cmap('plasma')
ax.plot_surface(X, Y, fzArr, cmap=cmap, label='Fz', alpha=0.7)
cmap = plt.get_cmap('inferno')
ax.plot_surface(X, Y, myArr, cmap=cmap, label='My', alpha=0.7)

# Label the axes
ax.set_xlabel('Tail Angle (tail_a)')
ax.set_ylabel('Angle of Attack (alph)')
ax.set_zlabel('Normalized Value')

# Show the plot
plt.show()

zero_my_alpha = []
zero_my_tail_a = []

# Populate the lists with values that meet the conditions
for alpha in alph:
    for tail_angle in tail_a:
        fx, fz, my = calculate_force(tail_angle, alpha)
        if np.abs(fz) <= 0.0001 and np.abs(my) <= 0.0001:
            zero_my_alpha.append(alpha)
            zero_my_tail_a.append(tail_angle)
            print("FOUND")

# Create a scatter plot for points where fz and my are both zero
plt.scatter(zero_my_tail_a, zero_my_alpha, marker='o', color='red', label='Fz = My = 0')

# Label the axes
plt.xlabel('Tail Angle (tail_a)')
plt.ylabel('Angle of Attack (alph)')

# Add a legend
plt.legend()

# Show the plot
plt.grid()
plt.show()

for i in range(len(zero_my_alpha)):
    fx, fz, my = calculate_force(zero_my_tail_a[i], zero_my_alpha[i])
    print(fx, fz, my)
    print(np.degrees(zero_my_tail_a[i]), np.degrees(zero_my_alpha[i]))
    print(50*'-')