m_airframe=.212e-3
#%m_transmission=.141e-3
m_transmission=.212e-3
m_motor=1.217e-3
m_wings=.083e-3*2
m_markers=.013e-3*4
#%m_legs=0.011e-3*4
m_legs=0
#%m_nose=0.02e-3
m_nose=0
#%m_launchertab=0.036e-3
m_launchertab=0.05e-3

m_board=.3e-3


m_robot=m_airframe+m_transmission+m_motor+m_wings+m_markers+m_nose+m_legs+m_launchertab+m_board;               #kg,

#m_bat=.415e-3
m_bat=1.2e-3
#m_board=.203e-3
m_wires=.07e-3
#m_boardmount=.075e-3
m_boardmount=0
m_batterymount=.1e-3

#m_battery=m_bat+m_board+m_wires+m_boardmount;         #kg, includes control board
m_battery=m_bat+m_wires+m_batterymount

m_tail=A_tail*tail_density_per_area+0.07e-3            #kg (estimate) including markers
m_rod=d*rod_density_per_length

m=m_robot + m_battery + m_tail + m_rod + additional_weight
print(m)

width_robot=13e-3           
height_robot=8e-3           
length_robot=28e-3          
length_wing=6e-2
width_wing=3e-2