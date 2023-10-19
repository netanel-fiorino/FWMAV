import numpy as np
import parameters.FWMAV_parameters as FWMAV

class MsgDelta:
    def __init__(self,
                 rudder=0.0,
                 throttle=0.5,
                 kappa=0.0,
                 flap=0,
                 tail_angle = FWMAV.tail_angle):
        self.rudder = rudder  # rudder command
        self.throttle = throttle  # throttle command
        self.kappa = kappa          #wing flap angle
        self.flap = flap            #flapping or glidding (1 or 0)
        self.tail_angle = tail_angle

    def to_array(self):
        return np.array([[self.rudder],
                         [self.throttle],
                         [self.kappa],
                         [self.flap],
                         [self.tail_angle]])

    def from_array(self, u):
        self.rudder = u.item(0)
        self.throttle = u.item(1)
        self.kappa = u.item(2)
        self.flap = u.item(3)
        self.tail_angle = u.item(4)

    def print(self):
        print('rudder=', self.rudder,
              'throttle=', self.throttle,
              'kappa=', self.kappa,
              'flap=', self.flap,
              'tail angle=', self.tail_angle)