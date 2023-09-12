import numpy as np


class MsgDelta:
    def __init__(self,
                 rudder=0.0,
                 throttle=0.5,
                 kappa=0.0,
                 flap=0):
        self.rudder = rudder  # rudder command
        self.throttle = throttle  # throttle command
        self.kappa = kappa          #wing flap angle
        self.flap = flap            #flappinf or glidding (1 or 0)

    def to_array(self):
        return np.array([[self.rudder],
                         [self.throttle],
                         [self.kappa],
                         [self.flap]])

    def from_array(self, u):
        self.rudder = u.item(0)
        self.throttle = u.item(1)
        self.kappa = u.item(2)
        self.flap = u.item(3)

    def print(self):
        print('rudder=', self.rudder,
              'throttle=', self.throttle,
              'kappa=', self.kappa)