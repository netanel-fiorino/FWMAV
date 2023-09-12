import numpy as np


class MsgDelta:
    def __init__(self,
                 elevator=0.0,
                 aileron=0.0,
                 rudder=0.0,
                 throttle=0.5,
                 kappa=0.0):
        self.elevator = elevator  # elevator command
        self.aileron = aileron  # aileron command
        self.rudder = rudder  # rudder command
        self.throttle = throttle  # throttle command
        self.kappa = kappa          #wing flap angle

    def to_array(self):
        return np.array([[self.elevator],
                         [self.aileron],
                         [self.rudder],
                         [self.throttle],
                         [self.kappa]])

    def from_array(self, u):
        self.elevator = u.item(0)
        self.aileron = u.item(1)
        self.rudder = u.item(2)
        self.throttle = u.item(3)
        self.kappa = u.item(4)

    def print(self):
        print('elevator=', self.elevator,
              'aileron=', self.aileron,
              'rudder=', self.rudder,
              'throttle=', self.throttle,
              'kappa=', self.kappa)