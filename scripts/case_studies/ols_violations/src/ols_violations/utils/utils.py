from dataclasses import dataclass

# @dataclass
class SimulationParams:
    # n_sims: int
    # n_samples: int
    # beta: float
    # alpha: float
    # rho: float
    # sigma : float

    def __init__(self, **kwargs):

        for (k, v) in kwargs.items():
            self.__setattr__(k, v)

    def __iter__(self):
        for k, v in self.__dict__.items():
            yield k, v
