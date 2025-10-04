import numpy as np


class Nuclide:
    """
    Represents a radioactive nuclide with decay and production properties.

    This class models the behavior of a radioactive nuclide, including its
    decay over time and external production rate. It is designed to handle
    various time units for half-life and production rate, converting them
    to a consistent unit (seconds) internally. The class supports
    initialization of attributes such as name, half-life, parent nuclides,
    and external production rate, and includes key calculations for decay
    constantly and tracking nucleotide quantities by time.

    :ivar name: Name of the nuclide.
    :type name: str
    :ivar parents: List of parent nuclides.
    :type parents: list[str]
    :ivar n_t: Current number of nuclide atoms at a given time.
    :type n_t: float
    :ivar n_0: Initial number of nuclide atoms at the start.
    :type n_0: float
    :ivar time: Current time for the nuclide's state tracking.
    :type time: float
    :ivar prod_rate_unit: Unit of measure for production rate time.
    :type prod_rate_unit: str
    :ivar half_life_unit: Unit of measure for half-life.
    :type half_life_unit: str
    :ivar half_life: Half-life of the nuclide in seconds.
    :type half_life: float
    :ivar decay_const: Decay constant of the nuclide, calculated from
        the half-life.
    :type decay_const: float
    :ivar external_prod_rate: External production rate of the nuclide
        in atoms per second.
    :type external_prod_rate: float
    :ivar decay_data: Dictionary storing time points and the quantities
        of the nuclide at those times.
    :type decay_data: dict[str, float]
    """
    def __init__(self,
                 name: str,
                 half_life: float,
                 parents: list[str]|None,
                 prod_rate_unit_time: str = "seconds",
                 half_life_unit: str = "seconds",
                 n_0: float = 0.0,
                 external_prod_rate: float = 0.0,
                 time: float = 0.0,
                )->None:

        # Conversion factors from non‑second units to seconds
        self.second_conversion : dict = {
            "minutes": 60,
            "hours": 3600,
            "days": 86400,
            "weeks": 604800,
            "months": 2.628*(10**6),
            "years": 3.154*(10**7),
        }

        self.name: str = name
        if parents is None:
            self.parents = []
        else:
            self.parents: list[str] = parents

        # sets n_t = to n_0 at startup of class, gets overwritten when in DecaySolver
        self.n_t: float = n_0
        self.n_0: float = n_0
        self.time: float = time

        self.prod_rate_unit: str = prod_rate_unit_time
        self.half_life_unit: str = half_life_unit

        self.half_life = self._calc_half_life(half_life)
        self.decay_const = self._calc_decay_constant()

        if self.prod_rate_unit == "seconds":
            self.external_prod_rate = external_prod_rate
        else:
            if self.prod_rate_unit not in self.second_conversion:
                raise ValueError(
                    f"Unsupported production rate unit '{self.prod_rate_unit}'."
                )
            self.external_prod_rate = float(external_prod_rate/self.second_conversion[self.prod_rate_unit])

        self.decay_data: dict[str, float] = {str(0.0): self.n_0}

    def _calc_half_life(self, half_life: float) -> float:
        units_list: list = ["seconds","minutes", "hours",
                            "days", "weeks", "months", "years"]


        if self.half_life_unit not in units_list:
            raise ValueError(f"Half life unit {self.half_life_unit} not in {units_list}")

        if self.half_life_unit == "seconds":
            return half_life

        if self.half_life_unit not in self.second_conversion:
            raise ValueError(
                f"Unsupported half‑life unit '{self.half_life_unit}'."
            )
        return float(half_life) * self.second_conversion[self.half_life_unit]

    def _calc_decay_constant(self):
        return np.log(2) / self.half_life


