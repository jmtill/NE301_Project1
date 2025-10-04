import pandas as pd

from decay_solver import DecaySolver
from nuclide_class import Nuclide
from pathlib import Path
import copy

from error_analysis_utils import check_convergence


def get_nuc_number(nuclide: Nuclide, target: float = 300.0):
    """
    Determines the nuclide key and its corresponding value from the decay data
    dictionary in the `nuclide` object that is closest to the given `target` value.
    If no specific `target` is provided, the default value of 300.0 is used. The method
    iterates over the decay data keys and computes the absolute difference between
    each key and the target, returning the key with the smallest difference.

    :param nuclide: Object used to access its `decay_data` dictionary for comparison.
    :type nuclide: Nuclide
    :param target: Numeric value to find the closest corresponding key in nuclide's
        decay data. Defaults to 300.0.
    :type target: float
    :return: A tuple containing the key closest to the `target` and its associated
        value from the nuclide's decay data.
    :rtype: tuple
    """

    # chooses the key whose numeric value is closest to target
    d = nuclide.decay_data
    nearest_key = min(d, key=lambda k: abs(float(k) - target))
    return nearest_key, d[nearest_key]


Se89 = Nuclide(name="Se89",
                    half_life=0.41,
                    half_life_unit="seconds",
                    parents = None,
                    external_prod_rate=1000000,
                    prod_rate_unit_time="seconds",
                    time=0.0,
                    n_0 = 0.0
                    )

Br89 = Nuclide(name="Br89",
                    parents=["Se89"],
                    half_life=4.41,
                    half_life_unit="seconds",
                    time=0.0,
                    n_0 = 0.0
                    )

Kr89 = Nuclide(name="Kr89",
                    half_life=189.0,
                    half_life_unit="seconds",
                    parents = ["Br89"],
                    time=0.0,
                    n_0 = 0.0
                    )

Rb89 = Nuclide(name="Rb89",
                    parents=["Kr89"],
                    half_life=909.0,
                    half_life_unit="seconds",
                    time=0.0,
                    n_0 = 0.0
                    )

Sr89 = Nuclide(name="Sr89",
                    half_life=4363200.0,
                    half_life_unit="seconds",
                    parents = ["Rb89"],
                    time=0.0,
                    n_0 = 0.0
                    )

nuclides:dict[str, Nuclide] = {str(Se89.name): Se89,
                               str(Br89.name): Br89,
                               str(Kr89.name): Kr89,
                               str(Rb89.name): Rb89,
                               str(Sr89.name): Sr89
                               }

timestep_i = 1.0
tolerance = 1.0


check_convergence(tolerance=tolerance, nuclides = nuclides, timestep_i=timestep_i, time_stop=10000)




# if you would like to check the decay solver alone, uncomment the code below
# for demo/speed purposes timestep is set to 0.5. Timestep is smaller for figures/data in the report
# decay_solver = DecaySolver(nuclides=nuclides,
#                            timestep=0.5,
#                            timestep_unit="seconds",
#                            duration=10000.0,
#                            duration_unit="seconds",
#                            plot_data=True,
#                            save_data=True,
#                            data_dir= "./data",
#                            save_name=f"p2_timestep_0.5"
#                            )
#
# decay_solver.run()
#
# parent_path = "./data"
#
# ts_1 = f"timestep_0.5_data"
#
# print(f"nuclides QOI at timestep: 0.5 seconds:")
#
#
# for nuclide in nuclides.values():
#     _, val = get_nuc_number(nuclide, target=3600.0)
#     print(f"nuclide {nuclide.name}: {val}")



