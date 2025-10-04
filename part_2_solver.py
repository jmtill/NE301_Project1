import pandas as pd

from decay_solver import DecaySolver
from nuclide_class import Nuclide
from pathlib import Path


def get_nuc_number(nuclide: Nuclide, target: float = 300.0):
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

decay_solver = DecaySolver(nuclides=nuclides,
                           timestep=0.001,
                           timestep_unit="seconds",
                           duration=10000.0,
                           duration_unit="seconds",
                           plot_data=True,
                           save_data=True,
                           data_dir= "/Users/jonathantill/PycharmProjects/NE301_Project1/data",
                           save_name="part_two_test"
                           )

decay_solver.run()
