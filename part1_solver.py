from decay_solver import DecaySolver
from nuclide_class import Nuclide


def get_nuc_number(nuclide: Nuclide, target: float = 300.0):
    # chooses the key whose numeric value is closest to target
    d = nuclide.decay_data
    nearest_key = min(d, key=lambda k: abs(float(k) - target))
    return nearest_key, d[nearest_key]

nuclide_a = Nuclide(name="A",
                    half_life=1.2,
                    half_life_unit="minutes",
                    parents = None,
                    external_prod_rate=20000,
                    prod_rate_unit_time="minutes",
                    time=0.0,
                    n_0 = 0.0
                    )

nuclide_b = Nuclide(name="B",
                    parents=["A"],
                    half_life=2.0,
                    half_life_unit="minutes",
                    time=0.0,
                    n_0 = 0
                    )


nuclides:dict[str, Nuclide] = {str(nuclide_a.name): nuclide_a,
                               str(nuclide_b.name): nuclide_b,
                               }
# [0.05,0.005,0.001, 0.0005, 0.0001,0.00005]
timestep = 0.05
decay_solver = DecaySolver(nuclides=nuclides,
                           timestep=timestep,
                           timestep_unit="seconds",
                           duration=20.0,
                           duration_unit="minutes",
                           plot_data=False,
                           save_data=True,
                           data_dir= "./data",
                           save_name="AB_nuclides"
                           )

# run solver
decay_solver.run()


print(nuclide_a.decay_data["300.0"])
print(nuclide_b.decay_data["300.0"])



