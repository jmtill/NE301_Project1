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
                           data_dir= "/Users/jonathantill/PycharmProjects/NE301_Project1/data",
                           save_name=f"timestep_{timestep}"
                           )

decay_solver.run()


print(nuclide_a.decay_data["300.0"])
print(nuclide_b.decay_data["300.0"])

tar = 300

key, value = get_nuc_number(nuclide_a, target=tar)

print(key, value)







"""
Write a computer program that will solve the radioactive decay equations numerically using the Euler
method presented in class. The program should be able to solve for an arbitrary number of coupled
nuclides.
"""

"""
Part 1: 
The purpose of Part 1 is to compare your program results to an analytic solution so you can show that
the program is working properly. This is known as “code verification”.
Solve the numerical decay equations for the following problem:
 2 nuclides, labeled “A” and “B”.
 The half-lives of A and B are 1.2 min, and 2.0 min, respectively.
 The initial quantities of A and B are zero.
 Nuclide A is being produced at a rate of 20,000 atom/min.
 Solve for the time interval 0 to 20 min.
 The “quantity of interest” (QOI) is the number of nuclides A and B at 5 min.


Plot the nuclide concentrations of A and B as a function of time as calculated by your program.
Include dashed lines of the equilibrium values
Run the problem with at least six timestep sizes and plot the absolute error between the
calculated QOI and analytic QOI. The plot should show the absolute error vs. the time step size.
Choose at least one time-step size so the absolute error is less than 1 atom.
e) What relationship between the error and time step size is observed? (generate a curve fit for
error vs. time step size.)
f) What time step (t) is necessary to obtain a maximum absolute error of the QOI to less than 1
atoms?
"""