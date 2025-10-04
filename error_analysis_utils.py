from nuclide_class import Nuclide
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from pandas import DataFrame
from decay_solver import DecaySolver
import copy


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


def NA_analytical(nuclide: Nuclide, time: float):
    """
    Computes the nuclide amount at a specific time using an analytical solution to the
    radioactive decay and production equation.

    The function calculates the remaining quantity of a radioactive nuclide over time
    considering both natural decay and an external production rate. The analytical formula
    takes into account the initial nuclide amount, decay constant, and external production
    rate to return the quantity of the nuclide at the given time.

    :param nuclide: Object representing the nuclide, which includes attributes such as
        initial quantity, decay constant, and external production rate.
    :type nuclide: Nuclide
    :param time: Time at which the nuclide amount is calculated.
    :type time: float

    :return: Calculated quantity of the nuclide at the specified time.
    :rtype: float
    """
    n0 = nuclide.n_0
    lam = nuclide.decay_const
    prod = nuclide.external_prod_rate

    nt = n0*np.exp(-lam * time) + (prod/lam)*(1.0 - np.exp(-lam * time))
    return float(nt)

def NB_analytical(nuclide: Nuclide, parent: Nuclide, time: float):
    """
    Calculate the analytical solution for the concentration of a nuclide
    accounting for its production and decay as well as the decay of its parent.

    This function computes the concentration of a nuclide over time,
    taking into consideration the decay constant of the nuclide, the decay
    constant of its parent, and the external production rate of the parent nuclide.

    :param nuclide: The daughter nuclide whose concentration is being calculated.
    :type nuclide: Nuclide
    :param parent: The parent nuclide contributing to the production of the daughter nuclide.
    :type parent: Nuclide
    :param time: The time over which the calculation is performed.
    :type time: float
    :return: The analytically calculated concentration of the daughter nuclide at the given time.
    :rtype: float
    """
    parent_prod = parent.external_prod_rate
    lam_b = nuclide.decay_const
    lam_a = parent.decay_const

    first_term = (1-np.exp(-lam_b*time))/lam_b

    n = np.exp(-lam_a*time) - np.exp(-lam_b*time)
    d = lam_b - lam_a

    n_t = parent_prod * (first_term - n/d)
    return float(n_t)


def plot_eq_vals(dataset: DataFrame,
                  nuclideA: Nuclide,
                  nuclideB: Nuclide):
    """
    Plots equilibrium values of concentrations for two nuclides along with concentration data over time.

    This function takes a dataset containing time-series nuclide concentration data, and two nuclides,
    calculates their equilibrium concentrations, and visualizes this data. The equilibrium points of the
    given nuclides are overlaid as horizontal dashed lines for comparison.

    :param dataset: A pandas DataFrame containing nuclide concentration data with time as a variable.
    :param nuclideA: An instance of the Nuclide class representing nuclide A.
    :param nuclideB: An instance of the Nuclide class representing nuclide B.
    :return: None. The function generates and displays a plot.
    """
    nA_eq = nuclideA.external_prod_rate/nuclideA.decay_const
    nB_eq = nuclideA.external_prod_rate/nuclideB.decay_const
    df = dataset.copy()


    nA_eq_idx = abs(df["A"]- nA_eq).idxmin()
    nB_eq_idx = abs(df["B"]- nB_eq).idxmin()

    nA_eq_point = df.iloc[nA_eq_idx]["A"]
    nB_eq_point = df.iloc[nB_eq_idx]["B"]


    fig, ax = plt.subplots()
    sns.set_theme(context="paper")


    df_long = df.melt(id_vars=["time (minutes)"],
                      value_vars=["A","B"],
                      var_name="nuclide",
                      value_name="nuclide concentration")

    ax = sns.lineplot(data=df_long, x=f"time (minutes)", y="nuclide concentration", hue="nuclide", ax=ax,
                 errorbar=None)


    ax.axhline(nA_eq_point, linestyle="--", color="red", label=f"Nuclide A EQ Value: {round(nA_eq_point,5)}")
    ax.axhline(nB_eq_point, linestyle="--", color="black", label=f"Nuclide B EQ Value: {round(nB_eq_point, 5)}")

    ax.legend()
    ax.set(xlabel=f"Time (minutes)", ylabel="Nuclide Concentration (count/unit volume)")

    sns.move_legend(ax, "best")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()


def plot_numer_error(time: float, frames: list, nuclideA: Nuclide, nuclideB: Nuclide):
    """
    Calculate numerical error for nuclide concentrations and generate a bar plot.

    This function computes the numerical errors in the calculated concentrations
    of two nuclides (nuclideA and nuclideB) over a given time interval. Errors are
    calculated by comparing numerical results to analytical solutions. The errors
    are then visualized using a bar plot to show the error dependence on time steps.

    :param time: Time of interest in minutes.
    :type time: float
    :param frames: List of data frames containing time step information and nuclide
        concentrations.
    :type frames: list
    :param nuclideA: First nuclide for which the numerical error will be calculated.
    :type nuclideA: Nuclide
    :param nuclideB: Second nuclide for which the numerical error will be calculated.
    :type nuclideB: Nuclide
    :return: None
    """
    n_A_an: float = NA_analytical(nuclideA, time*60)
    n_B_an: float = NB_analytical(nuclideB, nuclideA, time*60)


    steps: list = []
    nA_errors: list = []
    nB_errors: list = []

    for frame in frames:
        idx = abs(frame['time (minutes)'] - time).idxmin()
        n_A_num = frame.iloc[idx][1]
        n_B_num = frame.iloc[idx][2]


        timestep = abs(frame.iloc[1][0] - frame.iloc[0][0])*60 # converts to seconds timestep
        timestep =  round(timestep, 5)
        steps.append(timestep)


        a_error = abs(n_A_an - n_A_num)
        print("Nuclide A absolute error: ", a_error)

        b_error = abs(n_B_an - n_B_num)

        print("Nuclide B absolute error: ", b_error, "\n")

        nA_errors.append(a_error)
        nB_errors.append(b_error)

    errors: dict = {"time_steps": steps,
                    "nuclide_a": nA_errors,
                    "nuclide_b": nB_errors}

    error_df = pd.DataFrame(errors)

    error_long = error_df.melt(id_vars = "time_steps", value_vars = ["nuclide_a", "nuclide_b"], var_name = 'nuclide', value_name = "nuclide_error")
    fig, ax = plt.subplots()

    sns.set_theme(context="paper")


    sns.barplot(data=error_long, ax=ax, x="time_steps", y="nuclide_error", hue="nuclide", width=0.5)
    ax.set_xlabel("Time Step (seconds)")
    ax.set_ylabel("Absolute Error (nuclides/unit volume)")
    ax.set_xlabel("Time Steps (seconds)")
    ax.set_ylabel("Absolute Error (nuclides/unit volume)")



    plt.tight_layout()
    plt.show()
    return


def check_convergence(tolerance, nuclides: dict[str, Nuclide], timestep_i, time_stop):
    dt = timestep_i
    error_inf: float = 0

    while True:
        error_inf: float = 0

        chain1 = copy.deepcopy(nuclides)
        chain2 = copy.deepcopy(nuclides)

        dt_solver = DecaySolver(nuclides=chain1,
                                   timestep=dt,
                                   timestep_unit="seconds",
                                   duration=time_stop,
                                   duration_unit="seconds",
                                   plot_data=False,
                                   save_data=False,
                                   data_dir="./data",
                                   save_name=f"p2_timestep_{timestep_i}"
                                )

        dt_solver.run()

        half_dt_solver = DecaySolver(nuclides=chain2,
                                   timestep=(dt/2),
                                   timestep_unit="seconds",
                                   duration=time_stop,
                                   duration_unit="seconds",
                                   plot_data=False,
                                   save_data=False,
                                   data_dir="./data",
                                   save_name=f"p2_timestep_{(timestep_i)/2}")

        half_dt_solver.run()

        for name in (chain1.keys() & chain2.keys()):
            nuclide1 = chain1[name]
            nuclide2 = chain2[name]

            # Get last recorded populations
            k1 = next(reversed(nuclide1.decay_data))
            k2 = next(reversed(nuclide2.decay_data))

            n_1 = nuclide1.decay_data[k1]
            n_2 = nuclide2.decay_data[k2]

            # Update max absolute difference
            error_inf = max(error_inf, abs(n_1 - n_2))

        if error_inf < tolerance:
            print(f"Converged at dt={dt / 2}, error={error_inf}")
            for nuclide in chain2.values():
                _, val = get_nuc_number(nuclide, target=3600.0)
                print(f"nuclide {nuclide.name}: {val}")
            print("\n")

            break
        else:
            print(f"Did not converge at dt={dt / 2}, error={error_inf}")
            dt = dt / 2


