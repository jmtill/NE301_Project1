from nuclide_class import Nuclide
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from pandas import DataFrame


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
                    n_0 = 0.0
                    )
parent_path = "/Users/jonathantill/PycharmProjects/NE301_Project1/data"
ts_1 = "timestep_0.05_data"
ts_2 = "timestep_0.01_data"
ts_3 = "timestep_0.001_data"
ts_4 = "timestep_0.0001_data"
ts_5 = "timestep_0.005_data"
ts_6 = "timestep_0.0005_data"

df_1 = pd.read_csv(f"{parent_path}/{ts_1}.csv")
df_2 = pd.read_csv(f"{parent_path}/{ts_2}.csv")
df_3 = pd.read_csv(f"{parent_path}/{ts_3}.csv")
df_4 = pd.read_csv(f"{parent_path}/{ts_4}.csv")
df_5 = pd.read_csv(f"{parent_path}/{ts_5}.csv")
df_6 = pd.read_csv(f"{parent_path}/{ts_6}.csv")

frames = [df_1, df_2, df_3, df_4, df_5, df_6]



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
        b_error = abs(n_B_an - n_B_num)
        print(n_B_an, n_B_num, "\n")
        print(b_error)

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




plot_eq_vals(dataset=df_1, nuclideA=nuclide_a, nuclideB=nuclide_b)


plot_numer_error(time=5, frames=frames, nuclideA=nuclide_a, nuclideB=nuclide_b)




