from nuclide_class import Nuclide
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from pandas import DataFrame
from pathlib import Path



class DecaySolver:
    """
    Manages the simulation of radioactive nuclide decay using the Euler method.

    This class calculates the decay of multiple nuclides, including parent decay
    relationships and external production rates, using explicit Euler methods.
    It supports configuration of simulation parameters such as time steps, durations,
    and optional saving and plotting of resulting data.

    The solver facilitates adding nuclides dynamically and computes their concentrations
    over time.

    :ivar nuclides: Dictionary to store nuclide objects by their names.
    :type nuclides: dict[str, Nuclide]
    :ivar timestep: Converted timestep value in seconds for the decay simulation.
    :type timestep: float
    :ivar timestep_unit: Unit of the provided timestep, such as 'seconds', 'minutes', etc.
    :type timestep_unit: str
    :ivar duration: Total duration of the simulation, converted to seconds.
    :type duration: float
    :ivar duration_unit: Unit of the provided simulation duration, such as 'seconds',
        'minutes', etc.
    :type duration_unit: str
    :ivar iterations: Number of Euler method iterations to compute over the specified duration.
    :type iterations: int
    :ivar times_list: Precomputed time values (as strings) for each timestep in seconds.
    :type times_list: list[str]
    :ivar plot_data: Indicates whether plots of the nuclide decay data should be generated.
    :type plot_data: bool
    :ivar save_data: Indicates whether the decay data should be saved to a file.
    :type save_data: bool
    :ivar save_name: Filename for saving the decay data, if `save_data` is True.
    :type save_name: str or None
    :ivar save_path: Directory path where data will be saved, or `None` if not applicable.
    :type save_path: Path or None
    """
    def __init__(self,
                 nuclides: dict[str, Nuclide] | None = None,
                 timestep: float | None = None,
                 timestep_unit: str | None = None,
                 duration: float | None = None,
                 duration_unit: str | None = None,
                 plot_data: bool = True,
                 save_data: bool = False,
                 data_dir: str | Path | None = None,
                 save_name: str | None = None
                 ) -> None:

        # placeholder till data is compiled
        self.decay_data = None
        self.plot_data = plot_data
        self.save_data = save_data
        self.save_name = save_name

        if nuclides is None:
            self.nuclides = dict()
        else:
            self.nuclides = nuclides

        # Conversion factors from non‑second units to seconds for timestep and duration
        self.time_conversions: dict = {
            "minutes": 60,
            "hours": 3600,
            "days": 86400,
            "weeks": 604800,
            "months": 2.628*(10**6),
            "years": 3.154*(10**7),
        }

        # Validate and Convert timestep
        if timestep_unit is None:
            raise ValueError("Must provide timestep units")
        elif timestep_unit not in  ["seconds","minutes", "hours", "days", "weeks", "months", "years"]:
            raise ValueError("Timestep unit is not valid. Valid option examples: "
                             "seconds, minutes, hours, days, weeks, months, years")

        self.timestep_unit = timestep_unit

        if timestep is None:
            raise ValueError("Must provide timestep")
        else:
            self.timestep = self._convert_time(unit=timestep_unit, time=timestep)

        if duration_unit is None:
            raise ValueError("Must provide duration unit")
        if duration_unit not in ["seconds","minutes", "hours", "days", "weeks", "months", "years"]:
            raise ValueError("duration unit is not valid. Valid option examples:"
                             "seconds, minutes, hours, days, weeks, months, years")
        else:
            self.duration_unit = duration_unit

        # Validate and Convert Duration
        if duration is None:
            raise ValueError("Must provide duration")
        elif float(duration) <= 0.0:
            raise ValueError("Duration must be positive")

        self.duration = self._convert_time(unit=duration_unit, time=duration)

        # Number of Euler iterations (integer number of steps). Ensure at least one step.
        self.iterations = int(self.duration/self.timestep)

        # Precompute list of time strings for each sample including t=0.  These are used
        # as dictionary keys into each nuclide's decay_data.  Use seconds as the base and
        # convert to strings for exact matching.
        self.times_list: list[str] = [f"{i * self.timestep}" for i in range(self.iterations + 1)]

        if plot_data or save_data:
            if data_dir is None:
                self.save_path: Path = Path.cwd()
            else:
                self.save_path = Path(data_dir)
        else:
            self.save_path= None

    def add_nuclide(self, nuclide: str,
                    parent_nuclides: list[str],
                    half_life: float,
                    n_0: float = 0,
                    external_prod_rate: float = 0,
                    half_life_unit: str = "seconds",
                    prod_rate_unit: str = "seconds")-> None:
        """
        Adds a new nuclide entry to the collection of nuclides.

        The method creates an instance of a nuclide with the specified properties
        and adds it to an internal collection using the nuclide's name as a key.

        :param nuclide: The name of the nuclide to be added.
        :type nuclide: str
        :param parent_nuclides: A list of strings representing the nuclide's parent nuclides.
        :type parent_nuclides: list[str]
        :param half_life: The half-life value of the nuclide in the specified units.
        :type half_life: float
        :param n_0: The initial amount of nuclide. Defaults to 0.
        :type n_0: float
        :param external_prod_rate: The external production rate of the nuclide. Defaults to 0.
        :type external_prod_rate: float
        :param half_life_unit: The unit of measurement for the half-life. Defaults to "seconds".
        :type half_life_unit: str
        :param prod_rate_unit: The unit of time for production rate. Defaults to "seconds".
        :type prod_rate_unit: str
        :return: None
        :rtype: None
        """
        nuclide_entry: dict[str:Nuclide] = Nuclide(name= nuclide,
                                                   parents= parent_nuclides,
                                                   external_prod_rate= float(external_prod_rate),
                                                   prod_rate_unit_time= prod_rate_unit,
                                                   half_life=half_life,
                                                   half_life_unit= half_life_unit,
                                                   n_0=float(n_0)
                                                   )

        self.nuclides[nuclide] = nuclide_entry

        return None

    def _convert_time(self, unit:str, time: float):
        """
        Converts a given time value into seconds based on the provided unit.

        This method takes a time value and a unit of time and converts the
        value into seconds. The available units depend on the entries in
        `time_conversions`. If the unit is "seconds", the input value is
        returned as is.

        :param unit: The unit of the provided time value. It must be a key
                     in the `time_conversions` mapping or "seconds".
        :type unit: str
        :param time: The time value to convert to seconds.
        :type time: float
        :return: The converted time value in seconds.
        :rtype: float
        """
        if unit == "seconds":
            return float(time)
        else:
            time = time*self.time_conversions[unit]
            return float(time)

    def calc_euler_decay(self) -> None:
        """
        Compute decay of radionuclides over time using an explicit Euler method.

        This method implements an explicit Euler method for computing the decay and production
        of radionuclides over a specified number of iterations and time steps. For every iteration
        beyond the initial step, the production and decay of each nuclide are computed based on the
        contributions from its parent nuclides, external production rate, and its decay constant.
        The computed values are then stored for each time step.

        The updates for all nuclides at each time step are performed in isolation, ensuring that the
        results of each step are based on the state at the beginning of the step. Each nuclide's
        time-series data is updated in a dictionary indexed by the time in seconds.

        :rtype: None
        """
        # For each time step beyond t=0, compute the increment
        for idx in range(1, self.iterations + 1):
            # Build a temporary mapping of updated amounts so that updates
            # are based on the values from the previous step (explicit Euler)
            next_values: dict[str, float] = {}
            for name, nuclide in self.nuclides.items():
                parents: list[Nuclide] = [self.nuclides[p] for p in nuclide.parents if p in self.nuclides]


                parent_sum: float = sum(parent.decay_const * parent.n_t for parent in parents)
                dN_dt: float = parent_sum + nuclide.external_prod_rate - nuclide.decay_const * nuclide.n_t

                # Euler update
                next_n_t: float = nuclide.n_t + dN_dt * self.timestep
                next_values[name] = next_n_t

                # Append to the decay_data dictionary.  Key is the time in seconds as string.
                t_val = self.times_list[idx]
                nuclide.decay_data[t_val] = next_n_t

            # Now update the current values for the next step
            for name, value in next_values.items():
                self.nuclides[name].n_t = value

    def _make_dataframe(self):
        """
        Creates and returns a preallocated DataFrame initialized with necessary columns
        and appropriate data types for simulation iterations.

        The method constructs a DataFrame with rows corresponding to the number of
        iterations plus one (to include t=0) and columns for each nuclide and time
        values. The time column is populated with float values representing time steps,
        and converted to the specified duration unit if different from seconds.

        :return: A preallocated pandas DataFrame containing the initialized data for
            iterations and time steps.
        :rtype: pandas.DataFrame
        """
        # rows = N+1 to include t=0
        n_rows = self.iterations + 1
        nuclide_names = [nuclide.name for nuclide in self.nuclides.values()]
        cols = [f"time ({self.duration_unit})"] + list(nuclide_names)

        # Preallocate with float dtype; initialize to 0.0
        df = pd.DataFrame(0.0, index=np.arange(n_rows), columns=cols)

        # Fill the time column once
        df[f"time ({self.duration_unit})"] = np.arange(n_rows, dtype=float) * self.timestep

        if self.duration_unit != "seconds":
            df[f"time ({self.duration_unit})"] = df[f"time ({self.duration_unit})"] / self.time_conversions[self.duration_unit]
        return df

    def compile_data(self):
        """
        Compiles decay data from multiple nuclides into a single DataFrame.

        This method ensures that the decay data for all nuclides aligns with
        the predefined solver times and consolidates the data into a DataFrame.
        Inconsistent time keys or mismatches in data alignment may lead to errors.

        :raises ValueError: If time keys for any nuclide do not match the solver times.
        :return: None
        """
        decay_df: DataFrame = self._make_dataframe()

        # For each nuclide, ensure times match and populate the DataFrame
        for nuclide in self.nuclides.values():
            nuclide_times = list(nuclide.decay_data.keys())


            if nuclide_times != self.times_list:
                raise ValueError(
                    f"Time keys for nuclide '{nuclide.name}' do not match solver times."
                )
            decay_df[nuclide.name] = list(nuclide.decay_data.values())
        self.decay_data = decay_df

    def _make_decay_plot(self, ax):
        """
        Generate a decay plot visualization for nuclide concentration over time.

        This method creates a line plot that represents the variation of nucleide
        concentrations as a function of time based on the provided decay data. The
        data is first transformed into a long format compatible with seaborn's
        lineplot function. A theme suitable for publication context is applied,
        and the plot is customized with axis labels and legends.

        :param ax: A matplotlib.axes.Axes instance where the decay plot will be
            drawn. This should be provided as an existing matplotlib axes to
            allow precise placement and integration into larger figures.
        :return: The provided matplotlib.axes.Axes instance, modified to contain
            the decay plot visualization.
        """
        sns.set_theme(context="paper")

        df = self.decay_data.copy()
        df_long = df.melt(id_vars=[f"time ({self.duration_unit})"],
                               value_vars=list(self.nuclides.keys()),
                               var_name="nuclide",
                               value_name="nuclide concentration")

        sns.lineplot(data=df_long, x=f"time ({self.duration_unit})", y="nuclide concentration", hue="nuclide",ax=ax, errorbar=None)
        ax.set(xlabel= f"Time ({self.duration_unit})", ylabel="Nuclide Concentration (count/unit volume)")
        sns.move_legend(ax, "best")
        ax.grid(True, alpha=0.3)

        return ax

    def _make_plots(self):
        """
        Creates and displays plots for visualization purposes. It adjusts the layout
        for better presentation and ensures the plot is displayed once created.

        :raises RuntimeError: If there's an issue with displaying the plot.
        :returns: None
        """
        fig, ax = plt.subplots(figsize=(6, 4))
        ax = self._make_decay_plot(ax)

        plt.tight_layout()
        plt.show()

    def save_nuclide_data(self):
        """
        Saves nuclide data including decay data and associated plot paths to specified
        storage location. This method ensures that the directory exists before saving the
        data file and plot file.

        :raises ValueError: If ``save_path`` is not defined.
        :raises ValueError: If ``save_name`` is not specified.
        :return: None
        """
        if self.save_path is None:
            raise ValueError("save_dir is None; cannot save data without a valid directory")
        if self.save_name is None:
            raise ValueError("save_name must be provided to save data")

        # Make filenames
        data_filename: str = f"{self.save_name}_data.csv"

        self.save_path.mkdir(parents=True, exist_ok=True)

        data_path = self.save_path / data_filename

        self.decay_data.to_csv(path_or_buf=data_path, index=False, float_format="%.15g")

    def run(self):
        self.calc_euler_decay()
        self.compile_data()

        if self.plot_data:
            self._make_plots()
        if self.save_data:
            print("saving data...")
            self.save_nuclide_data()

        # make sure to compile data before running plots
        pass
