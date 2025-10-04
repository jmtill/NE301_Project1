from error_analysis_utils import *


parent_path = "./data"
ts_1 = "timestep_0.5_data"
ts_2 = "timestep_0.05_data"
ts_3 = "timestep_0.01_data"
ts_4 = "timestep_0.005_data"
ts_5 = "timestep_0.001_data"
ts_6 = "timestep_0.0005_data"

timesteps = [0.45,0.40,0.35, 0.30, 0.25,0.20, 0.15, 0.10]


df_1 = pd.read_csv(f"{parent_path}/{ts_1}.csv")
df_2 = pd.read_csv(f"{parent_path}/{ts_2}.csv")
df_3 = pd.read_csv(f"{parent_path}/{ts_3}.csv")
df_4 = pd.read_csv(f"{parent_path}/{ts_4}.csv")
df_5 = pd.read_csv(f"{parent_path}/{ts_5}.csv")
df_6 = pd.read_csv(f"{parent_path}/{ts_6}.csv")


frames = [df_1, df_2, df_3, df_4, df_5, df_6]

nuclide_a = Nuclide(name="A",
                    half_life=1.2,
                    half_life_unit="minutes",
                    parents=None,
                    external_prod_rate=20000,
                    prod_rate_unit_time="minutes",
                    time=0.0,
                    n_0=0.0
                    )

nuclide_b = Nuclide(name="B",
                    parents=["A"],
                    half_life=2.0,
                    half_life_unit="minutes",
                    time=0.0,
                    n_0=0
                    )


# function plots horizontal line at the numerical y value closest to analytical value. As such, only one df is used for this plot.
plot_eq_vals(df_5, nuclide_a, nuclide_b)

qoi_time = 5 # minutes

# note, warnings may appear in terminal but do not prevent code from running
plot_numer_error(time=qoi_time,nuclideA=nuclide_a, nuclideB=nuclide_b, frames=frames)