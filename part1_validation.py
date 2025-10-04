from error_analysis_utils import *



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