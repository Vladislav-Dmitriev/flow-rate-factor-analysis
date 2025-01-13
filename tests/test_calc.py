from src.entities.models import TaskTarget
from src.use_cases.calc import send_to_calc_flow_rate, send_to_calc_p_bhp
import pandas as pd
from pandas.api.types import is_numeric_dtype
from sklearn.metrics import mean_absolute_percentage_error


def eq_by_mape(df1: pd.DataFrame, df2: pd.DataFrame, skip={}, limit: float = 0.02):
    # assert list(df1.columns) == list(df2.columns)
    for column in df1.columns:
        if column in skip:
            continue
        if is_numeric_dtype(df1[column]):
            assert (
                    mean_absolute_percentage_error(df1[column], df2[column]) < limit
            ), f"""{mean_absolute_percentage_error(df1[column], df2[column])} < {limit}: {column}"""


#@pytest.mark.skip
def test_calc_flow_rate():
    input_data = {
        "base_prop": {
            "calc_type": "optimal",
            "segment_count": 20,
            "border_type_x": "n",
            "border_type_y": "n",
            "border_type_z_up": "n",
            "border_type_z_down": "n",
        },
        "target": {
            "time_step": 0.01,
            "number_of_steps": 231,
            "cumulative_work_time": 100,
            "target_values": {
                "q_liq": 50,
                "p_bhp": 150
            },
        },
        "unit": {
            "skin": 0,
            "h_eff": 10,
            "vertical_offset": 0.5,
            "perfres_ratio": 1,
            "wellbore_prop": {
                "wellbore_type": "vertical",
                "wellbore_r": 0.1,
                "horizontal_wellbore_length": None,
                "horizontal_wellbore_perf_ratio": None,
                "horizontal_perf_count": None,
                "permeability": None,
            },
            "layer_prop": {
                "permeability": 10,
                "kv_kh_ratio": 0.1,
                "compressibility": 0.00012,
                "porosity": 0.2,
                "p_res_init": 250,
                "viscosity_oil": 1,
                "b_oil": 1.2,
                "res_model_type": "Homogeneous",
                "f_compressibility": 0.0001,
                "f_porosity": 0.001,
                "lambda": 0.00001,
                "internal_r": None,
                "kmu_in_out_ratio": 2,
                "kmuphict_in_out_ratio": 2,
                "grp_flag": False,
                "grp_prop": None,
                "mgrp_flag": False,
                "xe": 500,
                "ye": 500,
                "lc_ratio": 0.5,
                "wc_rectangle_ratio": 0.5,
            },
        },
    }
    expected_df = pd.read_csv("expected_calc_flow_rate.csv")
    actual = (
        send_to_calc_flow_rate(TaskTarget(**input_data), fake_t=expected_df["t"])
    )["result"]
    actual_df = pd.DataFrame.from_dict(actual)
    actual_df.to_csv("actual_calc_flow_rate.csv", index=False)
    eq_by_mape(actual_df, expected_df)


#@pytest.mark.skip
def test_calc_p_bhp():
    input_data = {
        "base_prop": {
            "calc_type": "optimal",
            "segment_count": 20,
            "border_type_x": "n",
            "border_type_y": "n",
            "border_type_z_up": "n",
            "border_type_z_down": "n",
        },
        "target": {
            "time_step": 0.01,
            "number_of_steps": 262,
            "cumulative_work_time": 100,
            "target_values": {
                "q_liq": 50,
                "p_bhp": 150
            },
        },
        "unit": {
            "skin": 0,
            "h_eff": 10,
            "vertical_offset": 0.5,
            "perfres_ratio": 1,
            "wellbore_prop": {
                "wellbore_type": "vertical",
                "wellbore_r": 0.1,
                "horizontal_wellbore_length": None,
                "horizontal_wellbore_perf_ratio": None,
                "horizontal_perf_count": None,
                "permeability": None,
            },
            "layer_prop": {
                "permeability": 10,
                "kv_kh_ratio": 0.1,
                "compressibility": 0.00012,
                "porosity": 0.2,
                "p_res_init": 250,
                "viscosity_oil": 1,
                "b_oil": 1.2,
                "res_model_type": "Homogeneous",
                "f_compressibility": 0.0001,
                "f_porosity": 0.001,
                "lambda": 0.00001,
                "internal_r": None,
                "kmu_in_out_ratio": 2,
                "kmuphict_in_out_ratio": 2,
                "grp_flag": False,
                "grp_prop": None,
                "mgrp_flag": False,
                "xe": 500,
                "ye": 500,
                "lc_ratio": 0.5,
                "wc_rectangle_ratio": 0.5,
            },
        },
    }
    expected_df = pd.read_csv("expected_calc_p_bhp.csv")
    actual = (
        send_to_calc_p_bhp(TaskTarget(**input_data), fake_t=expected_df["t"])
    )["result"]
    actual_df = pd.DataFrame.from_dict(actual)
    actual_df.to_csv("actual_calc_p_bhp.csv", index=False)
    eq_by_mape(actual_df, expected_df)
