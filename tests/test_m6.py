from src.entities.models import TaskTarget
from src.use_cases.calc import send_to_calc_flow_rate, send_to_calc_p_bhp
import pandas as pd

from test_calc import eq_by_mape

_path = 'test_m6'
input_data = {
    "base_prop": {
        "calc_type": 'calc_type',
        "segment_count": 20,
        "border_type_x": "n",
        "border_type_y": "n",
        "border_type_z_up": "c",
        "border_type_z_down": "n",
    },
    "target": {
        "time_step": 0.01,
        "number_of_steps": 0,
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
            "horizontal_wellbore_length": 624,
            "horizontal_wellbore_perf_ratio": 1,
            "horizontal_perf_count": 5,
            "permeability": 100000000
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
            "xe": 1000,
            "ye": 1000,
            "lc_ratio": 0.5,
            "wc_rectangle_ratio": 0.5
        }
    }
}


# @pytest.mark.skip
def test_calc_flow_rate_m6():
    for calc_type in [
        'optimal',  # 42 sec ok
        'segmentation',  # 41 min ok
    ]:
        expected_df = pd.read_csv(f"{_path}/expected_m6_{calc_type}_flow_rate.csv")

        input_data['base_prop']['calc_type'] = calc_type
        input_data['target']['number_of_steps'] = len(expected_df)
        actual = (
            send_to_calc_flow_rate(TaskTarget(**input_data), fake_t=expected_df["t"])
        )["result"]

        actual_df = pd.DataFrame(actual)
        actual_df.to_csv(f"{_path}/actual_m6_{calc_type}_flow_rate.csv", index=False)
        eq_by_mape(actual_df, expected_df, skip={'t'})


# @pytest.mark.skip
def test_calc_p_bhp_m6():
    for calc_type in [
        'optimal',  # 42 sec ok
        'segmentation',  # 40 min ok
    ]:
        expected_df = pd.read_csv(f'{_path}/expected_m6_{calc_type}_pressure.csv')

        input_data['base_prop']['calc_type'] = calc_type
        input_data['target']['number_of_steps'] = len(expected_df)
        actual = (
            send_to_calc_p_bhp(TaskTarget(**input_data), fake_t=expected_df["t"])
        )["result"]

        actual_df = pd.DataFrame(actual)
        actual_df.to_csv(f'{_path}/actual_m6_{calc_type}_pressure.csv', index=False)
        eq_by_mape(actual_df, expected_df, skip={'t'})
