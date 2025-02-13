import pandas as pd
import pytest

from src.entities.models import TaskTarget
from src.use_cases.calc import send_to_calc_flow_rate, send_to_calc_p_bhp
from test_calc import eq_by_mape

input_data = {
    "base_prop": {
        "calc_type": '',
        "segment_count": 20,
        "border_type_x": "n",
        "border_type_y": "n",
        "border_type_z_up": "n",
        "border_type_z_down": "n",
    },
    "target": {
        "cumulative_time": 1000,
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
            "wellbore_type": "horizontal",
            "wellbore_r": 0.1,
            "horizontal_wellbore_length": 624,
            "horizontal_wellbore_perf_ratio": 1,
            "horizontal_perf_count": 5,
            "permeability": 100000000
        },
        "layer_prop": {
            "permeability": 10,
            "water_cut": 10,
            "p_bubble": 100,
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
def test_calc_flow_rate_horizont():
    for calc_type in [
        'optimal',
        'desuperposition',
        'segmentation'
    ]:
        expected_df = pd.read_csv(f"horisont_csv/Expected_Flow_{calc_type}.csv")
        input_data['base_prop']['calc_type'] = calc_type
        input_data['target']['number_of_steps'] = len(expected_df)
        model = TaskTarget(**input_data)
        actual = (
            send_to_calc_flow_rate(model, fake_t=expected_df["t"])
        )["result"]
        actual_df = pd.DataFrame({key: [value] for key, value in actual.items()})
        actual_df.to_csv(f"horisont_csv/actual_calc_flow_rate_{calc_type}.csv", index=False)
        expected_df = pd.DataFrame({'t': [2400.0], 'JD': [1.029817263], 'flow_rate_result': [41.76220505]})
        eq_by_mape(actual_df, expected_df, skip={'t'})


@pytest.mark.skip
def test_calc_p_bhp_horizont():
    for calc_type in [
        'optimal',
        'desuperposition',
        'segmentation'
    ]:
        expected_df = pd.read_csv(f'horisont_csv/Expected_Pres_{calc_type}.csv')
        input_data['base_prop']['calc_type'] = calc_type
        input_data['target']['number_of_steps'] = len(expected_df)

        actual = (
            send_to_calc_p_bhp(TaskTarget(**input_data), fake_t=expected_df["t"])
        )["result"]
        actual_df = pd.DataFrame(actual)
        # actual_df.to_csv(f"horisont_csv/actual_calc_p_bhp_{calc_type}.csv", index=False)
        eq_by_mape(actual_df, expected_df, skip={'t'})
