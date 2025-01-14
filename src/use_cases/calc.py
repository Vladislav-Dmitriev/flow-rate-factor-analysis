import numpy as np
import traceback

from src.entities.models import TaskTarget
from src.use_cases.builder import Builder
from logger import logger


def send_to_calc_flow_rate(input_data: TaskTarget, fake_t=None) -> dict:
    """Отправка параметров для расчета дебита жидкости

    :param input_data: Входные данные для расчета параметров
    :type input_data: src.entities.models.TaskTarget

    :param fake_t: Тестовое накопленное время работы

    :return: Дебит жидкости (в сутки)
    :rtype: dict
    """
    res: dict = {"flow_rate_result": None, "t": None, "JD": None}
    input_data_dict = input_data.model_dump(by_alias=True)
    builder = Builder(**input_data_dict)
    p_bhp = input_data_dict["target"]["target_values"]["p_bhp"]
    try:
        if fake_t is None:
            t = builder.calc_t()
        else:
            t = np.array(fake_t)
        S = builder.calc_S(t)
        q_d = builder.calc_q_d(S)
        S_new = builder.calc_double_porosity(S)
        p_d = builder.calc_p_d(S_new)
        delta_p = builder.calc_delta_p(p_bhp)
        flow_rate = builder.calc_flow_rate(S, q_d, p_d, delta_p)
        res["flow_rate_result"] = flow_rate
        res["t"] = t.tolist()
        cumulative_prod = builder.calc_cumulative_prod(S, q_d, p_d, p_bhp)
        JD = builder.calc_JD_for_flow_rate(p_bhp, flow_rate, cumulative_prod)
        res["JD"] = JD
    except Exception as e:
        msg = f"Error calc_flow_rate route: {e}!"
        logger.error(msg)
    else:
        msg = "calc_flow_rate was successful!"
        logger.info(msg)

    response = {
        "result": res,
        "message": msg,
    }
    return response


def send_to_calc_p_bhp(input_data: TaskTarget, fake_t=None) -> dict:

    """
    Используем эту функцию для расчета JD
    """

    """Отправка параметров для расчета давления

    :param input_data: Входные данные для расчета параметров
    :type input_data: src.entities.models.TaskTarget

    :param fake_t: Тестовое накопленное время работы

    :return: Результат расчета параметров
    :rtype: dict
    """

    res: dict = {"t": {}, "JD": {}, "p_wf_result": {}}

    input_data_dict = input_data.model_dump(by_alias=True)
    builder = Builder(**input_data_dict)
    ql = input_data_dict["target"]["target_values"]["q_liq"]
    try:
        if fake_t is None:
            t = builder.calc_t()
        else:
            t = np.array(fake_t)
        S = builder.calc_S(t)
        q_d = builder.calc_q_d(S)
        S_new = builder.calc_double_porosity(S)
        p_d = builder.calc_p_d(S_new)
        p_wf, delta_p = builder.calc_p_bhp(S, q_d, p_d, ql)
        res["p_wf_result"] = p_wf
        res["t"] = t.tolist()
        q_storage = builder.calc_q_storage(S, q_d, p_d, ql)
        JD = builder.calc_JD_for_pressure(t, ql, q_storage, delta_p)
        res["JD"] = JD
    except Exception as e:
        msg = f"Error calc_p_bhp route: {e}!"
        logger.error(traceback.format_exc(chain=False))
    else:
        msg = "calc_p_bhp was successful!"
        logger.info(msg)

    response = {
        "result": res,
        "message": msg,
    }

    return response
