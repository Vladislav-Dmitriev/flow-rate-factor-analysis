import copy

import numpy as np
import pandas as pd
from src.entities.models import TaskTarget
from src.use_cases.builder import Builder


class FactorAnalysis:

    def __init__(self, data_fact: dict, data_plan: dict, steps: int):
        # Разницы параметров
        self.k_delta = data_fact["unit"]["layer_prop"]["permeability"] - data_plan["unit"]["layer_prop"]["permeability"]
        self.h_delta = data_fact["unit"]["h_eff"] - data_plan["unit"]["h_eff"]
        self.Pwf_delta = data_fact["target"]["target_values"]["p_bhp"] - data_plan["target"]["target_values"]["p_bhp"]
        self.fw_delta = data_fact["unit"]["layer_prop"]["water_cut"] - data_plan["unit"]["layer_prop"]["water_cut"]
        self.Pr_delta = data_fact["unit"]["layer_prop"]["p_res_init"] - data_plan["unit"]["layer_prop"]["p_res_init"]
        self.ct_delta = (data_fact["unit"]["layer_prop"]["compressibility"] -
                         data_plan["unit"]["layer_prop"]["compressibility"])
        self.wf_delta = (data_fact["unit"]["layer_prop"]["grp_prop"]["wellbore_wf"] -
                         data_plan["unit"]["layer_prop"]["grp_prop"]["wellbore_wf"])
        self.kf_delta = (data_fact["unit"]["layer_prop"]["grp_prop"]["kf"] -
                         data_plan["unit"]["layer_prop"]["grp_prop"]["kf"])
        # Кол-во шагов
        self.steps = steps
        # Шаги по параметрам
        self.k_step = data_plan["unit"]["layer_prop"]["permeability"] / self.steps
        self.h_step = data_plan["unit"]["h_eff"] / self.steps
        self.Pwf_step = data_plan["target"]["target_values"]["p_bhp"] / self.steps
        self.fw_step = data_plan["unit"]["layer_prop"]["water_cut"] / self.steps
        self.Pr_step = data_plan["unit"]["layer_prop"]["p_res_init"] / self.steps
        self.ct_step = data_plan["unit"]["layer_prop"]["compressibility"] / self.steps
        self.wf_step = data_plan["unit"]["layer_prop"]["grp_prop"]["wellbore_wf"] / self.steps
        self.kf_step = data_plan["unit"]["layer_prop"]["grp_prop"]["kf"] / self.steps
        # Первые шаги для расчета
        self.k1 = data_plan["unit"]["layer_prop"]["permeability"] + self.k_step
        self.k2 = data_plan["unit"]["layer_prop"]["permeability"] - self.k_step
        self.h1 = data_plan["unit"]["h_eff"] + self.h_step
        self.h2 = data_plan["unit"]["h_eff"] - self.h_step
        self.P_wf1 = data_plan["target"]["target_values"]["p_bhp"] + self.Pwf_step
        self.P_wf2 = data_plan["target"]["target_values"]["p_bhp"] - self.Pwf_step
        self.f_w1 = data_plan["unit"]["layer_prop"]["water_cut"] + self.fw_step
        self.f_w2 = data_plan["unit"]["layer_prop"]["water_cut"] - self.fw_step
        self.Pr1 = data_plan["unit"]["layer_prop"]["p_res_init"] + self.Pr_step
        self.Pr2 = data_plan["unit"]["layer_prop"]["p_res_init"] - self.Pr_step
        self.ct1 = data_plan["unit"]["layer_prop"]["compressibility"] + self.ct_step
        self.ct2 = data_plan["unit"]["layer_prop"]["compressibility"] - self.ct_step
        self.wf1 = data_plan["unit"]["layer_prop"]["grp_prop"]["wellbore_wf"] + self.wf_step
        self.wf2 = data_plan["unit"]["layer_prop"]["grp_prop"]["wellbore_wf"] - self.wf_step
        self.kf1 = data_plan["unit"]["layer_prop"]["grp_prop"]["kf"] + self.kf_step
        self.kf2 = data_plan["unit"]["layer_prop"]["grp_prop"]["kf"] - self.kf_step
        # Входные данные
        self.data_plan = data_plan
        self.data_fact = data_fact
        # Параметры для факторного анализа
        self.dict_factor_params: dict = {"k": {"delta": self.k_delta, "step": self.k_step,
                                               "plan": self.data_plan["unit"]["layer_prop"]["permeability"]},
                                         "h": {"delta": self.h_delta, "step": self.h_step,
                                               "plan": self.data_plan["unit"]["h_eff"]}
                                         }

    def calc_flow_rate(self, input_data: dict):
        model = TaskTarget(**input_data)
        actual = (
            self.calculation_by_laplace(model)
        )["result"]

        return actual["flow_rate_result"]

    @staticmethod
    def calculation_by_laplace(input_data: TaskTarget) -> dict:
        """Отправка параметров для расчета дебита жидкости

        :param input_data: Входные данные для расчета параметров
        :type input_data: src.entities.models.TaskTarget

        :return: Дебит жидкости (в сутки)
        :rtype: dict
        """
        res: dict = {"t": None, "flow_rate_result": None}
        input_data_dict = input_data.model_dump(by_alias=True)
        res["t"] = input_data_dict["target"]["cumulative_time"]
        builder = Builder(**input_data_dict)
        p_bhp = input_data_dict["target"]["target_values"]["p_bhp"]
        try:
            S = builder.calc_S()
            q_d = builder.calc_q_d(S)
            S_new = builder.calc_double_porosity(S)
            p_d = builder.calc_p_d(S_new)
            delta_p = builder.calc_delta_p(p_bhp)
            flow_rate = builder.calc_flow_rate(S, q_d, p_d, delta_p)
            res["flow_rate_result"] = flow_rate
        except Exception as e:
            msg = f"Error calc_flow_rate route: {e}!"
        else:
            msg = "calc_flow_rate was successful!"

        response = {
            "result": res,
            "message": msg,
        }
        return response

    @staticmethod
    def change_values_in_data(k, h, Pwf, fw, Pr, ct, wf, kf, dict_for_calc) -> dict:
        """
        Изменение значений подаваемых параметров в словаре
        :param k: проницаемость, мДа
        :param h: эффективная мощность пласта, м
        :param Pwf: забойное давление, атм
        :param fw: обводненность, %
        :param Pr: пластовое давление, атм
        :param ct: общая сжимаемость, 1/атм
        :param wf: ширина трещины ГРП у ствола скважины, мм
        :param kf: проницаемость трещины ГРП, дарси
        :param dict_for_calc: словарь с входными данными
        :return: словарь с измененными значениями заданных параметров
        """
        dict_for_calc["unit"]["layer_prop"]["permeability"] = k
        dict_for_calc["unit"]["h_eff"] = h
        dict_for_calc["target"]["target_values"]["p_bhp"] = Pwf
        dict_for_calc["unit"]["layer_prop"]["water_cut"] = fw
        dict_for_calc["unit"]["layer_prop"]["p_res_init"] = Pr
        dict_for_calc["unit"]["layer_prop"]["compressibility"] = ct
        dict_for_calc["unit"]["layer_prop"]["grp_prop"]["wellbore_wf"] = wf
        dict_for_calc["unit"]["layer_prop"]["grp_prop"]["kf"] = kf
        return dict_for_calc

    def calc_qkliq(self):
        """
        Рассчитывает дебит жидкости Qkliq.

        :return: факторы по дебиту жидкости QKliq
        """
        q_delta_by_factor = dict.fromkeys(self.dict_factor_params.keys(), None)

        # for param in self.dict_factor_params.keys():

        Qk = np.zeros(self.steps)

        # Основной расчет Qk
        for i in range(self.steps):
            # Q1 =
            k1 = self.k1 + self.k_delta * i / self.steps
            k2 = self.k2 + self.k_delta * i / self.steps
            h1 = self.data_plan["unit"]["h_eff"] + self.h_delta * i / self.steps
            h2 = self.data_plan["unit"]["h_eff"] + self.h_delta * i / self.steps
            Pwf1 = self.data_plan["target"]["target_values"]["p_bhp"] + self.Pwf_delta * i / self.steps
            Pwf2 = self.data_plan["target"]["target_values"]["p_bhp"] + self.Pwf_delta * i / self.steps
            fw1 = self.data_plan["unit"]["layer_prop"]["water_cut"] + self.fw_delta * i / self.steps
            fw2 = self.data_plan["unit"]["layer_prop"]["water_cut"] + self.fw_delta * i / self.steps
            Pr1 = self.data_plan["unit"]["layer_prop"]["p_res_init"] + self.Pr_delta * i / self.steps
            Pr2 = self.data_plan["unit"]["layer_prop"]["p_res_init"] + self.Pr_delta * i / self.steps
            ct1 = self.data_plan["unit"]["layer_prop"]["compressibility"] + self.ct_delta * i / self.steps
            ct2 = self.data_plan["unit"]["layer_prop"]["compressibility"] + self.ct_delta * i / self.steps
            wf1 = self.data_plan["unit"]["layer_prop"]["grp_prop"]["wellbore_wf"] + self.wf_delta * i / self.steps
            wf2 = self.data_plan["unit"]["layer_prop"]["grp_prop"]["wellbore_wf"] + self.wf_delta * i / self.steps
            kf1 = self.data_plan["unit"]["layer_prop"]["grp_prop"]["kf"] + self.kf_delta * i / self.steps
            kf2 = self.data_plan["unit"]["layer_prop"]["grp_prop"]["kf"] + self.kf_delta * i / self.steps

            # Рассчитываем дебит для текущего шага
            Qk[i] = ((self.calc_flow_rate(self.change_values_in_data(k1, h1, Pwf1, fw1, Pr1, ct1, wf1, kf1,
                                                                     copy.deepcopy(self.data_plan))) -
                      self.calc_flow_rate(self.change_values_in_data(k2, h2, Pwf2, fw2, Pr2, ct2, wf2, kf2,
                                                                     copy.deepcopy(self.data_plan)))) / 2 /
                     self.dict_factor_params["k"]["step"])

        # Интеграл по Qk
        delta = 0

        for g in range(self.steps - 1):
            delta += (((Qk[g + 1]) + Qk[g]) / 2) * (
                    self.dict_factor_params["k"]["delta"] / self.steps)

        q_delta_by_factor["k"] = delta

        # for param in self.dict_factor_params.keys():
        Qk = np.zeros(self.steps)

        # Основной расчет Qk
        for i in range(self.steps):
            k1 = self.data_plan["unit"]["layer_prop"]["permeability"] + self.k_delta * i / self.steps
            k2 = self.data_plan["unit"]["layer_prop"]["permeability"] + self.k_delta * i / self.steps
            h1 = self.h1 + self.h_delta * i / self.steps
            h2 = self.h2 + self.h_delta * i / self.steps
            Pwf1 = self.data_plan["target"]["target_values"]["p_bhp"] + self.Pwf_delta * i / self.steps
            Pwf2 = self.data_plan["target"]["target_values"]["p_bhp"] + self.Pwf_delta * i / self.steps
            fw1 = self.data_plan["unit"]["layer_prop"]["water_cut"] + self.fw_delta * i / self.steps
            fw2 = self.data_plan["unit"]["layer_prop"]["water_cut"] + self.fw_delta * i / self.steps
            Pr1 = self.data_plan["unit"]["layer_prop"]["p_res_init"] + self.Pr_delta * i / self.steps
            Pr2 = self.data_plan["unit"]["layer_prop"]["p_res_init"] + self.Pr_delta * i / self.steps
            ct1 = self.data_plan["unit"]["layer_prop"]["compressibility"] + self.ct_delta * i / self.steps
            ct2 = self.data_plan["unit"]["layer_prop"]["compressibility"] + self.ct_delta * i / self.steps
            wf1 = self.data_plan["unit"]["layer_prop"]["grp_prop"]["wellbore_wf"] + self.wf_delta * i / self.steps
            wf2 = self.data_plan["unit"]["layer_prop"]["grp_prop"]["wellbore_wf"] + self.wf_delta * i / self.steps
            kf1 = self.data_plan["unit"]["layer_prop"]["grp_prop"]["kf"] + self.kf_delta * i / self.steps
            kf2 = self.data_plan["unit"]["layer_prop"]["grp_prop"]["kf"] + self.kf_delta * i / self.steps

            # Рассчитываем дебит для текущего шага
            Qk[i] = ((self.calc_flow_rate(self.change_values_in_data(k1, h1, Pwf1, fw1, Pr1, ct1, wf1, kf1,
                                                                     copy.deepcopy(self.data_fact))) -
                      self.calc_flow_rate(self.change_values_in_data(k2, h2, Pwf2, fw2, Pr2, ct2, wf2, kf2,
                                                                     copy.deepcopy(self.data_fact)))) / 2 /
                     self.dict_factor_params["h"]["step"])

        # Интеграл по Qk
        delta = 0

        for g in range(self.steps - 1):
            delta += (((Qk[g + 1]) + Qk[g]) / 2) * (
                    self.dict_factor_params["h"]["delta"] / self.steps)

        q_delta_by_factor["h"] = delta

        return q_delta_by_factor
