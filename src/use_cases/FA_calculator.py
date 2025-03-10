import copy
import numpy as np
from src.entities.models import TaskTarget
from src.use_cases.builder import Builder


class FactorAnalysis:

    def __init__(self, data_fact: dict, data_plan: dict, steps: int):
        # Разницы параметров
        self.k_delta = data_fact["unit"]["layer_prop"]["permeability"] - data_plan["unit"]["layer_prop"]["permeability"]
        self.h_delta = data_fact["unit"]["h_eff"] - data_plan["unit"]["h_eff"]
        self.mu_delta = (data_fact["unit"]["layer_prop"]["viscosity_eff"] -
                         data_plan["unit"]["layer_prop"]["viscosity_eff"])
        self.Pbhp_delta = data_fact["target"]["target_values"]["p_bhp"] - data_plan["target"]["target_values"]["p_bhp"]
        self.fw_delta = data_fact["unit"]["layer_prop"]["water_cut"] - data_plan["unit"]["layer_prop"]["water_cut"]
        self.Pr_delta = data_fact["unit"]["layer_prop"]["p_res_init"] - data_plan["unit"]["layer_prop"]["p_res_init"]
        self.ct_delta = (data_fact["unit"]["layer_prop"]["compressibility"] -
                         data_plan["unit"]["layer_prop"]["compressibility"])
        self.wf_delta = (data_fact["unit"]["layer_prop"]["grp_prop"]["wellbore_wf"] -
                         data_plan["unit"]["layer_prop"]["grp_prop"]["wellbore_wf"])
        self.kf_delta = (data_fact["unit"]["layer_prop"]["grp_prop"]["kf"] -
                         data_plan["unit"]["layer_prop"]["grp_prop"]["kf"])
        self.b_oil_delta = data_fact["unit"]["layer_prop"]["b_eff"] - data_plan["unit"]["layer_prop"]["b_eff"]
        self.xf_delta = (data_fact["unit"]["layer_prop"]["grp_prop"]["xf"] -
                         data_plan["unit"]["layer_prop"]["grp_prop"]["xf"])
        # Кол-во шагов
        self.steps = steps
        # Шаги по параметрам
        self.k_step = data_plan["unit"]["layer_prop"]["permeability"] / self.steps
        self.h_step = data_plan["unit"]["h_eff"] / self.steps
        self.mu_step = data_plan["unit"]["layer_prop"]["viscosity_eff"] / self.steps
        self.Pbhp_step = data_plan["target"]["target_values"]["p_bhp"] / self.steps
        self.fw_step = data_plan["unit"]["layer_prop"]["water_cut"] / self.steps
        self.Pr_step = data_plan["unit"]["layer_prop"]["p_res_init"] / self.steps
        self.ct_step = data_plan["unit"]["layer_prop"]["compressibility"] / self.steps
        self.wf_step = data_plan["unit"]["layer_prop"]["grp_prop"]["wellbore_wf"] / self.steps
        self.kf_step = data_plan["unit"]["layer_prop"]["grp_prop"]["kf"] / self.steps
        self.b_oil_step = data_plan["unit"]["layer_prop"]["b_eff"] / self.steps
        self.xf_step = data_plan["unit"]["layer_prop"]["grp_prop"]["xf"] / self.steps
        # Первые шаги для расчета
        self.k1 = data_plan["unit"]["layer_prop"]["permeability"] + self.k_step
        self.k2 = data_plan["unit"]["layer_prop"]["permeability"] - self.k_step
        self.h1 = data_plan["unit"]["h_eff"] + self.h_step
        self.h2 = data_plan["unit"]["h_eff"] - self.h_step
        self.mu1 = data_plan["unit"]["layer_prop"]["viscosity_eff"] + self.mu_step
        self.mu2 = data_plan["unit"]["layer_prop"]["viscosity_eff"] - self.mu_step
        self.P_bhp1 = data_plan["target"]["target_values"]["p_bhp"] + self.Pbhp_step
        self.P_bhp2 = data_plan["target"]["target_values"]["p_bhp"] - self.Pbhp_step
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
        self.b_oil1 = data_plan["unit"]["layer_prop"]["b_eff"] + self.b_oil_step
        self.b_oil2 = data_plan["unit"]["layer_prop"]["b_eff"] - self.b_oil_step
        self.xf1 = data_plan["unit"]["layer_prop"]["grp_prop"]["xf"] + self.xf_step
        self.xf2 = data_plan["unit"]["layer_prop"]["grp_prop"]["xf"] - self.xf_step
        # Входные данные
        self.data_plan = data_plan
        self.data_fact = data_fact
        # Параметры для факторного анализа
        self.dict_factor_params: dict = {"k": {"delta": self.k_delta, "step": self.k_step,
                                               "plan": self.data_plan["unit"]["layer_prop"]["permeability"],
                                               "path": ["unit", "layer_prop", "permeability"]},
                                         "h_eff": {"delta": self.h_delta, "step": self.h_step,
                                                   "plan": self.data_plan["unit"]["h_eff"], "path": ["unit", "h_eff"]},
                                         "mu": {"delta": self.mu_delta, "step": self.mu_step,
                                                "plan": self.data_plan["unit"]["layer_prop"]["viscosity_eff"],
                                                "path": ["unit", "layer_prop", "viscosity_eff"]},
                                         "P_bhp": {"delta": self.Pbhp_delta, "step": self.Pbhp_step,
                                                   "plan": data_plan["target"]["target_values"]["p_bhp"],
                                                   "path": ["target", "target_values", "p_bhp"]},
                                         "fw": {"delta": self.fw_delta, "step": self.fw_step,
                                                "plan": data_plan["unit"]["layer_prop"]["water_cut"],
                                                "path": ["unit", "layer_prop", "water_cut"]},
                                         "Pr": {"delta": self.Pr_delta, "step": self.Pr_step,
                                                "plan": data_plan["unit"]["layer_prop"]["p_res_init"],
                                                "path": ["unit", "layer_prop", "p_res_init"]},
                                         "ct": {"delta": self.ct_delta, "step": self.ct_step,
                                                "plan": data_plan["unit"]["layer_prop"]["compressibility"],
                                                "path": ["unit", "layer_prop", "compressibility"]},
                                         "wf": {"delta": self.wf_delta, "step": self.wf_step,
                                                "plan": data_plan["unit"]["layer_prop"]["grp_prop"]["wellbore_wf"],
                                                "path": ["unit", "layer_prop", "grp_prop", "wellbore_wf"]},
                                         "kf": {"delta": self.kf_delta, "step": self.kf_step,
                                                "plan": data_plan["unit"]["layer_prop"]["grp_prop"]["kf"],
                                                "path": ["unit", "layer_prop", "grp_prop", "kf"]},
                                         "b_eff": {"delta": self.b_oil_delta, "step": self.b_oil_step,
                                                   "plan": data_plan["unit"]["layer_prop"]["b_eff"],
                                                   "path": ["unit", "layer_prop", "b_eff"]},
                                         "xf": {"delta": self.xf_delta, "step": self.xf_step,
                                                "plan": data_plan["unit"]["layer_prop"]["grp_prop"]["xf"],
                                                "path": ["unit", "layer_prop", "grp_prop", "xf"]},
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

        :return: Дебит нефти (в сутки)
        :rtype: dict
        """
        res: dict = {"t": None, "flow_rate_result": None}
        input_data_dict = input_data.model_dump(by_alias=True)
        res["t"] = input_data_dict["target"]["cumulative_time"]
        builder = Builder(**input_data_dict)
        p_bhp = input_data_dict["target"]["target_values"]["p_bhp"]
        water_cut = input_data_dict["unit"]["layer_prop"]["water_cut"]
        try:
            S = builder.calc_S()
            q_d = builder.calc_q_d(S)
            S_new = builder.calc_double_porosity(S)
            p_d = builder.calc_p_d(S_new)
            delta_p = builder.calc_delta_p(p_bhp)
            flow_rate = builder.calc_flow_rate(S, q_d, p_d, delta_p) * (1 - water_cut / 100) * 0.832
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

    def set_value(self, param, val, current_dict):
        """
        Изменение значения по ключу словаря
        :param param: параметр, значение которого нужно изменить в current_dict
        :param val: новое значение параметра
        :param current_dict: текущий словарь с параметрами расчета
        :return: обновленный dict с параметрами расчета
        """
        for key in self.dict_factor_params[param]["path"][:-1]:
            current_dict = current_dict[key]
        current_dict[self.dict_factor_params[param]["path"][-1]] = val

    def change_values(self, i, param, dict_params1, dict_params2):
        """
        Итерация по факторам и изменение словарей входных параметров на текущий шаг расчета
        :param i: номер шага
        :param param: имя параметра
        :param dict_params1: dict набор параметров для расчета левой точки отрезка итерации по оси значений
        :param dict_params2: dict набор параметров для расчета правой точки отрезка итерации по оси значений
        :return: dict, dict - наборы параметров с обновленными значениями для текущего шага расчета
        """
        for p in self.dict_factor_params.keys():
            if p == param:
                val = (self.dict_factor_params[p]["plan"] - self.dict_factor_params[p]["plan"] / self.steps +
                       self.dict_factor_params[p]["delta"] * i / self.steps)
                val2 = (self.dict_factor_params[p]["plan"] + self.dict_factor_params[p]["plan"] / self.steps +
                        self.dict_factor_params[p]["delta"] * i / self.steps)
                self.set_value(p, val, dict_params1)
                self.set_value(p, val2, dict_params2)
            else:
                val = self.dict_factor_params[p]["plan"] + self.dict_factor_params[p]["delta"] * i / self.steps
                self.set_value(p, val, dict_params1)
                self.set_value(p, val, dict_params2)

    def calc_factors(self):
        """
        Расчет значений факторов
        :return: dict со значениями факторов
        """
        q_delta_by_factor = dict.fromkeys(self.dict_factor_params.keys(), None)

        dict_calc1 = copy.deepcopy(self.data_plan)
        dict_calc2 = copy.deepcopy(self.data_plan)
        for param in self.dict_factor_params.keys():
            if self.dict_factor_params[param]["delta"] == 0:
                q_delta_by_factor[param] = float(0)
                continue
            delta = 0
            Qk = np.zeros(self.steps)
            for i in range(self.steps):
                self.change_values(i, param, dict_calc1, dict_calc2)
                Qk[i] = ((self.calc_flow_rate(dict_calc1) - self.calc_flow_rate(dict_calc2)) /
                         (2 * self.dict_factor_params[param]["step"]))

            for g in range(self.steps - 1):
                delta += (((Qk[g + 1]) + Qk[g]) / 2) * (self.dict_factor_params[param]["delta"] / self.steps)

            q_delta_by_factor[param] = delta

        return q_delta_by_factor
