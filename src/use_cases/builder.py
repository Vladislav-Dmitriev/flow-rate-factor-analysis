import numpy as np
from math import log, exp, pi

from src.entities.Well_model import Wellbore
from src.entities.Layer_model import Layer
from src.entities.Vert_pp_calculator import Vert_pp_Calculator
from src.entities.Laplace_calculator import LaplaceCalculator
from src.entities.models import LayerProperty, LaplaceAccuracyOptions
from src.entities.Frac_calculator import Frac_Calculator
from src.entities.Multifrac_calculator import Multifrac_Calculator
from src.entities.Multilateral_calculator import Multilateral_Calculator


class Builder:

    def __init__(self, **input_data):
        """
        :param input_data: входные данные
        """
        k = input_data["unit"]["layer_prop"]["permeability"]
        h = input_data["unit"]["h_eff"]
        ct = input_data["unit"]["layer_prop"]["compressibility"]
        skin = input_data["unit"]["skin"]
        phi = input_data["unit"]["layer_prop"]["porosity"]
        w = input_data["unit"]["layer_prop"]["ye"]
        wf = input_data["unit"]["layer_prop"]["wc_rectangle_ratio"]
        l = input_data["unit"]["layer_prop"]["xe"]
        lf = input_data["unit"]["layer_prop"]["lc_ratio"]
        zw_h = input_data["unit"]["vertical_offset"]
        xbound = input_data["base_prop"]["border_type_x"]
        ybound = input_data["base_prop"]["border_type_y"]
        zbound_up = input_data["base_prop"]["border_type_z_up"]
        zbound_down = input_data["base_prop"]["border_type_z_down"]
        rw = input_data["unit"]["wellbore_prop"]["wellbore_r"]
        mu = input_data["unit"]["layer_prop"]["viscosity_oil"]
        bl = input_data["unit"]["layer_prop"]["b_oil"]
        reservoir_model = input_data["unit"]["layer_prop"]["res_model_type"]
        l_hor = input_data["unit"]["wellbore_prop"]["horizontal_wellbore_length"]
        kv_kh = input_data["unit"]["layer_prop"]["kv_kh_ratio"]
        wellbore_type = input_data["unit"]["wellbore_prop"]["wellbore_type"]
        perf = input_data["unit"]["wellbore_prop"]["horizontal_wellbore_perf_ratio"]
        n_perf = input_data["unit"]["wellbore_prop"]["horizontal_perf_count"]
        well_permeability = input_data["unit"]["wellbore_prop"]["permeability"]
        hw_f = input_data["unit"]["perfres_ratio"]
        self.num_segments = input_data["base_prop"]["segment_count"]
        self.calc_type = input_data["base_prop"]["calc_type"]
        self.p_res = input_data["unit"]["layer_prop"]["p_res_init"]
        self.water_cut = input_data["unit"]["layer_prop"]["water_cut"]
        self.p_b = input_data["unit"]["layer_prop"]["p_bubble"]
        self.first_time_step = input_data["target"]["time_step"]
        self.cumulative_work_time = input_data["target"]["cumulative_work_time"]
        self.number_of_steps = input_data["target"]["number_of_steps"]
        self.grp_flag = input_data["unit"]["layer_prop"]["grp_flag"]
        self.mgrp_flag = input_data["unit"]["layer_prop"]["mgrp_flag"]
        rel_M = input_data["unit"]["layer_prop"]["kmu_in_out_ratio"]
        rel_D = input_data["unit"]["layer_prop"]["kmuphict_in_out_ratio"]
        R_in = input_data["unit"]["layer_prop"]["internal_r"]
        f_compressibility = input_data["unit"]["layer_prop"]["f_compressibility"]
        f_porosity = input_data["unit"]["layer_prop"]["f_porosity"]
        lambda_ = input_data["unit"]["layer_prop"]["lambda"]
        # константы
        extend_reflections_mode = LayerProperty().extend_reflections_mode
        tiny = LaplaceAccuracyOptions().tiny
        tiny_2 = LaplaceAccuracyOptions().tiny_2
    
        large_s = LaplaceAccuracyOptions().large_s
        small = LaplaceAccuracyOptions().small
        part_sum_num = LaplaceAccuracyOptions().part_sum_num
        max_it = LaplaceAccuracyOptions().max_it
        small_bessel_arg = LaplaceAccuracyOptions().small_bessel_arg
        part_sum_num_b2_1 = LaplaceAccuracyOptions().part_sum_num_b2_2
        part_sum_num_f_4 = LaplaceAccuracyOptions().part_sum_num_f_4
        self.number_of_lapl_coeff = LaplaceAccuracyOptions().number_of_lapl_coeff
        self.laplace_calculator = LaplaceCalculator(extend_reflections_mode, tiny, tiny_2, large_s, small,
                                                    part_sum_num, max_it, small_bessel_arg, part_sum_num_b2_1,
                                                    part_sum_num_f_4)
        self.v = self.laplace_calculator.coef(self.number_of_lapl_coeff)
        self.vert_pp_calculator = Vert_pp_Calculator(extend_reflections_mode, tiny, tiny_2, large_s, small,
                                                     part_sum_num, max_it, small_bessel_arg, part_sum_num_b2_1,
                                                     part_sum_num_f_4)
        self.unit_length, rw, skin, r = self.calc_unit_length(rw, wellbore_type, l_hor, skin)
        array_well_size = np.zeros(2)
        if wellbore_type == "multilateral":
            n_lateral = input_data["unit"]["layer_prop"]["multilateral_prop"]["n_lateral"]
            l_lateral = input_data["unit"]["layer_prop"]["multilateral_prop"]["l_lateral"]
            psi_lateral = input_data["unit"]["layer_prop"]["multilateral_prop"]["psi_lateral"]
            self.multilateral_calculator = Multilateral_Calculator(extend_reflections_mode, tiny, tiny_2, large_s,
                                                                   small,
                                                                   part_sum_num, max_it, small_bessel_arg,
                                                                   part_sum_num_b2_1,
                                                                   part_sum_num_f_4, n_lateral, l_lateral, psi_lateral)
            array_well_size = self.multilateral_calculator.calc_multilateral_geometry(l_hor)
        if self.grp_flag:
            hf = input_data["unit"]["layer_prop"]["grp_prop"]["hf"]
            fracture_grow_t = input_data["unit"]["layer_prop"]["grp_prop"]["fracture_grow_t"]
            wellbore_wf = input_data["unit"]["layer_prop"]["grp_prop"]["wellbore_wf"]
            res_wf = input_data["unit"]["layer_prop"]["grp_prop"]["res_wf"]
            grade = input_data["unit"]["layer_prop"]["grp_prop"]["grade"]
            kf = input_data["unit"]["layer_prop"]["grp_prop"]["kf"]
            S_choke = input_data["unit"]["layer_prop"]["grp_prop"]["skin_ch"]
            sf = input_data["unit"]["layer_prop"]["grp_prop"]["skin_border"]
            # экземпляр класса ГРП
            self.frac = Frac_Calculator(wellbore_wf, res_wf, grade, kf, hf, fracture_grow_t, rel_M, rel_D, R_in,
                                        S_choke, sf,
                                        k, hw_f, extend_reflections_mode,
                                        tiny, tiny_2, large_s, small,
                                        part_sum_num, max_it, small_bessel_arg, part_sum_num_b2_1,
                                        part_sum_num_f_4)
            self.unit_length = hf
            array_well_size[0] = 2 * hf
            if self.mgrp_flag:
                n_frac = input_data["unit"]["layer_prop"]["mgrp_prop"]["grp_count"]
                f_direction = input_data["unit"]["layer_prop"]["mgrp_prop"]["f_direction"]
                self.multifrac_calculator = Multifrac_Calculator(wellbore_wf, res_wf, grade, kf, hf, fracture_grow_t,
                                                                 rel_M, rel_D, R_in,
                                                                 S_choke, sf,
                                                                 k, hw_f, extend_reflections_mode,
                                                                 tiny, tiny_2, large_s, small,
                                                                 part_sum_num, max_it, small_bessel_arg,
                                                                 part_sum_num_b2_1,
                                                                 part_sum_num_f_4, n_frac, f_direction)
                array_well_size = self.multifrac_calculator.calc_well_size(l_hor)
        # экземпляр пласта
        self.layer = Layer(k, h, ct, skin, phi, w, wf, l, lf, zw_h, xbound, ybound, zbound_up, zbound_down, kv_kh, rw,
                           hw_f, reservoir_model, self.unit_length, f_compressibility, f_porosity, lambda_, self.mgrp_flag, r)
        # экземпляр скважины
        self.wellbore = Wellbore(rw, mu, bl, l_hor, self.unit_length, h, phi, ct,
                                 wellbore_type, perf, n_perf, well_permeability, k)
        if self.wellbore.wellbore_type == "horizontal":
            array_well_size[0] = l_hor
        self.multT = self.calc_multT()
        if not self.wellbore.well_fits_rectangle(self.layer.l, self.layer.w, array_well_size):
            raise ValueError("Скважина не помещается в прямоугольную область дренирования. Введите другие "
                             "значения длины и ширины прямоугольника.")

    def calc_t(self):
        """
        Рассчитывает накопленное время работы
        :return: накопленное время работы
        """
        time_step_multiplier = exp(
            1 / (self.number_of_steps - 1) * log(self.cumulative_work_time / self.first_time_step)
        )
        t = np.zeros(self.number_of_steps)
        t[0] = self.first_time_step
        t[1:] = time_step_multiplier
        t = np.multiply.accumulate(t)
        return t

    def calc_multT(self):
        """
        Рассчитывает multT
        :return: множитель multT
        """
        if self.layer.reservoir_model == "DP Slab":
            multT = (self.layer.omega
                     * self.layer.ct
                     * self.layer.phi
                     * self.wellbore.mu
                     * self.unit_length ** 2
                     / 0.00036
                     / self.layer.k
                     )
        else:
            multT = (
                    self.layer.ct
                    * self.layer.phi
                    * self.wellbore.mu
                    * self.unit_length ** 2
                    / 0.00036
                    / self.layer.k
            )
        return multT

    def calc_unit_length(self, rw, wellbore_type, l_hor, skin):
        """
        Рассчитывает unit_length
        :param rw: радиус скважины
        :param wellbore_type: тип скважины
        :param l_hor: длина горизонтального ствола
        :param skin: скин-фактор
        :return: единица длины
        """
        if wellbore_type == "vertical":
            rw = rw * exp(-skin)
            unit_length = rw
            r = rw
        elif wellbore_type == "horizontal" or wellbore_type == "multilateral":
            rw = rw * exp(-skin)
            unit_length = l_hor / 2
            r = 0
        skin = 0
        return unit_length, rw, skin, r

    def calc_p_bhp(self, S, q_d, p_d, ql):
        """
        Рассчитывает забойное давление
        :param S: переменные пространства Лапласа
        :param q_d: вектор безразмерного дебита
        :param p_d: вектор безразмерного давления
        :param ql: целевой дебит жидкости, m3/день
        :return: забойное давление
        """
        res_mult = (
                18.42
                * self.wellbore.mu
                * ql
                * self.wellbore.bl
                / (self.layer.k * self.layer.h)
        )
        p_wf = []
        delta_p = []
        for i in range(0, self.number_of_steps):
            SumR = 0
            for j in range(1, self.number_of_lapl_coeff + 1):
                S_ = S[(i, j)]
                p_d_ = self.wellbore.calc_pd_for_stable_model(
                    S_, p_d[(i, j)], self.wellbore.cd
                )
                add = (self.v[j] / j) * (p_d_ * q_d[(i, j)] * S_)
                SumR += add
            p_bhp_t = SumR * res_mult
            p_bhp_t = float(p_bhp_t)
            delta_p.append(p_bhp_t)
            p_wf_t = self.p_res - p_bhp_t
            p_wf.append(p_wf_t)
        return p_wf, delta_p

    def calc_S(self, t):
        """
        Рассчитывает переменные пространства Лапласа
        :param t: накопленное время работы
        :return: переменные пространства Лапласа
        """
        DlogTW = log(2)
        td = t / self.multT
        S = np.zeros((self.number_of_steps, self.number_of_lapl_coeff + 1))
        for i in range(self.number_of_steps):
            for j in range(1, self.number_of_lapl_coeff + 1):
                S[(i, j)] = j * DlogTW / td[i]
        return S

    def calc_double_porosity(self, S):
        """
        Пересчитывает переменные пространства Лапласа в зависимости от модели пласта
        :param S: переменные пространства Лапласа
        :return: переменные пространства Лапласа
        """
        if self.layer.reservoir_model == "Homogeneous" or self.layer.reservoir_model == "composite":
            return S
        elif self.layer.reservoir_model == "DP PSS":
            S = S * (S * self.layer.omega * (1 - self.layer.omega) + self.layer.lambda_) / (
                    S * (1 - self.layer.omega) + self.layer.lambda_)
        elif self.layer.reservoir_model == "DP Slab":
            omega_ = (1 - self.layer.omega) / self.layer.omega
            S = S * (1 + np.sqrt(self.layer.lambda_ * omega_ / 3 / S) * np.tanh(
                np.sqrt(3 * omega_ * S / self.layer.lambda_)))

        return S

    def calc_q_d(self, S):
        """
        Рассчитывает вектор безразмерного дебита
        :param S: переменные пространства Лапласа
        :return: вектор безразмерного дебита
        """
        S[S == 0] = 1e-10
        q_d = (2 * pi) / S
        return q_d

    def calc_p_d(self, S):
        """
        Рассчитывает вектор безразмерного давления
        :param S: переменные пространства Лапласа
        :return: вектор безразмерного давления
        """
        pd = np.zeros((self.number_of_steps, self.number_of_lapl_coeff + 1))
        for i in range(self.number_of_steps):
            for j in range(1, self.number_of_lapl_coeff + 1):
                if self.wellbore.wellbore_type == 'vertical':
                    if self.grp_flag:
                        pd[(i, j)] = self.frac.calc_frac(self.calc_type, S[i, j], self.layer.xd, self.layer.yd,
                                                         self.layer.xwd, self.layer.ywd, self.layer.xed, self.layer.yed,
                                                         self.layer.xbound, self.layer.ybound,
                                                         self.layer.reservoir_model,
                                                         self.layer.zd, self.layer.zwd, self.layer.zed, self.layer.hwd,
                                                         self.layer.zbound_up, self.layer.zbound_down,
                                                         self.num_segments, self.multT, self.layer.k, self.layer.hw_f,
                                                         self.unit_length)
                    elif self.layer.hw_f != 1 or (
                            self.layer.hw_f == 1 and (self.layer.zbound_up == "c" or self.layer.zbound_down == "c")):
                        pd[(i, j)] = self.vert_pp_calculator.partial_penetration_vert_rect_lapl(S[(i, j)],
                                                                                                self.layer.xd,
                                                                                                self.layer.yd,
                                                                                                self.layer.zd,
                                                                                                self.layer.xwd,
                                                                                                self.layer.ywd,
                                                                                                self.layer.zwd,
                                                                                                self.layer.xed,
                                                                                                self.layer.yed,
                                                                                                self.layer.zed,
                                                                                                self.layer.hwd,
                                                                                                self.layer.xbound,
                                                                                                self.layer.ybound,
                                                                                                self.layer.zbound_up,
                                                                                                self.layer.zbound_down,
                                                                                                self.layer.rwd,
                                                                                                self.calc_type)
                    else:
                        pd[(i, j)] = self.laplace_calculator.calc_pd(
                            S[(i, j)],
                            self.layer.xd,
                            self.layer.yd,
                            self.layer.xwd,
                            self.layer.ywd,
                            self.layer.xed,
                            self.layer.yed,
                            self.layer.xbound,
                            self.layer.ybound,
                            self.layer.skin,
                        )
                elif self.wellbore.wellbore_type == 'horizontal':
                    if self.mgrp_flag:
                        xw, yw, xe, ye, xbound, ybound = self.multifrac_calculator.calc_multifrac_geometry(
                            self.layer.xe, self.layer.ye, self.layer.wf, self.layer.lf, self.wellbore.l_hor,
                            self.layer.xbound, self.layer.ybound)
                        xwd, ywd, zwd, xed, yed, zed, rwd = self.layer.calc_dimensionless_geometry(xw, yw,
                                                                                                   self.layer.zw, xe,
                                                                                                   ye, self.layer.ze,
                                                                                                   self.wellbore.rw,
                                                                                                   self.layer.s_kv_kh,
                                                                                                   self.unit_length)
                        pd[(i, j)] = self.multifrac_calculator.calc_multifrac(self.layer.reservoir_model, S[i, j],
                                                                              xwd, ywd,
                                                                              xed, yed,
                                                                              self.num_segments,
                                                                              xbound, ybound,
                                                                              self.calc_type, self.unit_length)
                    else:
                        pd[(i, j)] = self.laplace_calculator.horizontal_rect_lapl(
                            S[i, j], self.layer.rwd, self.layer.zed,
                            self.layer.yd, self.layer.ywd,
                            self.layer.zd, self.layer.zwd,
                            self.layer.xwd, self.layer.xd,
                            self.layer.xed, self.layer.yed,
                            self.layer.xbound, self.layer.ybound,
                            self.layer.zbound_up,
                            self.layer.zbound_down,
                            self.wellbore.fcd, self.layer.skin,
                            self.wellbore.perf, self.wellbore.n_perf,
                            self.num_segments, self.calc_type, j
                        )
                elif self.wellbore.wellbore_type == "multilateral":
                    pd[(i, j)] = self.multilateral_calculator.multilateral_rect_lapl(S[i, j], self.layer.xwd,
                                                                                     self.layer.ywd, self.layer.zwd,
                                                                                     self.layer.xed, self.layer.yed,
                                                                                     self.layer.zed, self.layer.rwd,
                                                                                     self.layer.xbound,
                                                                                     self.layer.ybound,
                                                                                     self.layer.zbound_up,
                                                                                     self.layer.zbound_down,
                                                                                     self.calc_type, self.num_segments,
                                                                                     self.unit_length)

        return pd

    def calc_flow_rate(self, S, q_d, p_d, delta_p, is_lift=False):
        """
        Рассчитывает дебит жидкости для заданного забойного давления
        :param S: переменные пространства Лапласа
        :param q_d: вектор безразмерного дебита
        :param p_d: вектор безразмерного давления
        :param delta_p: депрессия с учетом газовой фазы и обводненности продукции
        :param is_lift: учитывать лифт или нет
        :return: дебит жидкости
        """
        res_mult = (
                self.layer.k
                * self.layer.h
                * delta_p
                / (18.42 * self.wellbore.bl * self.wellbore.mu)
        )
        flow_rate = []
        if not is_lift:
            for i in range(0, self.number_of_steps):
                SumR = 0
                for j in range(1, self.number_of_lapl_coeff + 1):
                    add = (self.v[j] / j) / (p_d[(i, j)] * q_d[(i, j)] * S[(i, j)])
                    SumR += add
                flow_rate_t = SumR * res_mult
                flow_rate_t = float(flow_rate_t)
                flow_rate.append(flow_rate_t)
        return flow_rate

    '''Функция для расчета депрессии с учетом обводненности и давления насыщения'''
    def calc_delta_p(self, Pwf):
        """
        Рассчитывает депрессию с учетом газовой фазы и обводненности продукции (метод Вогеля с учетом обводненности).
        :param Pwf: Забойное давление (бар)
        :return: Депрессия (бар)
        """
        # Корректировка давления насыщения для пластов ниже давления насыщения
        if self.p_b > self.p_res:
            self.p_b = self.p_res

        # Если продукция полностью водяная
        if self.water_cut == 100:
            return self.p_res - Pwf

        # Инициализация переменных
        Pwf_G = (4 / 9) * (self.water_cut / 100) * self.p_b  # Давление для нефти с водой и газом

        if Pwf >= self.p_b:
            # Нефть или нефть + вода (нет газа)
            return self.p_res - Pwf

        elif Pwf < Pwf_G:
            # Нефть + вода + газ
            tgb_r = (81 - 80 * (0.999 * self.p_b - 0.0018 * (self.p_res - self.p_b)) / self.p_b) ** 0.5
            tgb = (
                    self.water_cut / 100 +
                    (0.125 * (1 - self.water_cut / 100) * self.p_b * (-1 + tgb_r)) / (0.001 * (self.p_res - (4 / 9) * self.p_b))
            )
            return (Pwf_G + (self.p_res - (4 / 9) * self.p_b) * tgb - Pwf) / tgb

        else:
            # Смешанная продукция (нефть + газ или нефть + вода + газ)
            a = (Pwf + 0.125 * (1 - self.water_cut / 100) * self.p_b - (self.water_cut / 100) * self.p_res) / (0.125 * (1 - self.water_cut / 100) * self.p_b)
            B = (self.water_cut / 100) / (0.125 * (1 - self.water_cut / 100) * self.p_b)
            c = 2 * a * B + 144 / self.p_b
            d = a ** 2 - 144 * (self.p_res - self.p_b) / self.p_b - 81

            if B == 0:
                # Нефть + газ (по Вогелю)
                return -d / c
            else:
                # Нефть + вода + газ
                return (-c + (c ** 2 - 4 * (B ** 2) * d) ** 0.5) / (2 * (B ** 2))

    def calc_q_storage(self, S, q_d, p_d, ql):
        """
        Рассчитывает дебит жидкости
        :param S: переменные пространства Лапласа
        :param q_d: вектор безразмерного дебита
        :param p_d: вектор безразмерного давления
        :param ql: целевой дебит жидкости, m3/день
        :return: дебит жидкости
        """
        mult_q_wbs = ql
        res_mult = mult_q_wbs
        q_storage = []
        for i in range(0, self.number_of_steps):
            SumR = 0
            for j in range(1, self.number_of_lapl_coeff + 1):
                S_ = S[(i, j)]
                p_d_ = self.wellbore.calc_pd_for_stable_model(
                    S_, p_d[(i, j)], self.wellbore.cd
                )
                wbs_f_d = 0
                add = (
                        (self.v[j] / j)
                        * (S_ ** 2)
                        * self.wellbore.cd
                        * (p_d_ * q_d[(i, j)] - wbs_f_d)
                )
                SumR += add
            q_storage_t = SumR * res_mult
            q_storage_t = float(q_storage_t)
            q_storage.append(q_storage_t)
        return q_storage

    def calc_JD_for_pressure(self, t, ql, q_storage, delta_p):
        """
        Рассчитывает коэффициент продуктивности (JD) для давления
        :param t: накопленное время работы
        :param ql: целевой дебит жидкости, m3/день
        :param q_storage: дебит жидкости, рассчитанный по режиму 10
        :param delta_p: давление
        :return: JD
        """
        JD_list = []
        for i in range(0, self.number_of_steps):
            Q = ql - q_storage[i]
            V_pore = self.layer.w * self.layer.l * self.layer.h * self.layer.phi
            np_reservoir = t[i] * Q / 24
            if self.layer.xbound == "c" or self.layer.ybound == "c" or self.layer.zbound_up == "c" or self.layer.zbound_down == "c":
                P_avg = self.p_res
            else:
                P_avg = (
                        self.p_res - np_reservoir * self.wellbore.bl / V_pore / self.layer.ct
                )
            Pwf = self.p_res - delta_p[i]
            JD = (
                    18.42
                    * Q
                    * self.wellbore.mu
                    * self.wellbore.bl
                    / (self.layer.k * self.layer.h * (P_avg - Pwf))
            )
            JD_list.append(JD)
        return JD_list

    def calc_cumulative_prod(self, S, q_d, p_d, p_bhp, is_lift=False):
        """
        Рассчитывает накопленную добычу жидкости для заданного забойного давления
        :param S: переменные пространства Лапласа
        :param q_d: вектор безразмерного дебита
        :param p_d: вектор безразмерного давления
        :param p_bhp: целевое забойное давление
        :param is_lift: учитывать лифт или нет
        :return: дебит жидкости
        """
        delta_p = self.p_res - p_bhp
        multq = (
                self.layer.k
                * self.layer.h
                * delta_p
                / (18.42 * self.wellbore.bl * self.wellbore.mu)
        )
        res_mult = multq * self.multT / 24
        cumulative_prod = []
        if not is_lift:
            for i in range(0, self.number_of_steps):
                SumR = 0
                for j in range(1, self.number_of_lapl_coeff + 1):
                    add = (self.v[j] / j) / (
                            p_d[(i, j)] * q_d[(i, j)] * S[(i, j)] ** 2
                    )
                    SumR += add
                flow_rate_t = SumR * res_mult
                flow_rate_t = float(flow_rate_t)
                cumulative_prod.append(flow_rate_t)
        return cumulative_prod

    def calc_JD_for_flow_rate(self, p_bhp, q_calc, cumulative_prod):
        """
        Рассчитывает коэффициент продуктивности (JD) для дебита
        :param p_bhp: целевое забойное давление
        :param q_calc: дебит жидкости
        :param cumulative_prod: накопленная добыча жидкости
        :return: JD
        """
        JD_list = []
        for i in range(0, self.number_of_steps):
            Q = q_calc[i]
            V_pore = self.layer.w * self.layer.l * self.layer.h * self.layer.phi
            np_reservoir = cumulative_prod[i]
            if self.layer.xbound == "c" or self.layer.ybound == "c" or self.layer.zbound_up == "c" or self.layer.zbound_down == "c":
                P_avg = self.p_res
            else:
                P_avg = (
                        self.p_res - np_reservoir * self.wellbore.bl / V_pore / self.layer.ct
                )
            Pwf = p_bhp
            JD = (
                    18.42
                    * Q
                    * self.wellbore.mu
                    * self.wellbore.bl
                    / (self.layer.k * self.layer.h * (P_avg - Pwf))
            )
            JD_list.append(JD)

        return JD_list
