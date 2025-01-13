import numpy as np
from math import pi, exp, sin, cos, tanh, fabs, sqrt
from src.entities.Vert_pp_calculator import Vert_pp_Calculator


class Frac_Calculator(Vert_pp_Calculator):
    """
    Класс расчёта параметров модели скважины с ГРП
    """
    def __init__(self, wellbore_wf, res_wf, grade, kf, hf, fracture_grow_t, rel_M, rel_D, R_in, S_choke, sf, k,
                 hw_f, extend_reflections_mode, tiny, tiny_2, large_s, small, part_sum_num, max_it, small_bessel_arg,
                 part_sum_num_b2_1, part_sum_num_f_4):
        """
        :param wellbore_wf: ширина трещины у ствола скважины, мм
        :param res_wf: ширина трещины на конце, мм
        :param grade: степень аппроксимационной функции
        :param kf: проницаемость трещины, Д
        :param hf: полудлина трещины, м
        :param fracture_grow_t: время роста трещины, ч
        :param rel_M: коэффициент подвижности
        :param rel_D: коэффициент пьезопроводности
        :param R_in: внутренний радиус для composit-модели, м
        :param S_choke: скин за счет схождения
        :param sf: скин на стенке трещины
        :param k: проницаемость пласта, Д
        :param hw_f: величина вскрытия пласта (в долях от мощности)
        """
        super().__init__(extend_reflections_mode, tiny, tiny_2, large_s, small, part_sum_num, max_it, small_bessel_arg,
                         part_sum_num_b2_1, part_sum_num_f_4)
        self.wellbore_wf = wellbore_wf
        self.res_wf = res_wf
        self.grade = grade
        self.kf = kf
        self.hf = hf
        self.fracture_grow_t = fracture_grow_t
        self.rel_M = rel_M
        self.rel_D = rel_D
        self.R_in = R_in
        self.S_choke = S_choke
        self.sf = sf
        self.cfd = self.calc_cfd(k, hw_f)

    def calc_cfd(self, k, hw_f):
        """
        :param k: проницаемость пласта, Д
        :param hw_f: величина вскрытия пласта (в долях от мощности)
        :return: средняя безразмерная проводимость
        """
        if self.wellbore_wf == self.res_wf:
            Afrac_Width = self.wellbore_wf
        else:
            Afrac_Width = self.grade / (self.grade + 1) / (self.wellbore_wf ** self.grade - self.res_wf ** self.grade) * \
                          ((self.wellbore_wf ** (self.grade + 1)) - (self.res_wf ** (self.grade + 1)))
        cfd = self.kf * Afrac_Width / (k * self.hf) * hw_f

        return cfd

    def calc_params(self, multT, n_seg, k, hw_f, unit_length):
        """
        :param multT: множитель по времени
        :param n_seg: количество сегментов
        :param k: проницаемость пласта, Д
        :param hw_f: величина вскрытия пласта (в долях от мощности)
        :param unit_length: единица длины
        :return:
                 Fcd_shape: параметр, характеризующий производительность трещины
                 growth_velocity_d: безразмерная скорость роста трещины
        """
        Frac_Width = self.fracture_shape_width(n_seg, self.hf, self.wellbore_wf, self.res_wf, self.grade)
        Fcd_shape = np.zeros(shape=Frac_Width.shape)
        for i in range(Fcd_shape.shape[0]):
            Fcd_shape[i] = self.kf * Frac_Width[i] / (k * self.hf) * hw_f
        if self.fracture_grow_t == 0:
            growth_velocity = 0
        else:
            growth_velocity = self.hf / self.fracture_grow_t
        growth_velocity_d = growth_velocity * multT / unit_length

        return Fcd_shape, growth_velocity_d

    def calc_frac(self, calc_type, S, xd, yd, xwd, ywd, xed, yed, xbound, ybound, reservoir_model,
                  zd, zwd, zed, hwd, zbound_up, zbound_down, num_segments, multT, k, hw_f, unit_length):
        """
        :param calc_type: тип расчёта
        :param S: переменная пространства Лапласа
        :param xd: безразмерная координата в направлении x
        :param yd: безразмерная координата в направлении y
        :param xwd: безразмерная координата скважины в направлении x
        :param ywd: безразмерная координата скважины в направлении y
        :param xed: безразмерная координата границы x
        :param yed: безразмерная координата границы y
        :param xbound: тип границы по x
        :param ybound: тип границы по y
        :param reservoir_model: модель пласта
        :param zd: безразмерная координата в направлении z
        :param zwd: безразмерная координата скважины в направлении z
        :param zed: безразмерная координата границы z
        :param hwd: безразмерная величина вскрытия пласта
        :param zbound_up: тип границы верхний по z
        :param zbound_down: тип границы нижний по z
        :param num_segments: количество сегментов
        :param multT: множитель по времени
        :param k: проницаемость пласта, Д
        :param hw_f: величина вскрытия пласта (в долях от мощности)
        :param unit_length: единица длины
        :return: p_d
        """
        if reservoir_model != "composite":
            Fcd_shape, growth_velocity_d = self.calc_params(multT, num_segments, k, hw_f, unit_length)
            p_d = self.partial_penetration_frac_rect_lapl(S, xd, yd, zd, xwd, ywd, zwd, xed, yed, zed, hwd,
                                                          xbound, ybound, zbound_up, zbound_down, growth_velocity_d,
                                                          calc_type, self.cfd, Fcd_shape, num_segments, 10000000, self.sf, self.S_choke)
        else:
            R_in_d = self.R_in / unit_length
            p_d = self.PD_Frac_Rect_Lapl_composite(S, xd, yd, xwd, ywd, xed, yed, xbound, ybound, self.cfd,
                                                   R_in_d, self.rel_M,
                                                   self.rel_D, 10000000)

        return p_d

    def partial_penetration_frac_rect_lapl(self, S, xd, yd, zd, xwd, ywd, zwd, xed, yed, zed, hwd,
                                           xbound, ybound, zbound_up, zbound_down, growth_velocity_d,
                                           calc_type, cfd, Fcd_shape, num_segments,
                                           etad, sf=0, S_choke=0):
        """
        :param S: переменная пространства Лапласа
        :param xd: безразмерная координата в направлении x
        :param yd: безразмерная координата в направлении y
        :param zd: безразмерная координата в направлении z
        :param xwd: безразмерная координата скважины в направлении x
        :param ywd: безразмерная координата скважины в направлении y
        :param zwd: безразмерная координата скважины в направлении z
        :param xed: безразмерная координата границы x
        :param yed: безразмерная координата границы y
        :param zed: безразмерная координата границы z
        :param hwd: безразмерная величина вскрытия пласта
        :param xbound: тип границы по x
        :param ybound: тип границы по y
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z
        :param growth_velocity_d: безразмерная скорость роста трещины
        :param calc_type: тип расчёта
        :param cfd: безразмерный поправочный коэффициент на эффект ствола относительно полудлины трещины
        :param Fcd_shape: параметр, характеризующий производительность трещины
        :param num_segments: количество сегментов
        :param etad: безразмерный коэффициент диффузии трещины
        :param sf: скин на стенке трещины
        :param S_choke: скин за счет схождения
        :return: pwd
        """

        pwd = 0
        Ld = 1 / zed

        compl_type = "frac"

        if cfd < self.tiny_2:
            Uniform_Flux = True
        else:
            Uniform_Flux = False

        if Uniform_Flux:
            dxd = 0
        else:
            dxd = self.calcxd(cfd)

        if calc_type == "optimal":
            if zbound_up == "n" and zbound_down == "n":
                if S > 0:
                    if growth_velocity_d == 0:
                        pwd = (self.pd_lapl_rect(S, xwd + dxd, xwd, xed, yd, ywd, yed, xbound, ybound,
                                                 compl_type) +
                               self.sum_pp_1(S, xwd + dxd, xwd, xed, yd, ywd, yed, zd, zwd, zed, hwd, Ld,
                                             xbound,
                                             ybound,
                                             zbound_up, zbound_down, compl_type) +
                               self.dpd_fincond_sf(S, cfd, sf, etad, Uniform_Flux))
                    else:
                        pwd = self.pd_lapl_rect_frac_growth(S, xwd, xwd, xed, yd, ywd, yed, xbound, ybound,
                                                            growth_velocity_d)
                else:
                    pwd = (self.pd_rect_BD(xwd + dxd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type) +
                           self.sum_pp(0, xwd + dxd, xwd, xed, yd, ywd, yed, zd, zwd, zed, hwd, Ld, xbound,
                                       ybound,
                                       zbound_up, zbound_down, compl_type) +
                           self.dpd_fincond_sf_BD(cfd, sf, Uniform_Flux))

                if abs(xd - xwd) < self.tiny_2:
                    pwd += sf / (2 * pi)
            else:
                if S > 0:
                    pwd = self.sum_pp_1(S, xwd + dxd, xwd, xed, yd, ywd, yed, zd, zwd, zed, hwd, Ld, xbound,
                                        ybound,
                                        zbound_up, zbound_down, compl_type)
                else:
                    pwd = self.sum_pp(0, xwd + dxd, xwd, xed, yd, ywd, yed, zd, zwd, zed, hwd, Ld, xbound,
                                      ybound,
                                      zbound_up,
                                      zbound_down, compl_type)

        elif calc_type == "desuperposition":
            if zbound_up == "n" and zbound_down == "n":
                if S > 0:
                    pwd = (self.pd_lapl_rect(S, xwd + dxd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type) +
                           self.dpd_frac_conv(S, ywd, yed, zd, zwd, zed, hwd, ybound) +
                           self.dpd_fincond_sf(S, cfd, sf, etad, Uniform_Flux))
                else:
                    pwd = (self.pd_rect_BD(xwd + dxd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type) +
                           self.dpd_frac_conv(S, ywd, yed, zd, zwd, zed, hwd, ybound) +
                           self.dpd_fincond_sf_BD(cfd, sf, Uniform_Flux))
        elif calc_type == "segmentation":
            n_seg = num_segments // 2
            pwd = self.pd_cinco_frac_pp_lapl_rect(n_seg, S, zd, xwd, ywd, zwd, xed, yed, zed, hwd,
                                                  xbound, ybound, Fcd_shape, sf, S_choke, 1)

        return pwd

    def pd_cinco_frac_pp_lapl_rect(self, N, S, zd, xwd, ywd, zwd, xed, yed, zed, hwd, xbound, ybound, Fcd_shape,
                                   sf, S_choke,
                                   Qd=1):
        """
        :param N: размерность вектора
        :param S: переменная пространства Лапласа
        :param zd: безразмерная координата в направлении z
        :param xwd: безразмерная координата скважины в направлении x
        :param ywd: безразмерная координата скважины в направлении y
        :param zwd: безразмерная координата скважины в направлении z
        :param xed: безразмерная координата границы x
        :param yed: безразмерная координата границы y
        :param zed: безразмерная координата границы z
        :param hwd: безразмерная величина вскрытия пласта
        :param xbound: тип границы по x
        :param ybound: тип границы по y
        :param Fcd_shape: параметр, характеризующий производительность трещины
        :param sf: скин на стенке трещины
        :param S_choke: скин за счет схождения
        :param Qd: общий дебит скважины
        :return: ql
        """
        ql = self.cinco_frac_pp_lapl_rect(N, S, zd, xwd, ywd, zwd, xed, yed, zed, hwd, xbound, ybound,
                                          Fcd_shape,
                                          sf, S_choke, Qd)

        return ql[-1]

    def cinco_frac_pp_lapl_rect(self, N, S, zd, xwd, ywd, zwd, xed, yed, zed, hwd, xbound, ybound, Fcd_shape,
                                sf, S_choke,
                                Qd=1):
        """
        :param N: размерность вектора
        :param S: переменная пространства Лапласа
        :param zd: безразмерная координата в направлении z
        :param xwd: безразмерная координата скважины в направлении x
        :param ywd: безразмерная координата скважины в направлении y
        :param zwd: безразмерная координата скважины в направлении z
        :param xed: безразмерная координата границы x
        :param yed: безразмерная координата границы y
        :param zed: безразмерная координата границы z
        :param hwd: безразмерная величина вскрытия пласта
        :param xbound: тип границы по x
        :param ybound: тип границы по y
        :param Fcd_shape: параметр, характеризующий производительность трещины
        :param sf: скин на стенке трещины
        :param S_choke: скин за счет схождения
        :param Qd: Общий дебит скважины
        :return: rght
        """
        flow_res_matr = np.zeros((2 * N + 1, 2 * N + 1))
        infl_matr = np.zeros((2 * N + 1, 2 * N + 1))
        main_matr = np.zeros((2 * N + 2, 2 * N + 2))

        flow_res_matr = self.flow_resistance_asymmetrical_pp(flow_res_matr, N, Fcd_shape)
        rght = self.rhs(2 * N, Qd)
        infl_matr = self.influence_matrix_frac_pp_rect(S, infl_matr, N, zd, xwd, ywd, zwd, xed, yed, zed, hwd, xbound,
                                                       ybound,
                                                       sf, S_choke)
        main_matr = self.main_matrix(S, main_matr, flow_res_matr, infl_matr, 2 * N)
        rght = self.LU_(main_matr, rght, 2 * N + 1)

        return rght

    def dpd_fincond_sf(self, S, cfd, sf, eta, Uniform_Flux):
        """
        :param S: переменная пространства Лапласа
        :param cfd: безразмерный поправочный коэффициент на эффект ствола относительно полудлины трещины
        :param sf: скин на стенке трещины
        :param eta: безразмерный коэффициент диффузии гидроразрыва
        :param Uniform_Flux: Модель равномерного притока к скважине (булева операция)
        :return: dpd
        """

        mult = 100000

        if Uniform_Flux:
            dpd = sf / (2 * pi)
        else:
            dpd = self.p_tril_infinit_res(S, cfd, sf, eta) - self.p_tril_infinit_res(S, mult * cfd, 0 * sf,
                                                                                     mult * eta)
        return dpd

    def dpd_fincond_sf_BD(self, cfd, sf, Uniform_Flux):
        """
        :param cfd: безразмерный поправочный коэффициент на эффект ствола относительно полудлины трещины
        :param sf: скин на стенке трещины
        :param Uniform_Flux: Модель равномерного притока к скважине (булева операция)
        :return: dpd
        """

        if Uniform_Flux:
            dpd = sf / (2 * pi)
        else:
            dpd = (1 / (6 * cfd) + sf / (2 * pi))

        return dpd

    def pd_lapl_rect_frac_growth(self, S, xd, xwd, xed, yd, ywd, yed, xbound, ybound, alpha):
        """
        :param S: переменная пространства Лапласа
        :param xd: безразмерная координата в направлении x
        :param xwd: безразмерная координата скважины в направлении x
        :param xed: безразмерная координата границы x
        :param yd: безразмерная координата в направлении y
        :param ywd: безразмерная координата скважины в направлении y
        :param yed: безразмерная координата границы y
        :param xbound: тип границы по x
        :param ybound: тип границы по y
        :param alpha: коэффициент (перфорированная часть ствола)
        :return: pd_lapl_rect_frac_growth
        """

        pd_lapl_rect_frac_growth = self.pd_b1_frac_growth(S, xwd, yd, ywd, yed, xed, xbound, ybound, alpha, False) + \
                                   self.pd_b2_frac_growth(S, xd, xwd, xed, yd, ywd, yed, xbound, ybound, alpha, False)

        return pd_lapl_rect_frac_growth

    def pd_b1_frac_growth(self, S, xwd, yd, ywd, yed, xed, xbound, ybound, alpha, subtract_inf=True):
        """
        :param S: переменная пространства Лапласа
        :param xwd: безразмерная координата скважины в направлении x
        :param yd: безразмерная координата в направлении y
        :param ywd: безразмерная координата скважины в направлении y
        :param yed: безразмерная координата границы y
        :param xed: безразмерная координата границы x
        :param xbound: тип границы по x
        :param ybound: тип границы по y
        :param alpha: коэффициент (перфорированная часть ствола)
        :param subtract_inf: параметр, определяющий, нужно ли считать часть уравнения
        :return: pd
        """

        if xbound == "c":
            return 0

        if ybound == "n":
            signum = 1
        elif ybound == "c":
            signum = -1
        else:
            signum = 1

        if subtract_inf:
            sbtr_inf = 0
        else:
            sbtr_inf = 1

        if S < self.small / yed ** 2:
            if ybound == "n":
                pd = yed / xed * (
                        1 / 3 - 1 / 2 * (fabs(yd - ywd) + (yd + ywd)) / yed + (yd ** 2 + ywd ** 2) / (
                        2 * yed ** 2))
                if S > 0:
                    pd = pd + 1 / (S * yed * xed)
            else:
                pd = yed / xed * (1 / 2 * ((yd + ywd) - fabs(yd - ywd)) / yed - (yd * ywd) / yed ** 2)
        else:
            u = S ** 2
            sum_ = self.sumexp(u, yed)
            yd1 = yed - fabs(yd - ywd)
            yd2 = yed - (yd + ywd)
            pd = 1 / 2 / u / xed * (exp(-u * fabs(yd - ywd)) * (sbtr_inf + sum_) +
                                    (signum * exp(-u * (yd + ywd)) + signum * exp(-u * (yed + yd2)) + exp(
                                        -u * (yed + yd1))) * (1 + sum_)) * \
                 (exp(S / alpha * (1 - xwd)) - exp(-S * (1 + xwd) / alpha)) / S

        return pd

    def pd_b2_frac_growth(self, S, xd, xwd, xed, yd, ywd, yed, xbound, ybound, alpha, subtract_inf=True):
        """
        :param S: переменная пространства Лапласа
        :param xd: безразмерная координата в направлении x
        :param xwd: безразмерная координата скважины в направлении x
        :param xed: безразмерная координата границы x
        :param yd: безразмерная координата в направлении y
        :param ywd: безразмерная координата скважины в направлении y
        :param yed: безразмерная координата границы y
        :param xbound: тип границы по x
        :param ybound: тип границы по y
        :param alpha: коэффициент (перфорированная часть ствола)
        :param subtract_inf: параметр, определяющий, нужно ли считать часть уравнения
        :return: sum_val
        """

        sum_val = 0
        k = 1
        while True:
            psum = 0
            psumabs = 0
            for i in range(1, self.part_sum_num + 1):
                add = self.pd_b2_k_frac_growth(k, S, xd, xwd, xed, yd, ywd, yed, xbound, ybound, alpha, subtract_inf)
                psum += add
                psumabs += fabs(add)
                k += 1
            sum_val += psum
            if not abs(psumabs) / self.part_sum_num >= self.tiny * (fabs(sum_val) + self.tiny) and k < self.max_it:
                break

        return sum_val

    def pd_b2_k_frac_growth(self, k, S, xd, xwd, xed, yd, ywd, yed, xbound, ybound, alpha, subtract_inf):
        """
        :param k: проницаемость пласта, Д
        :param S: переменная пространства Лапласа
        :param xd: безразмерная координата в направлении x
        :param xwd: безразмерная координата скважины в направлении x
        :param xed: безразмерная координата границы x
        :param yd: безразмерная координата в направлении y
        :param ywd: безразмерная координата скважины в направлении y
        :param yed: безразмерная координата границы y
        :param xbound: тип границы по x
        :param ybound: тип границы по y
        :param alpha: коэффициент (перфорированная часть ствола)
        :param subtract_inf: параметр, определяющий, нужно ли считать часть уравнения
        :return: part__1 * part__2
        """

        if xbound == "n":
            part__1 = 2 * xed * alpha * cos(k * pi * xwd / xed) * (S * xed + exp(-S / alpha) *
                                                                   (k * pi * alpha * sin(
                                                                       k * pi / xed) - S * xed * cos(
                                                                       k * pi / xed))) / \
                      ((S * xed) ** 2 + (k * pi * alpha) ** 2) * cos(k * pi * xd / xed)
        elif xbound == "c":
            part__1 = 2 * xed * alpha * sin(k * pi * xwd / xed) * (S * xed + exp(-S / alpha) *
                                                                   (k * pi * alpha * sin(
                                                                       k * pi / xed) - S * xed * cos(

                                                                       k * pi / xed))) / \
                      ((S * xed) ** 2 + (k * pi * alpha) ** 2) * sin(k * pi * xd / xed)
        else:
            part__1 = 0

        if ybound == "n":
            signum = 1
        elif ybound == "c":
            signum = -1
        else:
            signum = 1

        if subtract_inf:
            sbtr_inf = 0
        else:
            sbtr_inf = 1

        ek = sqrt(S + k ** 2 * pi ** 2 / xed ** 2)

        smexp = self.sumexp(ek, yed)

        part__2 = exp(-ek * fabs(yd - ywd)) * (sbtr_inf + smexp) + \
                  (signum * exp(-ek * (yd + ywd)) + exp(-ek * (2 * yed - fabs(yd - ywd))) +
                   signum * exp(-ek * (2 * yed - (yd + ywd)))) * (1 + smexp)

        part__2 = 1 / ek * part__2 * 1 / (2 * xed)

        return part__1 * part__2

    def influence_matrix_frac_pp_rect(self, S, nn, N, zd, xwd, ywd, zwd, xed, yed, zed, hwd, xbound, ybound,
                                      sf,
                                      S_choke):
        """
        :param S: переменная пространства Лапласа
        :param nn: матрица влияния
        :param N: размерность вектора
        :param zd: безразмерная координата в направлении z
        :param xwd: безразмерная координата скважины в направлении x
        :param ywd: безразмерная координата скважины в направлении y
        :param zwd: безразмерная координата скважины в направлении z
        :param xed: безразмерная координата границы x
        :param yed: безразмерная координата границы y
        :param zed: безразмерная координата границы z
        :param hwd: безразмерная величина вскрытия пласта
        :param xbound: тип границы по x
        :param ybound: тип границы по y
        :param sf: скин на стенке трещины
        :return: nn
        """
        dx = 1 / N
        unit_length = dx / 2
        Sd = S * unit_length ** 2
        xedd = xed / unit_length
        yedd = yed / unit_length
        ywdd = ywd / unit_length
        S_frac_dd = sf / unit_length

        dpd = self.dpd_frac_conv(S, ywd, yed, zd, zwd, zed, hwd, ybound)

        ydd = ywdd

        for i in range(1, 2 * N + 1):
            xdd = (xwd + (-1 + dx / 2 + dx * (i - 1))) / unit_length

            for j in range(1, 2 * N + 1):
                xwdd = (xwd + (-1 + dx / 2 + dx * (j - 1))) / unit_length

                if Sd > 0:
                    nn[i, j] = self.pd_lapl_rect(Sd, xdd, xwdd, xedd, ydd, ywdd, yedd, xbound, ybound,
                                                 "frac") + dpd + S_choke / (pi * 2)
                else:
                    nn[i, j] = self.pd_rect_BD(xdd, xwdd, xedd, ydd, ywdd, yedd, xbound, ybound,
                                               "frac") + dpd + S_choke / (pi * 2)

                if i == j:
                    nn[i, j] += S_frac_dd / (pi * 2)

        return nn

    def flow_resistance_asymmetrical_pp(self, mm, N, Fcd_shape):
        """
        :param mm: матрица сопротивлений
        :param N: размерность вектора
        :param Fcd_shape: параметр, характеризующий производительность трещины
        :return: mm
        """

        dx = 1 / N
        a = np.zeros(2 * N + 1)

        for i in range(1, 2 * N + 1):
            a[i] = dx / Fcd_shape[i]

        for i in range(1, N + 1):
            for k in range(1, i):
                mm[i, k] = a[k] * (N - (i - 1 / 2))

            mm[i, i] = a[i] * (N - (i - 3 / 8))

            for k in range(i + 1, N + 1):
                mm[i, k] = a[k] * (N - (k - 1 / 2))

        for i in range(1, N + 1):
            for k in range(1, i):
                mm[i + N, k + N] = a[k + N] * (k - 1 / 2)

            mm[i + N, i + N] = a[i + N] * (i - 5 / 8)

            for k in range(i + 1, N + 1):
                mm[i + N, k + N] = a[k + N] * (i - 1 / 2)

        return mm

    def dpd_frac_conv(self, S, ywd, yed, zd, zwd, zed, hwd, ybound):
        """
        :param S: переменная пространства Лапласа
        :param ywd: безразмерная координата скважины в направлении y
        :param yed: безразмерная координата границы y
        :param zd: безразмерная координата в направлении z
        :param zwd: безразмерная координата скважины в направлении z
        :param zed: безразмерная координата границы z
        :param hwd: безразмерная величина вскрытия пласта
        :param ybound: тип границы по y
        :return: dpd_frac_conv
        """

        compl_type_frac = "frac"
        vert_bound = "n"

        pd_mult = zed / 2

        Increased_speed = False

        yd = ywd
        zd = zwd

        unit_length = hwd / 2

        ydd = yd / unit_length
        zdd = zd / unit_length

        ywdd = ywd / unit_length
        zwdd = zwd / unit_length

        yedd = yed / unit_length
        zedd = zed / unit_length

        Sdd = S * unit_length ** 2

        if Sdd > 0:
            if Increased_speed:
                dpd = self.pd_lapl_rect(Sdd, ydd, ywdd, yedd, zdd, zwdd, zedd, ybound, vert_bound,
                                        compl_type_frac)
            else:
                dpd = self.pd_lapl_rect(Sdd, zdd, zwdd, zedd, ydd, ywdd, yedd, vert_bound, ybound,
                                        compl_type_frac)
        else:
            if Increased_speed:
                dpd = self.pd_rect_BD(ydd, ywdd, yedd, zdd, zwdd, zedd, ybound, vert_bound, compl_type_frac)
            else:
                dpd = self.pd_rect_BD(zdd, zwdd, zedd, ydd, ywdd, yedd, vert_bound, ybound, compl_type_frac)

        unit_length = zed / 2

        ydd = yd / unit_length
        zdd = zd / unit_length

        ywdd = ywd / unit_length
        zwdd = zwd / unit_length

        yedd = yed / unit_length
        zedd = zed / unit_length

        Sdd = S * unit_length ** 2

        if Sdd > 0:
            dpd = dpd - self.pd_lapl_rect(Sdd, zdd, zwdd, zedd, ydd, ywdd, yedd, vert_bound, ybound,
                                          compl_type_frac)
        else:
            dpd = dpd - self.pd_rect_BD(zdd, zwdd, zedd, ydd, ywdd, yedd, vert_bound, ybound, compl_type_frac)

        dpd_frac_conv = dpd * pd_mult

        return dpd_frac_conv

    def PD_Frac_Rect_Lapl_composite(self, S, xd, yd, xwd, ywd, xed, yed, xbound, ybound, cfd, R_in_d, rel_M,
                                    rel_D,
                                    etad):
        """
        :param S: переменная пространства Лапласа
        :param xd: безразмерная координата в направлении x
        :param yd: безразмерная координата в направлении y
        :param xwd: безразмерная координата скважины в направлении x
        :param ywd: безразмерная координата скважины в направлении y
        :param xed: безразмерная координата границы x
        :param yed: безразмерная координата границы y
        :param xbound: тип границы по x
        :param ybound: тип границы по y
        :param cfd: безразмерный поправочный коэффициент на эффект ствола относительно полудлины трещины
        :param R_in_d: безразмерный внутренний диаметр
        :param rel_M: коэффициент подвижности
        :param rel_D: коэффициент пьезопроводности
        :param etad: безразмерный коэффициент диффузии трещины
        :return: pwd
        """
        compl_type = "frac"

        if cfd < self.tiny:
            Uniform_Flux = True
        else:
            Uniform_Flux = False

        if Uniform_Flux:
            dxd = 0
        else:
            dxd = self.calcxd(cfd)

        if (abs(xd - xwd) < 1) and (abs(yd - ywd) < self.tiny):
            if S > 0:
                pwd = (self.pd_lapl_rect(S, xwd + dxd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type) +
                       self.dpd_fincond_composite(S, cfd, None, etad, Uniform_Flux, R_in_d, rel_M,
                                                  rel_D))
            else:
                pwd = 0
        else:
            if S > 0:
                pwd = self.pd_lapl_rect(S, xd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type)
            else:
                pwd = self.pd_rect_BD(xd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type)

        return pwd

    def p_tril_infinit_res(self, S, cfd, sf, eta):
        """
        :param S: переменная пространства Лапласа
        :param cfd: безразмерный поправочный коэффициент на эффект ствола относительно полудлины трещины
        :param sf: скин на стенке трещины
        :param eta: безразмерный коэффициент диффузии гидроразрыва
        :return: p_tril_infinit_res
        """
        u = (S + S ** 0.5) ** 0.5
        Psi = (1 * (S / eta) + 2 / cfd * u / (1 + 2 / pi * sf * u)) ** 0.5

        p_tril_infinit_res = 1 / (2 * cfd) / (Psi * tanh(Psi))

        return p_tril_infinit_res

    def p_tril_infinit_res_composite(self, S, cfd, eta, R_in_d, rel_M, rel_D):
        """
        :param S: переменная пространства Лапласа
        :param cfd: безразмерный поправочный коэффициент на эффект ствола относительно полудлины трещины
        :param eta: безразмерный коэффициент диффузии гидроразрыва
        :param R_in_d: безразмерный внутренний диаметр
        :param rel_M: коэффициент подвижности
        :param rel_D: коэффициент пьезопроводности
        :return: p_tril_infinit_res_composite
        """

        alpha_o_1 = (S ** 0.5) + S
        alpha_o_2 = (S ** 0.5) / rel_M + S / rel_D
        gamma_o = (alpha_o_1 ** 0.5) / rel_M

        beta_f_arg = exp(-2 * (alpha_o_2 ** 0.5) * R_in_d)
        beta_f = (alpha_o_2 * (1 - beta_f_arg) + gamma_o * (alpha_o_2 ** 0.5) * (1 + beta_f_arg)) / \
                 ((alpha_o_2 ** 0.5) * (1 + beta_f_arg) + gamma_o * (1 - beta_f_arg))

        Psi = (S / eta + 2 / cfd * rel_M * beta_f) ** 0.5

        p_tril_infinit_res_composite = 1 / (2 * cfd) / (Psi * tanh(Psi))

        return p_tril_infinit_res_composite

    def dpd_fincond_composite(self, S, cfd, sf, eta, Uniform_Flux, R_in_d, rel_M, rel_D):
        """
        :param S: переменная пространства Лапласа
        :param cfd: безразмерный поправочный коэффициент на эффект ствола относительно полудлины трещины
        :param sf: скин на стенке трещины
        :param eta: безразмерный коэффициент диффузии гидроразрыва
        :param Uniform_Flux: Модель равномерного притока к скважине (булева операция)
        :param R_in_d: безразмерный внутренний диаметр
        :param rel_M: коэффициент подвижности
        :param rel_D: коэффициент пьезопроводности
        :return: dpd
        """
        mult = 100000

        if Uniform_Flux:
            dpd = 0
        else:
            dpd = self.p_tril_infinit_res_composite(S, cfd, eta, R_in_d, rel_M, rel_D) - self.p_tril_infinit_res(S,
                                                                                                                     mult * cfd,
                                                                                                                     0,
                                                                                                                     mult * eta)

        return dpd

    def fracture_shape_width(self, n_seg, xf, wellbore_wf, res_wf, grade):
        """
        :param n_seg: количество сегментов
        :param xf: полудлина трещины, м
        :param wellbore_wf: ширина трещины у ствола скважины, мм
        :param res_wf: ширина трещины на конце, мм
        :param grade: степень аппроксимационной функции
        :return: arr_width
        """
        n_segments_profile_max = 100
        arr_width = np.zeros(shape=n_segments_profile_max)
        arr_coord = self.coordinates_of_segments(xf, n_seg)
        if wellbore_wf == res_wf:
            for i in range(n_seg // 2):
                arr_width[i] = wellbore_wf
                arr_width[i + n_seg // 2] = wellbore_wf
        else:
            a = xf * 2 ** grade / (res_wf ** grade - wellbore_wf ** grade)
            b = -a * (wellbore_wf / 2) ** grade
            for i in range(n_seg // 2):
                arr_width[i] = ((arr_coord[i] + b) / -a) ** (1 / grade) * 2
                arr_width[i + n_seg // 2] = ((arr_coord[i + n_seg // 2] - b) / a) ** (1 / grade) * 2
        for i in range(n_seg, n_segments_profile_max):
            arr_width[i] = arr_width[n_seg - 1]

        return arr_width

    def coordinates_of_segments(self, x_f, n_seg):
        """
        :param x_f: полудлина трещины, м
        :param n_seg: количество сегментов
        :return: arr_coord
        """
        n_segments_profile_max = 100

        arr_coord = np.zeros(shape=n_segments_profile_max)
        dx = x_f * 2 / n_seg
        arr_coord[0] = -x_f + dx / 2

        for i in range(1, n_seg):
            arr_coord[i] = arr_coord[n_seg - 1] + dx
        for i in range(n_seg, n_segments_profile_max):
            arr_coord[i] = arr_coord[n_seg - 1]

        return arr_coord
