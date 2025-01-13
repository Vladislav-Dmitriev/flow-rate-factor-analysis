import numpy as np
from math import exp, sqrt, cos, sin, log, pi

from scipy import special
from src.entities import series_fn


class LaplaceCalculator:
    def __init__(self, extend_reflections_mode, tiny, tiny_2, large_s, small,
                 part_sum_num, max_it, small_bessel_arg, part_sum_num_b2_1, part_sum_num_f_4):
        """
        :param extend_reflections_mode: активировать расширенный режим отражений
        :param tiny: параметр управления сходимостью рядов
        :param large_s: параметр улучшения сходимости
        :param small: параметр аппроксимации времени
        :param part_sum_num: минимальное приращение
        :param max_it: максимальное количество итераций при суммировании рядов
        :param small_bessel_arg: промежуток времени для функции Бесселя
        """
        self.gamma = 0.577215664901533
        self.extend_reflections_mode = extend_reflections_mode
        self.tiny = tiny
        self.tiny_2 = tiny_2
        self.large_s = large_s
        self.small = small
        self.part_sum_num = part_sum_num
        self.part_sum_num_b1 = part_sum_num
        self.part_sum_num_b3 = part_sum_num
        self.part_sum_num_f_3 = part_sum_num
        self.max_it = max_it
        self.small_bessel_arg = small_bessel_arg
        self.part_sum_num_b2_1 = part_sum_num_b2_1
        self.part_sum_num_b2_2 = part_sum_num_b2_1
        self.part_sum_num_f_4 = part_sum_num_f_4

    def coef(self, N):
        """
        :param N: количество коэффициентов Стефеста
        :return: коэффициенты Стефеста
        """
        v = np.zeros(N + 1)
        g = np.zeros(N + 1)
        h = np.zeros(N + 1)
        m = 0
        if m != N:
            g[1] = 1
            NH = int(N / 2)
            for i in range(2, N + 1):
                g[i] = g[i - 1] * i
            h[1] = 2 / g[NH - 1]
            for i in range(2, NH + 1):
                if i != NH:
                    h[i] = (i ** NH) * g[2 * i] / (g[NH - i] * g[i] * g[i - 1])
                else:
                    h[i] = (i ** NH) * g[2 * i] / (g[i] * g[i - 1])
            SN = 2 * (NH - (NH // 2) * 2) - 1
            for i in range(1, N + 1):
                v[i] = 0
                K1 = (i + 1) // 2
                K2 = i
                if K2 > NH:
                    K2 = NH
                for k in range(K1, K2 + 1):
                    if 2 * k - i == 0:
                        v[i] = v[i] + h[k] / (g[i - k])
                        continue
                    if i == k:
                        v[i] = v[i] + h[k] / g[2 * k - i]
                        continue
                    v[i] = v[i] + h[k] / (g[i - k] * g[2 * k - i])
                v[i] = SN * v[i]
                SN = -SN
        return v

    def calc_distance(self, x0, y0, x1, y1):
        """
        Рассчитывает расстояние между двумя точками
        :param x0: x координата точки
        :param y0: y координата точки
        :param x1: x координата второй точки
        :param y1: y координата второй точки
        :return: расстояние
        """
        distance = sqrt((x0 - x1) ** 2 + (y0 - y1) ** 2)
        return distance

    def calc_pd(self, S, xd, yd, xwd, ywd, xed, yed, xbound, ybound, skin):
        """
        Рассчитывает безразмерное давление для единичной скорости Лапласа
        :param S: переменная пространства Лапласа
        :param xd: безразмерная координата точки, в которой рассчитывается давление
        :param yd: безразмерная координата точки, в которой рассчитывается давление
        :param xwd: безразмерная координата скважины
        :param ywd: безразмерная координата скважины
        :param xed: безразмерная координата границы x
        :param yed: безразмерная координата границы y
        :param xbound: тип границы
        :param ybound: тип границы
        :param skin: скин-фактор
        :return: безразмерное давление для единичной скорости Лапласа
        """
        pd = self.calc_pwd(S, xd, yd, xwd, ywd, xed, yed, xbound, ybound, skin)
        in_wellbore = self.calc_distance(xd, yd, xwd, ywd) < (1 + self.tiny)
        if in_wellbore and pd > 0:
            pd = pd + skin / (2 * pi)
        return pd

    def calc_pwd(self, S, xd, yd, xwd, ywd, xed, yed, xbound, ybound, skin):
        """
        :param S: переменная пространства Лапласа
        :param xd: безразмерная координата точки, в которой рассчитывается давление
        :param yd: безразмерная координата точки, в которой рассчитывается давление
        :param xwd: безразмерная координата скважины
        :param ywd: безразмерная координата скважины
        :param xed: безразмерная координата границы x
        :param yed: безразмерная координата границы y
        :param xbound: тип границы
        :param ybound: тип границы
        :param skin: скин-фактор
        :return: давление
        """
        rho = sqrt((xd - xwd) ** 2 + (yd - ywd) ** 2)
        in_wellbore = rho <= 1 + self.tiny
        if in_wellbore:
            xd1 = xwd + 1
            yd1 = ywd
        else:
            xd1 = xd
            yd1 = yd
        compl_type = 'vert'
        if S > 0:
            pwd = self.pd_lapl_rect(S, xd1, xwd, xed, yd1, ywd, yed, xbound, ybound, compl_type)
        else:
            pwd = self.pd_rect_BD(xd1, xwd, xed, yd1, ywd, yed, xbound, ybound, compl_type)
        if in_wellbore:
            pwd = pwd + skin / (2 * pi)
        return pwd

    def pd_rect_BD(self, xd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type):
        """
        :param xd: безразмерная координата точки, в которой рассчитывается давление
        :param xwd: безразмерная координата скважины
        :param xed: безразмерная координата границы x
        :param yd: безразмерная координата точки, в которой рассчитывается давление
        :param ywd: безразмерная координата скважины
        :param yed: безразмерная координата границы y
        :param xbound: тип границы
        :param ybound: тип границы
        :param compl_type: тип скважины
        :return: давление
        """
        if abs(xed - 2) < self.tiny:
            pwd = self.pd_b1(-1, yd, ywd, yed, xed, xbound, ybound, False)
            return pwd
        else:
            pwd = (
                    self.pd_b1(-1, yd, ywd, yed, xed, xbound, ybound, True)
                    + self.pd_b2_2(
                0,
                xd,
                xwd,
                xed,
                yd,
                ywd,
                yed,
                xbound,
                ybound,
                compl_type,
                True,
            )
                    + self.pd_b3_BD(xd, xwd, xed, yd, ywd, xbound, compl_type)
            )
            return pwd

    def pd_b3_BD(self, xd, xwd, xed, yd, ywd, xbound, compl_type):
        """
        :param xd: безразмерная координата точки, в которой рассчитывается давление
        :param xwd: безразмерная координата скважины
        :param xed: безразмерная координата границы x
        :param yd: безразмерная координата точки, в которой рассчитывается давление
        :param ywd: безразмерная координата скважины
        :param xbound: тип границы
        :param compl_type: тип скважины
        :return: давление
        """
        signum = 1
        if xbound == "c":
            signum = -1
        X = pi * xd / xed
        Y = pi * xwd / xed
        t = pi * abs(yd - ywd) / xed
        if compl_type == 'vert':
            dpd1 = log(4 * ((sin((X - Y) / 2)) ** 2 + (np.sinh(t / 2)) ** 2))
            dpd2 = log(4 * ((sin((X + Y) / 2)) ** 2 + (np.sinh(t / 2)) ** 2))
            dpd = -1 / (2 * pi) * 1 / 2 * (dpd1 + signum * dpd2 - (1 + signum) * t)
        else:
            sin2a = 1 / (np.cosh(t / 2)) ** 2
            a = pi / 2 * (1 - 1 / xed)
            b = pi / 2 * (1 + 1 / xed)
            alpha = (X - Y) / 2
            betta = (X + Y) / 2
            dpd = (1 + signum) * log(4 / sin2a)
            dpd1 = self.L_m(b - alpha, sin2a) - self.L_m(a - alpha, sin2a)
            dpd2 = self.L_m(b - betta, sin2a) - self.L_m(a - betta, sin2a)
            dpd = -1 / (2 * pi) * 1 / 2 * (dpd + xed / pi * (dpd1 + signum * dpd2) - (1 + signum) * t)

        return dpd

    def L_m(self, u, sin2a, sum_number=100):
        """
        :param u: f(x)
        :param sin2a: синус двойного угла
        :param sum_number: число элементов суммы
        :return: интеграл
        """
        tiny = 1E-20
        if sin2a < tiny or sin2a > 1:
            L_m = 0
            return L_m
        signum = np.signum(u)
        u = abs(u)
        if (1 - sin2a) < tiny:
            integral = -2 * self.l(u, sum_number)
        else:
            k = int(u / pi)
            u = u - k * pi
            sina = sqrt(sin2a)
            cosa = sqrt(1 - sin2a)
            ctga2 = sina / (1 - cosa)
            if abs(u - pi / 2) < tiny:
                Theta = 0
            elif u < tiny:
                Theta = pi / 2
            else:
                Theta = np.arctan(1 / (cosa * np.tan(u)))
            integral = 2 * pi * k * log((cosa + 1) / 2)
            integral = integral + (pi - 2 * Theta) * log(ctga2) + 2 * u * log(1 / 2 * sina) - pi / 2 * log(2) + self.l(
                Theta + u, sum_number) - self.l(Theta - u, sum_number) + self.l(pi / 2 - 2 * u, sum_number)
        L_m = signum * integral

        return L_m

    def l(self, X, sum_number=100):
        """
        :param X: аргумент
        :param sum_number: число элементов суммы
        :return: интеграл
        """
        summa = 0
        for k in range(1, sum_number + 1):
            summa = summa + (-1) ** (k - 1) * sin(2 * k * X) / k ** 2
        summa = X * log(2) - 1 / 2 * summa
        l = summa

        return l

    def pd_b1(self, S, yd, ywd, yed, xed, xbound, ybound, subtract_inf=False):
        """
        :param S: переменная пространства Лапласа
        :param yd: безразмерная координата точки, в которой рассчитывается давление
        :param ywd: безразмерная координата скважины
        :param yed: безразмерная координата границы y
        :param xed: безразмерная координата границы x
        :param xbound: тип границы
        :param ybound: тип границы
        :param subtract_inf: параметр, определяющий, нужно ли считать часть уравнения
        :return: давление
        """
        if xbound == "c":
            return 0
        signum = 1
        if ybound == "c":
            signum = -1
        sbtr_inf = 1
        if subtract_inf:
            sbtr_inf = 0
        sbtr_inf_ext = 1
        if self.extend_reflections_mode and subtract_inf:
            sbtr_inf_ext = 0
        if S < self.small / yed ** 2:
            if ybound == "n":
                pd = (
                        yed
                        / xed
                        * (
                                1 / 3
                                - 1 / 2 * (abs(yd - ywd) + yd + ywd) / yed
                                + (yd ** 2 + ywd ** 2) / (2 * yed ** 2)
                        )
                )
                if S > 0:
                    pd = pd + 1 / (S * yed * xed)
                return pd
            else:
                pd = (
                        yed
                        / xed
                        * (1 / 2 * ((yd + ywd) - abs(yd - ywd)) / yed - (yd * ywd) / yed ** 2)
                )
                return pd
        else:
            u = sqrt(S)
            summa = self.sumexp(u, yed)
            pd = (
                    1
                    / 2
                    / u
                    / xed
                    * (
                            exp(-u * abs(yd - ywd)) * sbtr_inf
                            + (
                                    signum * exp(-u * (yd + ywd))
                                    + signum * exp(-u * (2 * yed - (yd + ywd)))
                                    + exp(-u * (2 * yed - (ywd - yd)))
                                    + exp(-u * (2 * yed - (yd - ywd)))
                            )
                            * (sbtr_inf_ext + summa)
                    )
            )
            return pd

    def pd_b3(self, S, xd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type):
        """
        :param S: переменная пространства Лапласа
        :param xd: безразмерная координата точки, в которой рассчитывается давление
        :param xwd:безразмерная координата скважины
        :param xed: безразмерная координата границы x
        :param yd: безразмерная координата точки, в которой рассчитывается давление
        :param ywd: безразмерная координата скважины
        :param yed: безразмерная координата границы y
        :param xbound: тип границы
        :param ybound: тип границы
        :param compl_type: тип скважины
        :return: давление
        """
        if self.extend_reflections_mode:
            signum = 1
            if ybound == "c":
                signum = -1
            pd = (
                    self.pd_b3_one_row(S, xd, xwd, xed, yd, ywd, xbound, compl_type)
                    + signum
                    * self.pd_b3_one_row(S, xd, xwd, xed, yd, 2 * yed - ywd, xbound, compl_type)
                    + signum * self.pd_b3_one_row(S, xd, xwd, xed, yd, -ywd, xbound, compl_type)
                    + self.pd_b3_one_row(S, xd, xwd, xed, yd, 2 * yed + ywd, xbound, compl_type)
                    + self.pd_b3_one_row(S, xd, xwd, xed, yd, -2 * yed + ywd, xbound, compl_type)
            )
            return pd
        else:
            pd = self.pd_b3_one_row(S, xd, xwd, xed, yd, ywd, xbound, compl_type)
            return pd

    def pd_lapl_rect(self, S, xd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type):
        """
        :param S: переменная пространства Лапласа
        :param xd: безразмерная координата точки, в которой рассчитывается давление
        :param xwd: безразмерная координата скважины
        :param xed: безразмерная координата границы x
        :param yd: безразмерная координата точки, в которой рассчитывается давление
        :param ywd: безразмерная координата скважины
        :param yed: безразмерная координата границы y
        :param xbound: тип границы
        :param ybound: тип границы
        :param compl_type: тип скважины
        :return: давление
        """
        if abs(xed - 2) < self.tiny_2:
            pd_lapl_rect = self.pd_b1(S, yd, ywd, yed, xed, xbound, ybound, False)
            return pd_lapl_rect
        else:
            if S > self.large_s * xed ** 2 and abs(yd - ywd) < self.tiny_2:
                pd_lapl_rect = (
                        self.pd_b1(S, yd, ywd, yed, xed, xbound, ybound, True)
                        + self.pd_b2_2(
                    S,
                    xd,
                    xwd,
                    xed,
                    yd,
                    ywd,
                    yed,
                    xbound,
                    ybound,
                    compl_type,
                    True,
                )
                        + self.pd_b3(S, xd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type)
                )
                return pd_lapl_rect
            else:
                pd_lapl_rect = self.pd_b1(
                    S, yd, ywd, yed, xed, xbound, ybound, False
                ) + self.pd_b2_2(
                    S,
                    xd,
                    xwd,
                    xed,
                    yd,
                    ywd,
                    yed,
                    xbound,
                    ybound,
                    compl_type,
                    False,
                )
                return pd_lapl_rect

    def pd_b3_one_row(self, S, xd, xwd, xed, yd, ywd, xbound, compl_type):
        """
        :param S: переменная пространства Лапласа
        :param xd: безразмерная координата точки, в которой рассчитывается давление
        :param xwd: безразмерная координата скважины
        :param xed: безразмерная координата границы x
        :param yd: безразмерная координата точки, в которой рассчитывается давление
        :param ywd: безразмерная координата скважины
        :param xbound: тип границы
        :param compl_type: тип скважины
        :return: давление
        """
        signum = 1
        if xbound == "c":
            signum = -1
        if compl_type == 'vert':
            xd1m00 = (xd - xwd) ** 2
            xd1p00 = (xd + xwd) ** 2
            yd1 = (yd - ywd) ** 2
            sum_ = (self.unit_cylinder_source(S, sqrt(xd1m00 + yd1)) +
                    signum * self.unit_cylinder_source(S, sqrt(xd1p00 + yd1)))
            sum_ = sum_ + series_fn.nsum(
                lambda k: (
                        signum
                        * self.unit_cylinder_source(S, sqrt((xd + xwd + 2 * xed * k) ** 2 + yd1))
                        + self.unit_cylinder_source(S, sqrt((xd - xwd + 2 * xed * k) ** 2 + yd1))
                        + signum
                        * self.unit_cylinder_source(S, sqrt((xd + xwd - 2 * xed * k) ** 2 + yd1))
                        + self.unit_cylinder_source(S, sqrt((xd - xwd - 2 * xed * k) ** 2 + yd1))
                ), 1, np.inf,
            )
        elif compl_type == "frac":
            yd1 = abs(yd - ywd)
            sum_ = (self.unit_fracture_func(S, abs(xd - xwd), yd1) +
                    signum * self.unit_fracture_func(S, abs(xd + xwd), yd1))
            sum_ = sum_ + series_fn.nsum(
                lambda k: (
                        signum
                        * self.unit_fracture_func(S, abs(xd + xwd + 2 * k * xed), yd1)
                        + self.unit_fracture_func(S, abs(xd - xwd + 2 * k * xed), yd1)
                        + signum * self.unit_fracture_func(S, abs(xd + xwd - 2 * k * xed), yd1)
                        + self.unit_fracture_func(S, abs(xd - xwd - 2 * k * xed), yd1)
                ), 1, np.inf,
            )
        else:
            sum_ = 0
        return sum_

    def unit_fracture_func(self, S, xd, yd):
        """
        :param S: переменная пространства Лапласа
        :param xd: безразмерная координата точки, в которой рассчитывается давление
        :param yd: безразмерная координата точки, в которой рассчитывается давление
        :return: pd однородной трещины
        """
        if S > 0:
            u = sqrt(S)
        else:
            u = -1
        if (sqrt(xd ** 2 + yd ** 2) + 1) * u / 2 > self.small_bessel_arg:
            if abs(yd) < self.tiny:
                abs_xd = abs(xd)
                if abs_xd <= 1:
                    lim_1 = u * (1 + xd)
                    lim_2 = u * (1 - xd)
                    sign = 1
                else:
                    lim_1 = u * (1 + abs_xd)
                    lim_2 = u * (abs_xd - 1)
                    sign = -1
                dpd = 1 / 2 * 1 / (2 * pi) * 1 / u * (self.fast_ik0(lim_1) + sign * self.fast_ik0(lim_2))
            else:
                dpd = 1 / 2 * 1 / (2 * pi) * self.int_k0_gauss_DW_(u, xd, yd)
        else:
            if u < 0:
                u = 1
            if abs(yd) > self.tiny:
                dpd = 1 / 4 * ((xd - 1) * log((xd - 1) ** 2 + yd ** 2) - (xd + 1) * log((xd + 1) ** 2 + yd ** 2)) + \
                      1 / 2 * yd * (np.arctan((xd - 1) / yd) - np.arctan((xd + 1) / yd)) + 1
            else:
                dpd = 1 / 4 * ((xd - 1) * log((xd - 1) ** 2) - (xd + 1) * log((xd + 1) ** 2)) + 1
            dpd = 1 / (2 * pi) * (dpd + log(2) - self.gamma - log(u))

        unit_fracture_func = dpd

        return unit_fracture_func

    def horizontal_rect_lapl(self, S, rwd, zed, yd, ywd, zd, zwd, xwd, xd, xed, yed, xbound, ybound, zbound_up,
                             zbound_down, Fcd, skin, perf, N_perf, num_segments, calc_type, j):
        """
        :param S: переменная пространства Лапласа
        :param rwd: безразмерный радиус скважины
        :param zed: безразмерная координата границы z
        :param yd: безразмерная координата точки, в которой рассчитывается давление
        :param ywd: безразмерная координата скважины
        :param zd: безразмерная координата точки, в которой рассчитывается давление
        :param zwd: безразмерная координата скважины
        :param xwd: безразмерная координата скважины
        :param xd: безразмерная координата точки, в которой рассчитывается давление
        :param xed: безразмерная координата границы x
        :param yed: безразмерная координата границы y
        :param xbound: тип границы
        :param ybound: тип границы
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z
        :param Fcd: безразмерная проводимость
        :param skin: скин-фактор
        :param perf: доля перфораций горизонтальной скважины
        :param N_perf: количество перфорированных интервалов
        :param num_segments: число сегментов
        :param calc_type: тип расчета
        :param j: номер итерации
        :return: давление
        """
        compl_type_hor = "frac"
        if Fcd > 0:
            dxd = 0.68
        else:
            dxd = 0
        if skin < 0 or skin > 0:
            rwd = rwd * exp(-skin)

        if calc_type == "optimal":
            if sqrt((yd - ywd) ** 2 + (zd - zwd) ** 2) < rwd:
                yd = ywd
                zd = zwd + rwd
            xd = xwd + dxd
            pwd = self.horizontal_well_for_cinco(S, xd, yd, zd, xwd, ywd, zwd, xed, yed, zed, xbound, ybound, zbound_up,
                                                 zbound_down)

        elif calc_type == "segmentation":
            if sqrt((yd - ywd) ** 2 + (zd - zwd) ** 2) < rwd:
                yd = ywd
                zd = zwd + rwd
            pwd = self.pd_cinco_hor_well_lapl_rect(S, num_segments, zd, xwd, ywd, zwd, xed, yed, zed,
                                                   rwd, Fcd, perf, N_perf, xbound, ybound, zbound_up, zbound_down)

        elif calc_type == "desuperposition":
            if abs(xd - xwd) < 1:
                if (j != 1) and (sqrt((yd - ywd) ** 2 + (zd - zwd) ** 2) < rwd):
                    yd = ywd + rwd
                    zd = zwd
                if S > 0:
                    pwd = self.pd_lapl_rect(S, xwd + dxd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type_hor)
                else:
                    pwd = self.pd_rect_BD(xwd + dxd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type_hor)
                if (j == 1) and (sqrt((yd - ywd) ** 2 + (zd - zwd) ** 2) < rwd):
                    yd = ywd + rwd
                    zd = zwd
                dpd_vert_conv = self.dpd_vert_conv(S, yd, ywd, yed, zd, zwd, zed, rwd, ybound)
                pwd += float(dpd_vert_conv)
            else:
                if S > 0:
                    pwd = self.pd_lapl_rect(S, xd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type_hor)
                else:
                    pwd = self.pd_rect_BD(xd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type_hor)

        horizontal_rect_lapl = pwd

        return horizontal_rect_lapl

    def dpd_vert_conv(self, S, yd, ywd, yed, zd, zwd, zed, rwd, ybound):
        """
        :param S: переменная пространства Лапласа
        :param yd: безразмерная координата точки, в которой рассчитывается давление
        :param ywd: безразмерная координата скважины
        :param yed: безразмерная координата границы y
        :param zd: безразмерная координата точки, в которой рассчитывается давление
        :param zwd: безразмерная координата скважины
        :param zed: безразмерная координата границы z
        :param rwd: безразмерный радиус скважины
        :param ybound: тип границы
        :return:
        """
        compl_type_vert = "vert"
        compl_type_frac = "frac"
        vert_bound = "n"
        pd_mult = zed / 2

        unit_length = rwd
        ydd = yd / unit_length
        zdd = zd / unit_length
        ywdd = ywd / unit_length
        zwdd = zwd / unit_length
        yedd = yed / unit_length
        zedd = zed / unit_length
        Sdd = S * unit_length ** 2

        if Sdd > 0:
            dpd = self.pd_lapl_rect(Sdd, ydd, ywdd, yedd, zdd, zwdd, zedd, ybound, vert_bound, compl_type_vert)
        else:
            dpd = self.pd_rect_BD(zdd, zwdd, zedd, ydd, ywdd, yedd, vert_bound, ybound, compl_type_vert)

        unit_length = zed / 2
        ydd = yd / unit_length
        zdd = zd / unit_length
        ywdd = ywd / unit_length
        zwdd = zwd / unit_length
        yedd = yed / unit_length
        zedd = zed / unit_length
        Sdd = S * unit_length ** 2
        if Sdd > 0:
            dpd = dpd - self.pd_lapl_rect(Sdd, zdd, zwdd, zedd, ydd, ywdd, yedd, vert_bound, ybound, compl_type_frac)
        else:
            dpd = dpd - self.pd_rect_BD(zdd, zwdd, zedd, ydd, ywdd, yedd, vert_bound, ybound, compl_type_frac)

        dpd_vert_conv = dpd * pd_mult

        return dpd_vert_conv

    def f(self, S, xd, xwd, zd, zwd, zed, Ld, yd, ywd, zbound_up, zbound_down):
        """
        :param S: переменная пространства Лапласа
        :param xd: безразмерная координата точки, в которой рассчитывается давление
        :param xwd: безразмерная координата скважины
        :param zd: безразмерная координата точки, в которой рассчитывается давление
        :param zwd: безразмерная координата скважины
        :param zed: безразмерная координата границы z
        :param Ld: безразмерная длина
        :param yd: безразмерная координата точки, в которой рассчитывается давление
        :param ywd: безразмерная координата скважины
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z
        :return:
        """
        if abs(xd - xwd) <= 1:
            if S <= 0:
                f = self.f_4_2(xd, xwd, zd, zwd, zed, Ld, zbound_up, zbound_down)
            else:
                f = 1 / 4 / pi * (
                        self.f_2(S, zd, zwd, zed, Ld, zbound_up, zbound_down) - self.f_3_2(S, xd, xwd, zd, zwd, zed, Ld,
                                                                                           zbound_up, zbound_down))
        else:
            if S <= 0:
                f = self.f_4_2(xd, xwd, zd, zwd, zed, Ld, zbound_up, zbound_down)
            else:
                f = 1 / 4 / pi * self.f_3_2(S, xd, xwd, zd, zwd, zed, Ld, zbound_up, zbound_down)
        return f

    def int_k0_gauss_DW_(self, u, X, Y):
        """
        :param u: f(s)
        :param X: координата точки
        :param Y: координата точки
        :return:
        """
        N = 8
        a = -1
        b = 1
        m = 5
        Points_xi = np.array([-0.96028986, -0.79666648, -0.52553242, -0.18343464, 0.18343464, 0.52553242, 0.79666648, 0.96028986])
        Coeff_ci = np.array([0.10122854, 0.22238104, 0.31370664, 0.36268378, 0.36268378, 0.31370664, 0.22238104, 0.10122854])
        summa = 0
        for j in range(0, m):
            for i in range(0, N):
                xi = (a + (b - a) / m * j + a + (b - a) / m * (j + 1)) / 2 + (b - a) / m / 2 * Points_xi[i]
                summa += Coeff_ci[i] * special.kn(0, u * sqrt((X - xi) ** 2 + Y ** 2))
        int_k0_gauss_DW_ = (b - a) * summa / 2 / m

        return int_k0_gauss_DW_

    def pd_cinco_hor_well_lapl_rect(self, S, n_seg, zd, xwd, ywd, zwd, xed, yed, zed, rwd, cfd, perf, N_perf, xbound,
                                    ybound, zbound_up, zbound_down):
        """
        :param S: переменная пространства Лапласа
        :param n_seg: число сегментов
        :param zd: безразмерная координата точки, в которой рассчитывается давление
        :param xwd: безразмерная координата скважины
        :param ywd: безразмерная координата скважины
        :param zwd: безразмерная координата скважины
        :param xed: безразмерная координата границы x
        :param yed: безразмерная координата границы y
        :param zed: безразмерная координата границы z
        :param rwd: безразмерный радиус скважины
        :param cfd: безразмерный поправочный коэффициент на эффект ствола относительно полудлины трещины
        :param perf: доля перфораций горизонтальной скважины
        :param N_perf: количество перфорированных интервалов
        :param xbound: тип границы
        :param ybound: тип границы
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z
        :return:
        """
        ql = self.cinco_hor_well_lapl_rect(S, n_seg, zd, xwd, ywd, zwd, xed, yed, zed, rwd, cfd, perf, N_perf, xbound,
                                           ybound, zbound_up, zbound_down)
        pd_cinco_hor_well_lapl_rect = ql[n_seg + 1]

        return pd_cinco_hor_well_lapl_rect

    def cinco_hor_well_lapl_rect(self, S, n_seg, zd, xwd, ywd, zwd, xed, yed, zed, rwd, cfd, perf, N_perf, xbound,
                                 ybound, zbound_up, zbound_down):
        """
        :param S: переменная пространства Лапласа
        :param n_seg: число сегментов
        :param zd: безразмерная координата точки, в которой рассчитывается давление
        :param xwd: безразмерная координата скважины
        :param ywd: безразмерная координата скважины
        :param zwd: безразмерная координата скважины
        :param xed: безразмерная координата границы x
        :param yed: безразмерная координата границы y
        :param zed: безразмерная координата границы z
        :param rwd: безразмерный радиус скважины
        :param cfd: безразмерный поправочный коэффициент на эффект ствола относительно полудлины трещины
        :param perf: доля перфораций горизонтальной скважины
        :param N_perf: количество перфорированных интервалов
        :param xbound: тип границы
        :param ybound: тип границы
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z
        :return:
        """
        flow_res_matr = np.zeros(shape=(n_seg + 1, n_seg + 1))
        infl_matr = np.zeros(shape=(n_seg + 1, n_seg + 1))
        main_matr = np.zeros(shape=(n_seg + 2, n_seg + 2))
        flow_res_matr = self.flow_resistance_hor(flow_res_matr, n_seg, cfd)
        rght = self.rhs(n_seg, 1)
        infl_matr = self.influence_matrix_hor_well_rect(S, infl_matr, n_seg, zd, xwd, ywd, zwd, xed, yed, zed, rwd,
                                                        perf, N_perf, xbound,
                                                        ybound, zbound_up, zbound_down)
        main_matr = self.main_matrix(S, main_matr, flow_res_matr, infl_matr, n_seg)
        rght = self.LU_(main_matr, rght, n_seg + 1)
        cinco_hor_well_lapl_rect = rght

        return cinco_hor_well_lapl_rect

    def flow_resistance_hor(self, mm, N, cfd):
        """
        :param mm: матрица сопротивлений
        :param N: число сегментов
        :param cfd: безразмерный поправочный коэффициент на эффект ствола относительно полудлины трещины
        :return: матрица сопротивлений
        """
        dx = 2 / N
        a = dx / cfd
        for i in range(1, N + 1):
            for k in range(1, i):
                mm[i, k] = a * (k - 1 / 2)
            mm[i, i] = a * (i - 5 / 8)
            for k in range(i + 1, N + 1):
                mm[i, k] = a * (i - 1 / 2)
        return mm

    def rhs(self, N, Qd=1):
        """
        :param N: размерность вектора
        :param Qd: общий дебит
        :return:
        """
        rght = np.zeros(shape=(N + 2))
        rght[-1] = Qd

        return rght

    def influence_matrix_hor_well_rect(self, S, nn, n_seg, zd, xwd, ywd, zwd, xed, yed, zed, rwd, perf, N_perf, xbound,
                                       ybound, zbound_up, zbound_down):
        """
        :param S: переменная пространства Лапласа
        :param nn: матрица влияния
        :param n_seg: число сегментов
        :param zd: безразмерная координата точки, в которой рассчитывается давление
        :param xwd: безразмерная координата скважины
        :param ywd: безразмерная координата скважины
        :param zwd: безразмерная координата скважины
        :param xed: безразмерная координата границы x
        :param yed: безразмерная координата границы y
        :param zed: безразмерная координата границы z
        :param rwd: безразмерный радиус скважины
        :param perf: доля перфораций горизонтальной скважины
        :param N_perf: количество перфорированных интервалов
        :param xbound: тип границы
        :param ybound: тип границы
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z
        :return: матрица влияний
        """
        alpha = perf
        dx = 2 * alpha / n_seg
        dx1 = 2 * (1 - alpha) / N_perf
        unit_length = dx / 2
        Sd = S * unit_length ** 2
        xedd = xed / unit_length
        yedd = yed / unit_length
        ywdd = ywd / unit_length
        ydd = ywdd
        zedd = zed / unit_length
        zwdd = zwd / unit_length
        zdd = zd / unit_length
        perf_coeff = n_seg / N_perf
        arr_coord = np.zeros(shape=n_seg)
        arr_coord[0] = (xwd - 1 + dx / 2) / unit_length
        if perf_coeff == 1:
            for i in range(1, n_seg):
                arr_coord[i] = arr_coord[i - 1] + (dx + dx1) / unit_length
        else:
            for i in range(1, n_seg):
                if i % perf_coeff == 0:
                    delta = dx1 + dx
                else:
                    delta = dx
                arr_coord[i] = arr_coord[i - 1] + delta / unit_length
        for i in range(1, n_seg + 1):
            for j in range(1, n_seg + 1):
                nn[i, j] = self.horizontal_well_for_cinco(Sd, arr_coord[i - 1], ydd, zdd, arr_coord[j - 1], ywdd, zwdd,
                                                          xedd, yedd, zedd, xbound, ybound, zbound_up, zbound_down)
        return nn

    def main_matrix(self, S, cc, mm, nn, N):
        """
        :param S: переменная пространства Лапласа
        :param cc: главная матрица
        :param mm: матрица сопротивлений
        :param nn: матрица влияния
        :param N: размерность матриц
        :return: главная матрица
        """
        cc[1:N + 1, N + 1] = -1
        cc[:N + 1, :N + 1] = cc[:N + 1, :N + 1] + mm[:N + 1, :N + 1] + nn[:N + 1, :N + 1]
        cc[N + 1, 1:N + 1] = 1

        return cc

    def LU_(self, a, b, N):
        """
        :param a: матрица
        :param b: вектор
        :param N: размерность вектора
        :return:
        """
        indx = np.zeros(N + 1)
        a = self.LU_decomposition(a, indx, N)
        b = self.LU_solve(a, b, indx, N)

        return b

    def LU_decomposition(self, a, indx, N):
        """
        :param a: матрица
        :param indx: вспомогательный вектор
        :param N: размерность вектора
        :return:
        """
        tiny = 1E-20
        vv = np.zeros(shape=N + 1)
        d = 1
        for i in range(1, N + 1):
            big = 0
            for j in range(1, N + 1):
                temp = abs(a[i, j])
                if temp > big:
                    big = temp
            if big == 0:
                return a
            vv[i] = 1 / big
        for j in range(1, N + 1):
            for i in range(1, j):
                summ = a[i, j]
                for k in range(1, i):
                    summ = summ - a[i, k] * a[k, j]
                a[i, j] = summ
            big = 0
            for i in range(j, N + 1):
                summ = a[i, j]
                for k in range(1, j):
                    summ = summ - a[i, k] * a[k, j]
                a[i, j] = summ
                dum = vv[i] * abs(summ)
                if dum >= big:
                    big = dum
                    imax = i
            if j != imax:
                for k in range(1, N + 1):
                    dum = a[imax, k]
                    a[imax, k] = a[j, k]
                    a[j, k] = dum
                d = -d
                vv[imax] = vv[j]
            indx[j] = imax
            if a[j, j] == 0:
                a[j, j] = tiny
            if j != N:
                dum = 1 / a[j, j]
                for i in range(j + 1, N + 1):
                    a[i, j] = a[i, j] * dum

        return a

    def LU_solve(self, a, b, indx, N):
        """
        :param a: матрица
        :param b: вектор
        :param indx: вспомогательный вектор
        :param N: размерность вектора
        :return:
        """
        ii = 0
        for i in range(1, N + 1):
            ip = int(indx[i])
            summ = b[ip]
            b[ip] = b[i]
            if ii:
                for j in range(ii, i):
                    summ = summ - a[i, j] * b[j]
            elif summ:
                ii = i
            b[i] = summ
        for i in range(N, 0, -1):
            summ = b[i]
            for j in range(i + 1, N + 1):
                summ = summ - a[i, j] * b[j]
            b[i] = summ / a[i, i]

        return b

    def horizontal_well_for_cinco(self, S, xd, yd, zd, xwd, ywd, zwd, xed, yed, zed, xbound, ybound, zbound_up,
                                  zbound_down):
        """
        :param S: переменная пространства Лапласа
        :param xd: безразмерная координата точки, в которой рассчитывается давление
        :param yd: безразмерная координата точки, в которой рассчитывается давление
        :param zd: безразмерная координата точки, в которой рассчитывается давление
        :param xwd: безразмерная координата скважины
        :param ywd: безразмерная координата скважины
        :param zwd: безразмерная координата скважины
        :param xed: безразмерная координата границы x
        :param yed: безразмерная координата границы y
        :param zed: безразмерная координата границы z
        :param xbound: тип границы
        :param ybound: тип границы
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z
        :return:
        """
        Ld = 1 / zed
        compl_type_hor = "frac"
        if abs(yd - ywd) < self.tiny_2:
            if zbound_up == "n" and zbound_down == "n":
                if S > 0:
                    pwd = self.pd_lapl_rect(S, xd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type_hor) + \
                          self.f(S, xd, xwd, zd, zwd, zed, Ld, yd, ywd, zbound_up, zbound_down) + \
                          self.F_b1_2(S, xed, yd, ywd, yed, zd, zwd, zed, Ld, zbound_up, zbound_down, xbound, ybound, True) + \
                          self.F_b2_2(S, yd, ywd, yed, xd, xwd, xed, zd, zwd, zed, Ld, zbound_up, zbound_down, xbound, ybound, True) - \
                          self.F_b3_2(S, xd, xwd, xed, zd, zwd, zed, yd, ywd, Ld, zbound_up, zbound_down, xbound)
                else:
                    pwd = self.pd_lapl_rect(S, xd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type_hor) + self.f(S,
                                                                                                                    xd,
                                                                                                                    xwd,
                                                                                                                    zd,
                                                                                                                    zwd,
                                                                                                                    zed,
                                                                                                                    Ld,
                                                                                                                    yd,
                                                                                                                    ywd,
                                                                                                                    zbound_up,
                                                                                                                    zbound_down) + self.F_b1_2(
                        S, xed, yd, ywd, yed, zd, zwd, zed, Ld, zbound_up, zbound_down, xbound, ybound, True) + \
                          self.F_b2_2(S, yd, ywd, yed, xd, xwd, xed, zd, zwd, zed, Ld, zbound_up, zbound_down, xbound,
                                      ybound,
                                      True) - self.F_b3_2(S, xd, xwd, xed, zd, zwd, zed, yd, ywd, Ld, zbound_up,
                                                          zbound_down, xbound)
            else:
                if S > 0:
                    pwd = self.f(S, xd, xwd, zd, zwd, zed, Ld, yd, ywd, zbound_up, zbound_down) + self.F_b1_2(S, xed,
                                                                                                              yd, ywd,
                                                                                                              yed, zd,
                                                                                                              zwd, zed,
                                                                                                              Ld,
                                                                                                              zbound_up,
                                                                                                              zbound_down,
                                                                                                              xbound,
                                                                                                              ybound,
                                                                                                              True) + \
                          self.F_b2_2(S, yd, ywd, yed, xd, xwd, xed, zd, zwd, zed, Ld, zbound_up, zbound_down, xbound,
                                      ybound,
                                      True) - self.F_b3_2(S, xd, xwd, xed, zd, zwd, zed, yd, ywd, Ld, zbound_up,
                                                          zbound_down, xbound)
                else:
                    pwd = self.f(-1, xd, xwd, zd, zwd, zed, Ld, yd, ywd, zbound_up, zbound_down) + self.F_b1_2(0, xed,
                                                                                                               yd, ywd,
                                                                                                               yed, zd,
                                                                                                               zwd, zed,
                                                                                                               Ld,
                                                                                                               zbound_up,
                                                                                                               zbound_down,
                                                                                                               xbound,
                                                                                                               ybound,
                                                                                                               True) + \
                          self.F_b2_2(0, yd, ywd, yed, xd, xwd, xed, zd, zwd, zed, Ld, zbound_up, zbound_down, xbound,
                                      ybound,
                                      True) - self.F_b3_2(0, xd, xwd, xed, zd, zwd, zed, yd, ywd, Ld, zbound_up,
                                                          zbound_down, xbound)
        else:
            if zbound_up == "n" and zbound_down == "n":
                if S > 0:
                    pwd = self.pd_lapl_rect(S, xd, xwd, xed, yd, ywd, yed, xbound, ybound,
                                            compl_type_hor) + self.F_b1_2(S, xed, yd,
                                                                          ywd, yed, zd,
                                                                          zwd, zed, Ld,
                                                                          zbound_up,
                                                                          zbound_down,
                                                                          xbound, ybound,
                                                                          False) + \
                          self.F_b2_2(S, yd, ywd, yed, xd, xwd, xed, zd, zwd, zed, Ld, zbound_up, zbound_down, xbound,
                                      ybound,
                                      False)
                else:
                    pwd = self.pd_rect_BD(xd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type_hor) + self.F_b1_2(0,
                                                                                                                    xed,
                                                                                                                    yd,
                                                                                                                    ywd,
                                                                                                                    yed,
                                                                                                                    zd,
                                                                                                                    zwd,
                                                                                                                    zed,
                                                                                                                    Ld,
                                                                                                                    zbound_up,
                                                                                                                    zbound_down,
                                                                                                                    xbound,
                                                                                                                    ybound,
                                                                                                                    False) + \
                          self.F_b2_2(0, yd, ywd, yed, xd, xwd, xed, zd, zwd, zed, Ld, zbound_up, zbound_down, xbound,
                                      ybound,
                                      False)
            else:
                if S > 0:
                    pwd = self.F_b1_2(S, xed, yd, ywd, yed, zd, zwd, zed, Ld, zbound_up, zbound_down, xbound, ybound,
                                      False) + self.F_b2_2(S, yd, ywd, yed, xd, xwd, xed, zd, zwd, zed, Ld, zbound_up,
                                                           zbound_down, xbound, ybound,
                                                           False)
                else:
                    pwd = self.F_b1_2(0, xed, yd, ywd, yed, zd, zwd, zed, Ld, zbound_up, zbound_down, xbound, ybound,
                                      False) + self.F_b2_2(0, yd, ywd, yed, xd, xwd, xed, zd, zwd, zed, Ld, zbound_up,
                                                           zbound_down, xbound, ybound,
                                                           False)

        horizontal_well_for_cinco = pwd

        return horizontal_well_for_cinco

    def unit_cylinder_source(self, S, rd):
        """
        :param S: переменная пространства Лапласа
        :param rd: безразмерное расстояние от центра скважины
        :return: давление
        """
        if S > 0:
            u = S ** (1 / 2)
            lapl_rd = rd * u
            if rd >= 5 and lapl_rd >= 100:
                return 0
            if lapl_rd / 2 > self.small_bessel_arg:
                result = (
                        1 / (2 * pi) * special.kn(0, lapl_rd) / (u * special.kn(1, u))
                )
                return result
        else:
            lapl_rd = rd
        result = 1 / (2 * pi) * (-log(lapl_rd / 2) - self.gamma)
        return result

    def sumexp(self, arg, yed):
        """
        Суммирование экспоненциального ряда
        :param arg: аргумент
        :param yed: безразмерная координата границы y
        :return: сумма экспоненциального ряда
        """
        if arg * yed <= 0:
            raise ValueError("sumexp: arg * yed <= 0")
        x = np.clip(2 * arg * yed, None, 700)
        result = 1 / (np.exp(x) - 1)
        return result

    def pd_b2_k(self, k, S, xd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type, subtract_inf):
        """
        :param k: номер слагаемого в функции pd_b3
        :param S: переменная пространства Лапласа
        :param xd: безразмерная координата точки, в которой рассчитывается давление
        :param xwd: безразмерная координата скважины
        :param xed: безразмерная координата границы x
        :param yd: безразмерная координата точки, в которой рассчитывается давление
        :param ywd: безразмерная координата скважины
        :param yed: безразмерная координата границы y
        :param xbound: тип границы
        :param ybound: тип границы
        :param compl_type: тип скважины
        :param subtract_inf: параметр, определяющий, нужно ли считать часть уравнения
        :return: давление
        """
        if xbound == "n":
            part_1 = 2 / xed * cos(k * pi * xd / xed) * cos(k * pi * xwd / xed)
        elif xbound == "c":
            part_1 = 2 / xed * sin(k * pi * xd / xed) * sin(k * pi * xwd / xed)
        else:
            part_1 = 0
        if compl_type == "frac":
            part_1 = 2 / 2 * xed / pi / k * sin(k * pi / xed) * part_1
        elif compl_type == "vert":
            part_1 = part_1
        else:
            pd_b2_k = 0
            return pd_b2_k
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
        if self.extend_reflections_mode and subtract_inf:
            sbtr_inf_ext = 0
        else:
            sbtr_inf_ext = 1
        ek = sqrt(S + k ** 2 * pi ** 2 / xed ** 2)
        smexp = self.sumexp(ek, yed)
        part_2 = exp(-ek * abs(yd - ywd)) * sbtr_inf + (
                signum * exp(-ek * (yd + ywd)) + signum * exp(-ek * (2 * yed - (yd + ywd))) +
                exp(-ek * (2 * yed - (ywd - yd))) + exp(-ek * (2 * yed - (yd - ywd)))) * (sbtr_inf_ext + smexp)
        part_2 = 1 / (2 * ek) * part_2
        pd_b2_k = part_1 * part_2

        return pd_b2_k

    def pd_b2_2(self, S, xd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type, subtract_inf=True):
        """
        :param S: переменная пространства Лапласа
        :param xd: безразмерная координата точки, в которой рассчитывается давление
        :param xwd: безразмерная координата скважины
        :param xed: безразмерная координата границы x
        :param yd: безразмерная координата точки, в которой рассчитывается давление
        :param ywd: безразмерная координата скважины
        :param yed: безразмерная координата границы y
        :param xbound: тип границы
        :param ybound: тип границы
        :param compl_type: тип скважины
        :param subtract_inf: параметр, определяющий, нужно ли считать часть уравнения
        :return: давление
        """
        summa = 0
        k = 1
        while True:
            psum = 0
            psumabs = 0
            for i in range(1, self.part_sum_num + 1):
                add = self.pd_b2_k(k, S, xd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type, subtract_inf)
                psum = psum + add
                psumabs = psumabs + abs(add)
                k = k + 1
            summa = summa + psum
            if not (abs(psumabs) / self.part_sum_num >= self.tiny_2 * (abs(summa) + self.tiny_2) and k < self.max_it):
                break

        return summa

    def F_b1_2(self, S, xed, yd, ywd, yed, zd, zwd, zed, Ld, zbound_up,
               zbound_down, xbound: str, ybound, subtract_inf):
        """
        :param S: переменная пространства Лапласа
        :param xed: безразмерная координата границы x
        :param yd: безразмерная координата точки, в которой рассчитывается давление
        :param ywd: безразмерная координата скважины
        :param yed: безразмерная координата границы y
        :param zd: безразмерная координата точки, в которой рассчитывается давление
        :param zwd: безразмерная координата скважины
        :param zed: безразмерная координата границы z
        :param Ld: безразмерная длина
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z
        :param xbound: тип границы
        :param ybound: тип границы
        :param subtract_inf: параметр, определяющий, нужно ли считать часть уравнения
        :return:
        """
        if xbound == "c":
            F_b1_2 = 0
        else:
            summa = series_fn.nsum(
                lambda k: self.F_b1_k(k, S, yd, ywd, yed, zd, zwd, zed, Ld, zbound_up, zbound_down, ybound,
                                      subtract_inf),
                1,
                self.max_it - 1,
            )
            F_b1_2 = 1 / xed * summa

        return F_b1_2

    def F_b2_2(self, S, yd, ywd, yed, xd, xwd, xed, zd, zwd, zed, Ld, zbound_up, zbound_down, xbound, ybound, subtract_inf):
        """
        :param S: переменная пространства Лапласа
        :param yd: безразмерная координата точки, в которой рассчитывается давление
        :param ywd: безразмерная координата скважины
        :param yed: безразмерная координата границы y
        :param xd: безразмерная координата точки, в которой рассчитывается давление
        :param xwd: безразмерная координата скважины
        :param xed: безразмерная координата границы x
        :param zd: безразмерная координата точки, в которой рассчитывается давление
        :param zwd: безразмерная координата скважины
        :param zed: безразмерная координата границы z
        :param Ld: безразмерная длина
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z
        :param xbound: тип границы
        :param ybound: тип границы
        :param subtract_inf: параметр, определяющий, нужно ли считать часть уравнения
        :return:
        """
        sum_ = series_fn.nsum(
            lambda k: self.F_b2_k_3(
                k, S, yd, ywd, yed, xd, xwd, xed, zd, zwd, zed, Ld, zbound_up, zbound_down,
                ybound, subtract_inf),
            1,
            self.max_it - 1,
        )
        return sum_

    def F_b2_k_3(self, p, S, yd, ywd, yed, xd, xwd, xed, zd, zwd, zed, Ld, zbound_up, zbound_down, ybound,
                 subtract_inf):
        """
        :param p: индекс суммирования
        :param S: переменная пространства Лапласа
        :param yd: безразмерная координата точки, в которой рассчитывается давление
        :param ywd: безразмерная координата скважины
        :param yed: безразмерная координата границы y
        :param xd: безразмерная координата точки, в которой рассчитывается давление
        :param xwd: безразмерная координата скважины
        :param xed: безразмерная координата границы x
        :param zd: безразмерная координата точки, в которой рассчитывается давление
        :param zwd: безразмерная координата скважины
        :param zed: безразмерная координата границы z
        :param Ld: безразмерная длина
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z
        :param ybound: тип границы
        :param subtract_inf: параметр, определяющий, нужно ли считать часть уравнения
        :return:
        """
        if (zbound_up == "c" and zbound_down == "n") or (zbound_up == "n" and zbound_down == "c"):
            N = (2 * p - 1) / 2
        else:
            N = p
        if zbound_up == "c" and zbound_down == "c":
            part_1 = sin(N * pi * zd / zed) * sin(N * pi * zwd / zed)
        else:
            part_1 = cos(N * pi * zd / zed) * cos(N * pi * zwd / zed)
        F_b2_k_3 = part_1 * self.F_b2_k_2(p, S, yd, ywd, yed, xd, xwd, xed, Ld, zbound_up, zbound_down,
                                          ybound, subtract_inf)
        return F_b2_k_3

    def F_b2_k_2(self, p, S, yd, ywd, yed, xd, xwd, xed, Ld, zbound_up, zbound_down, ybound, subtract_inf):
        """
        :param p: индекс суммирования
        :param S: переменная пространства Лапласа
        :param yd: безразмерная координата точки, в которой рассчитывается давление
        :param ywd: безразмерная координата скважины
        :param yed: безразмерная координата границы y
        :param xd: безразмерная координата точки, в которой рассчитывается давление
        :param xwd: безразмерная координата скважины
        :param xed: безразмерная координата границы x
        :param Ld: безразмерная длина
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z
        :param ybound: тип границы
        :param subtract_inf: параметр, определяющий, нужно ли считать часть уравнения
        :return:
        """
        sum_ = series_fn.nsum(
            lambda k: self.F_b2_k_1(p, k, S, yd, ywd, yed, xd, xwd, xed, Ld, zbound_up, zbound_down, ybound,
                                    subtract_inf),
            1,
            self.max_it - 1,
        )
        return sum_

    def F_b2_k_1(self, p, m, S, yd, ywd, yed, xd, xwd, xed, Ld, zbound_up, zbound_down, ybound, subtract_inf):
        """
        :param p: индекс внешнего суммирования
        :param m: индекс внутреннего суммирования
        :param S: переменная пространства Лапласа
        :param yd: безразмерная координата точки, в которой рассчитывается давление
        :param ywd: безразмерная координата скважины
        :param yed: безразмерная координата границы y
        :param xd: безразмерная координата точки, в которой рассчитывается давление
        :param xwd: безразмерная координата скважины
        :param xed: безразмерная координата границы x
        :param Ld: безразмерная длина
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z
        :param ybound: тип границы
        :param subtract_inf: параметр, определяющий, нужно ли считать часть уравнения
        :return:
        """
        sbtr_inf = 0 if subtract_inf else 1
        k = m

        N = (2 * p - 1) / 2 if (zbound_up, zbound_down) in [("c", "n"), ("n", "c")] else p

        k_pi_over_xed = k * pi / xed
        if ybound == "n":
            part_1 = cos(k_pi_over_xed * xd) * cos(k_pi_over_xed * xwd)
            sign = 1
        elif ybound == "c":
            part_1 = sin(k_pi_over_xed * xd) * sin(k_pi_over_xed * xwd)
            sign = -1
        else:
            part_1 = 0
            sign = 1

        e_n_k = sqrt(S + (N * pi * Ld) ** 2 + k_pi_over_xed ** 2)

        exp_e_n_k_yd_ywd = exp(-e_n_k * (yd + ywd))
        exp_e_n_k_2yed_yd_ywd = exp(-e_n_k * (2 * yed - yd - ywd))
        exp_e_n_k_2yed_abs_yd_ywd = exp(-e_n_k * (2 * yed - abs(yd - ywd)))
        exp_e_n_k_abs_yd_ywd = exp(-e_n_k * abs(yd - ywd))

        sumexp_e_n_k_yed = self.sumexp(e_n_k, yed)

        exp_sum = sign * (exp_e_n_k_yd_ywd + exp_e_n_k_2yed_yd_ywd) + exp_e_n_k_2yed_abs_yd_ywd
        factor = 1 + sumexp_e_n_k_yed

        return sin(k_pi_over_xed) * part_1 / e_n_k / k * (
                exp_sum * factor + exp_e_n_k_abs_yd_ywd * (sbtr_inf + sumexp_e_n_k_yed))

    def F_b1_k(self, k, S, yd, ywd, yed, zd, zwd, zed, Ld, zbound_up, zbound_down, ybound, subtract_inf):
        """
        :param k: индекс суммирования
        :param S: переменная пространства Лапласа
        :param yd: безразмерная координата точки, в которой рассчитывается давление
        :param ywd: безразмерная координата скважины
        :param yed: безразмерная координата границы y
        :param zd: безразмерная координата точки, в которой рассчитывается давление
        :param zwd: безразмерная координата скважины
        :param zed: безразмерная координата границы z
        :param Ld: безразмерная длина
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z
        :param ybound: тип границы
        :param subtract_inf: параметр, определяющий, нужно ли считать часть уравнения
        :return:
        """
        sbtr_inf = 0 if subtract_inf else 1

        N = (2 * k - 1) / 2 if (zbound_up, zbound_down) in [("c", "n"), ("n", "c")] else k

        zd_zed_ratio = N * pi * zd / zed
        zwd_zed_ratio = N * pi * zwd / zed

        if zbound_up == "c" and zbound_down == "c":
            part_1 = sin(zd_zed_ratio) * sin(zwd_zed_ratio)
        else:
            part_1 = cos(zd_zed_ratio) * cos(zwd_zed_ratio)

        sign = -1 if ybound == "c" else 1

        N_pi_Ld = N * pi * Ld
        e_n = sqrt(S + N_pi_Ld ** 2)

        exp_neg_e_n_yd_ywd = exp(-e_n * (yd + ywd))
        exp_neg_e_n_2yed_yd_ywd = exp(-e_n * (2 * yed - yd - ywd))
        exp_neg_e_n_2yed_abs_yd_ywd = exp(-e_n * (2 * yed - abs(yd - ywd)))
        exp_neg_e_n_abs_yd_ywd = exp(-e_n * abs(yd - ywd))

        sumexp_e_n_yed = self.sumexp(e_n, yed)

        exp_sum = sign * (exp_neg_e_n_yd_ywd + exp_neg_e_n_2yed_yd_ywd) + exp_neg_e_n_2yed_abs_yd_ywd
        factor = 1 + sumexp_e_n_yed

        F_b1_k = (part_1 / e_n) * (exp_sum * factor + exp_neg_e_n_abs_yd_ywd * (sbtr_inf + sumexp_e_n_yed))

        return F_b1_k

    def F_b3_2(self, S, xd, xwd, xed, zd, zwd, zed, yd, ywd, Ld, zbound_up, zbound_down, xbound):
        """
        :param S: переменная пространства Лапласа
        :param xd: безразмерная координата точки, в которой рассчитывается давление
        :param xwd: безразмерная координата скважины
        :param xed: безразмерная координата границы x
        :param zd: безразмерная координата точки, в которой рассчитывается давление
        :param zwd: безразмерная координата скважины
        :param zed: безразмерная координата границы z
        :param yd: безразмерная координата точки, в которой рассчитывается давление
        :param ywd: безразмерная координата скважины
        :param Ld: безразмерная длина
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z
        :param xbound: тип границы
        :return:
        """
        sum_ = series_fn.nsum(
            lambda k: self.F_b3_k(k, S, xd, xwd, xed, zd, zwd, zed, Ld, zbound_up, zbound_down, xbound),
            1,
            self.max_it - 1,
        )
        result = 1 / 2 / pi * sum_
        return result

    def F_b3_k(self, k, S, xd, xwd, xed, zd, zwd, zed, Ld, zbound_up, zbound_down, xbound):
        """
        :param k: индекс суммирования
        :param S: переменная пространства Лапласа
        :param xd: безразмерная координата точки, в которой рассчитывается давление
        :param xwd: безразмерная координата скважины
        :param xed: безразмерная координата границы x
        :param zd: безразмерная координата точки, в которой рассчитывается давление
        :param zwd: безразмерная координата скважины
        :param zed: безразмерная координата границы z
        :param Ld: безразмерная длина
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z
        :param xbound: тип границы
        :return:
        """
        N = (2 * k - 1) / 2 if (zbound_up, zbound_down) in [("c", "n"), ("n", "c")] else k

        zd_ratio = N * pi * zd / zed
        zwd_ratio = N * pi * zwd / zed
        part_1 = (sin(zd_ratio) * sin(zwd_ratio) if zbound_up == "c" and zbound_down == "c"
                  else cos(zd_ratio) * cos(zwd_ratio))

        sign = -1 if xbound == "c" else 1
        e_n = sqrt(S + (N * pi * Ld) ** 2)

        summa = series_fn.nsum(
            lambda p: (self.F_b3_k_calc(xd, xwd, xed, p, e_n, sign)), 1, np.inf,
        )

        final_term1 = e_n * (xd + xwd + 1)
        final_term2 = e_n * (xd + xwd - 1)
        F_b3_k = part_1 / e_n * (sign * self.ki1(final_term1) - sign * self.ki1(final_term2) + summa)

        return F_b3_k

    def F_b3_k_calc(self, xd, xwd, xed, p, e_n, sign):
        tmp = 2 * p * xed + xd
        exp_term1 = e_n * (tmp - xwd + 1)
        exp_term2 = e_n * (tmp - xwd - 1)
        exp_term3 = e_n * (tmp + xwd + 1)
        exp_term4 = e_n * (tmp + xwd - 1)
        res = (-self.integrals(p, xd, xwd, xed, e_n, 1) + self.ki1(exp_term1) - self.ki1(exp_term2)
               - sign * self.integrals(p, xd, xwd, xed, e_n, 2) + sign * self.ki1(exp_term3) - sign * self.ki1(
                    exp_term4))
        return res

    def integrals(self, p, xd, xwd, xed, e_n, a):
        """
        :param p: индекс суммирования
        :param xd: безразмерная координата точки, в которой рассчитывается давление
        :param xwd: безразмерная координата скважины
        :param xed: безразмерная координата границы x
        :param e_n: эпсилон_n
        :param a: переменная, отвечающая за выбор случая
        :return: интеграл
        """
        if a == 1:
            if abs(xd - xwd - 2 * p * xed) <= 1:
                integrals = pi - self.ki1(e_n * (xd - xwd - 2 * p * xed + 1)) - self.ki1(
                    e_n * (1 - xd + xwd + 2 * p * xed))
            else:
                integrals = self.ki1(e_n * (abs(xd - xwd - 2 * p * xed) - 1)) - self.ki1(
                    e_n * (abs(xd - xwd - 2 * p * xed) + 1))
        elif a == 2:
            if abs(xd - xwd - 2 * p * xed) <= 1:
                integrals = pi - self.ki1(e_n * (xd + xwd - 2 * p * xed + 1)) - self.ki1(
                    e_n * (1 - xd - xwd + 2 * p * xed))
            else:
                integrals = self.ki1(e_n * (abs(xd + xwd - 2 * p * xed) - 1)) - self.ki1(
                    e_n * (abs(xd + xwd - 2 * p * xed) + 1))

        return integrals

    def ki1(self, X):
        """
        :param X: аргумент
        :return:
        """
        ki1 = pi / 2 - self.fast_ik0(X)
        return ki1

    def fast_ik0(self, X):
        """
        :param X: аргумент
        :return:
        """
        tiny = 1E-17
        if X <= 0:
            return 0
        if X > 18:
            return pi * 0.5

        x2 = X * 0.5
        x2_k = 1
        kfct = 1
        part1 = 1
        part2 = 1
        part3 = 0
        k = 0
        sum1n = 0

        while True:
            k += 1
            kfct *= k
            x2_k *= x2
            inc = (x2_k / kfct) ** 2
            div2k = 1 / (k + k + 1)
            inc *= div2k
            part1 = part1 + inc
            part2 = part2 + inc * div2k
            sum1n = sum1n + 1 / k
            part3 = part3 + inc * sum1n
            if not ((30 * inc / part3) >= tiny):
                break

        ans = -(log(x2) + self.gamma) * part1 * X + part2 * X + part3 * X
        return ans

    def f_4_2(self, xd, xwd, zd, zwd, zed, Ld, zbound_up, zbound_down):
        """
        :param xd: безразмерная координата точки, в которой рассчитывается давление
        :param xwd: безразмерная координата скважины
        :param zd: безразмерная координата точки, в которой рассчитывается давление
        :param zwd: безразмерная координата скважины
        :param zed: безразмерная координата границы z
        :param Ld: безразмерная длина
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z
        :return:
        """
        summa = series_fn.nsum(
            lambda k: self.f_4_k(k, xd, xwd, zd, zwd, zed, Ld, zbound_up, zbound_down),
            1,
            self.max_it - 1,
        )
        f_4_2 = summa / 2 / pi

        return f_4_2

    def f_4_k(self, k, xd, xwd, zd, zwd, zed, Ld, zbound_up, zbound_down):
        """
        :param k: индекс суммирования
        :param xd: безразмерная координата точки, в которой рассчитывается давление
        :param xwd: безразмерная координата скважины
        :param zd: безразмерная координата точки, в которой рассчитывается давление
        :param zwd: безразмерная координата скважины
        :param zed: безразмерная координата границы z
        :param Ld: безразмерная длина
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z
        :return:
        """
        if (zbound_up == "c" and zbound_down == "n") or (zbound_up == "n" and zbound_down == "c"):
            N = (2 * k - 1) / 2
        else:
            N = k
        if zbound_up == "c" and zbound_down == "c":
            part_1 = sin(N * pi * zd / zed) * sin(N * pi * zwd / zed)
        else:
            part_1 = cos(N * pi * zd / zed) * cos(N * pi * zwd / zed)
        e_n = N * pi * Ld
        if abs(xd - xwd) <= 1:
            add = part_1 / e_n * (pi - (self.ki1(e_n * (1 + xd - xwd)) + self.ki1(e_n * (1 - xd + xwd))))
        else:
            add = part_1 / e_n * (self.ki1(e_n * (abs(xd - xwd) - 1)) - self.ki1(e_n * (abs(xd - xwd) + 1)))
        f_4_k = add

        return f_4_k

    def f_3_2(self, S, xd, xwd, zd, zwd, zed, Ld, zbound_up, zbound_down):
        """
        :param S: переменная пространства Лапласа
        :param xd: безразмерная координата точки, в которой рассчитывается давление
        :param xwd: безразмерная координата скважины
        :param zd: безразмерная координата точки, в которой рассчитывается давление
        :param zwd: безразмерная координата скважины
        :param zed: безразмерная координата границы z
        :param Ld: безразмерная длина
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z
        :return:
        """
        summa = series_fn.nsum(
            lambda k: self.f_3_k(k, S, xd, xwd, zd, zwd, zed, Ld, zbound_up, zbound_down),
            1,
            self.part_sum_num_f_3 + 1,
        )
        if (zbound_up == "c" and zbound_down == "n") or (zbound_up == "n" and zbound_down == "c") or (
                zbound_up == "c" and zbound_down == "c"):
            f_3_2 = 2 * summa
        else:
            if abs(xd - xwd) <= 1:
                f_3_2 = pi / sqrt(S) + 2 * summa
            else:
                f_3_2 = 2 * summa

        return f_3_2

    def f_3_k(self, k, S, xd, xwd, zd, zwd, zed, Ld, zbound_up, zbound_down):
        """
        :param k: индекс суммирования
        :param S: переменная пространства Лапласа
        :param xd: безразмерная координата точки, в которой рассчитывается давление
        :param xwd: безразмерная координата скважины
        :param zd: безразмерная координата точки, в которой рассчитывается давление
        :param zwd: безразмерная координата скважины
        :param zed: безразмерная координата границы z
        :param Ld: безразмерная длина
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z
        :return:
        """
        if (zbound_up == "c" and zbound_down == "n") or (zbound_up == "n" and zbound_down == "c"):
            N = (2 * k - 1) / 2
        else:
            N = k
        if zbound_up == "c" and zbound_down == "c":
            part_1 = sin(N * pi * zd / zed) * sin(N * pi * zwd / zed)
        else:
            part_1 = cos(N * pi * zd / zed) * cos(N * pi * zwd / zed)
        e_n = sqrt(S + (N * pi * Ld) ** 2)
        if abs(xd - xwd) <= 1:
            add = part_1 / e_n * (self.ki1(e_n * (1 + xd - xwd)) +
                                  self.ki1(e_n * (1 - xd + xwd)))
        else:
            add = part_1 / e_n * (self.ki1(e_n * (abs(xd - xwd) - 1)) - self.ki1(e_n * (1 + abs(xd - xwd))))
        f_3_k = add

        return f_3_k

    def f_2(self, S, zd, zwd, zed, Ld, zbound_up, zbound_down):
        """
        :param S: переменная пространства Лапласа
        :param zd: безразмерная координата точки, в которой рассчитывается давление
        :param zwd: безразмерная координата скважины
        :param zed: безразмерная координата границы z
        :param Ld: безразмерная длина
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z
        :return:
        """
        add1 = 1
        add2 = 1
        summa = 0
        u = sqrt(S)
        k = 1
        if zbound_up == "n" and zbound_down == "n":
            sign1 = 1
        elif zbound_up == "c" and zbound_down == "c":
            sign1 = -1

        zd_zed = zd / zed
        zwd_zed = zwd / zed
        common_term1 = zd_zed + zwd_zed
        common_term2 = zd_zed - zwd_zed

        if (zbound_up == "c" and zbound_down == "n") or (zbound_up == "n" and zbound_down == "c"):
            while add1 + add2 >= 1E-20:
                for sign in [-1, 1]:
                    arg1 = abs(common_term1 + sign * 2 * k) * u / Ld
                    arg2 = abs(common_term2 / zed + sign * 2 * k) * u / Ld
                    add1 = special.kn(0, arg1)
                    add2 = special.kn(0, arg2)
                    summa += (add1 + add2) * (-1) ** k
                k = k + 1
            summa += special.kn(0, abs(common_term1) * u / Ld) + special.kn(0, abs(common_term2) * u / Ld)
        else:
            while add1 + add2 >= 1E-20:
                for sign in [-1, 1]:
                    arg1 = abs(common_term1 + sign * 2 * k) * u / Ld
                    arg2 = abs(common_term2 + sign * 2 * k) * u / Ld
                    add1 = special.kn(0, arg1)
                    add2 = special.kn(0, arg2)
                    summa += sign1 * add1 + add2
                k += 1
            summa += sign1 * special.kn(0, abs(common_term1) * u / Ld) + special.kn(0, abs(common_term2) * u / Ld)
        f_2 = 1 / Ld * summa

        return f_2
