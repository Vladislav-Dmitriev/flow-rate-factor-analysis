import numpy as np
from math import pi, log
from src.entities.Frac_calculator import Frac_Calculator


class Multifrac_Calculator(Frac_Calculator):
    """
    Класс расчёта параметров модели ГС с МГРП
    """
    def __init__(self, wellbore_wf, res_wf, grade, kf, hf, fracture_grow_t, rel_M, rel_D, R_in, S_choke, sf,
                 k, hw_f, extend_reflections_mode, tiny, tiny_2, large_s, small, part_sum_num, max_it, small_bessel_arg,
                 part_sum_num_b2_1, part_sum_num_f_4, n_frac, f_direction):
        super().__init__(wellbore_wf, res_wf, grade, kf, hf, fracture_grow_t, rel_M, rel_D, R_in, S_choke, sf, k, hw_f,
                         extend_reflections_mode, tiny, tiny_2, large_s, small, part_sum_num, max_it, small_bessel_arg,
                         part_sum_num_b2_1, part_sum_num_f_4)
        self.num_div_hor = 1
        self.n_frac = n_frac
        self.frac_dir = f_direction

    def calc_multifrac_geometry(self, xe, ye, wf, lf, l_hor, xbound, ybound):
        """
        :param xe: ширина прямоугольника
        :param ye: длина прямоугольника
        :param wf: относительное расстояние от длинной стороны прямоугольника до скважины, доля ширины прямоугольника
        :param lf: относительное расстояние от короткой стороны прямоугольника до скважины, доля длины прямоугольника
        :param l_hor: длина горизонтального ствола
        :param xbound: тип границы по x
        :param ybound: тип границы по y
        :return:
        """
        xw = np.zeros(self.n_frac + 1)
        xw[0] = xe * wf
        yw = np.zeros(self.n_frac + 1)
        yw[0] = ye * lf
        if self.n_frac == 1:
            xw[1] = xw[0]
            yw[1] = yw[0]
        else:
            xmin = xw[0] - l_hor / 2
            delta = l_hor / (self.n_frac - 1)
            for i in range(1, self.n_frac + 1):
                xw[i] = xmin + delta * (i - 1)
                yw[i] = yw[0]
        if self.frac_dir == "transverse":
            dummy = xe
            xe = ye
            ye = dummy

            dummyA = xw
            xw = yw
            yw = dummyA

            Dummy_bound = xbound
            xbound = ybound
            ybound = Dummy_bound

        return xw, yw, xe, ye, xbound, ybound

    def calc_well_size(self, l_hor):
        """
        :param l_hor: длина горизонтального ствола
        :return: размеры скважины
        """
        array_well_size = np.zeros(2)
        if self.frac_dir == "longitudal":
            array_well_size[0] = l_hor + 2 * self.hf
            array_well_size[1] = 0
        elif self.frac_dir == "transverse":
            array_well_size[0] = l_hor
            array_well_size[1] = 2 * self.hf

        return array_well_size

    def calc_multifrac(self, reservoir_model, S, xwd, ywd,
                       xed, yed, n_segments,
                       xbound, ybound,
                       calc_type, unit_length):
        """
        :param reservoir_model: тип пласта
        :param S: переменная пространства Лапласа
        :param xwd: безразмерная координата скважины
        :param ywd: безразмерная координата скважины
        :param xed: безразмерная координата границы x
        :param yed: безразмерная координата границы y
        :param n_segments: количество сегментов
        :param xbound: тип границы по x
        :param ybound: тип границы по y
        :param calc_type: тип расчёта
        :param unit_length: единица длины
        :return:
        """
        if reservoir_model != "composite":
            pd = self.multifrac_response_lapl(S, xwd, ywd,
                                              xed, yed, self.cfd, n_segments,
                                              xbound, ybound,
                                              calc_type)
        else:
            R_in_d = self.R_in / unit_length
            pd = self.multifrac_response_lapl_composite(S, xwd, ywd, xed, yed, self.cfd,
                                                        R_in_d, self.rel_M, self.rel_D, xbound, ybound
                                                        )

        return pd

    def multifrac_response_lapl(self, S, xwd, ywd, xed, yed, Fcd, n_segments, xbound, ybound,
                                calc_type):
        """
        :param S: переменная пространства Лапласа
        :param xwd: безразмерная координата скважины
        :param ywd: безразмерная координата скважины
        :param xed: безразмерная координата границы x
        :param yed: безразмерная координата границы y
        :param Fcd: безразмерный поправочный коэффициент на эффект ствола относительно полудлины трещины
        :param n_segments: количество сегментов
        :param xbound: тип границы по x
        :param ybound: тип границы по y
        :param calc_type: тип расчёта
        :return:
        """
        if calc_type == "segmentation":
            n_seg = n_segments // 2
            arr_Rates_and_Pressure, N_seg_total = self.pd_cinco_multifrac_rect(S, n_seg, xwd, ywd, xed, yed, xbound,
                                                                               ybound,
                                                                               Fcd, 1
                                                                               )
        else:
            arr_Rates_and_Pressure, N_seg_total = self.multifrac_response_lapl_desuperposition(S, xwd, ywd, xed, yed,
                                                                                               Fcd, xbound, ybound,
                                                                                               calc_type)

        multifrac_response_lapl = arr_Rates_and_Pressure[N_seg_total + 1]

        return multifrac_response_lapl

    def flow_resistance_asymmetrical(self, mm, N, cfd):
        """
        :param mm: матрица сопротивлений
        :param N: размерность вектора
        :param cfd: безразмерный поправочный коэффициент на эффект ствола относительно полудлины трещины
        :return: mm
        """
        dx = 1 / N
        a = dx / cfd
        for i in range(1, N + 1):
            for k in range(1, i):
                mm[i, k] = a * (N - (i - 1 / 2))
            mm[i, i] = a * (N - (i - 3 / 8))
            for k in range(1 + i, N + 1):
                mm[i, k] = a * (N - (k - 1 / 2))
        for i in range(1, N + 1):
            for k in range(1, i):
                mm[i + N, k + N] = a * (k - 1 / 2)
            mm[i + N, i + N] = a * (i - 5 / 8)
            for k in range(1 + i, N + 1):
                mm[i + N, k + N] = a * (i - 1 / 2)

        return mm

    def pd_cinco_multifrac_rect(self, S, n_seg, xwd, ywd, xed, yed, xbound,
                                ybound,
                                Fcd, Qd
                                ):
        """
        :param S: переменная пространства Лапласа
        :param n_seg: количество сегментов
        :param xwd: безразмерная координата скважины
        :param ywd: безразмерная координата скважины
        :param xed: безразмерная координата границы x
        :param yed: безразмерная координата границы y
        :param xbound: тип границы по x
        :param ybound: тип границы по y
        :param Fcd: безразмерный поправочный коэффициент на эффект ствола относительно полудлины трещины
        :param Qd: общий дебит скважины
        :return:
        """
        N_seg_total = self.n_frac * 2 * n_seg

        ql = self.cinco_multifrac_rect(S, n_seg, xwd, ywd, xed, yed, xbound, ybound, Fcd,
                                       Qd)

        pd_cinco_multifrac_rect = ql

        return pd_cinco_multifrac_rect, N_seg_total

    def cinco_multifrac_rect(self, S, n_seg, xwd, ywd, xed, yed, xbound, ybound, Fcd, Qd):
        """
        :param S: переменная пространства Лапласа
        :param n_seg: количество сегментов
        :param xwd: безразмерная координата скважины
        :param ywd: безразмерная координата скважины
        :param xed: безразмерная координата границы x
        :param yed: безразмерная координата границы y
        :param xbound: тип границы по x
        :param ybound: тип границы по y
        :param Fcd: безразмерный поправочный коэффициент на эффект ствола относительно полудлины трещины
        :param Qd: общий дебит скважины
        :return:
        """
        N_seg_total = self.n_frac * (2 * n_seg)
        flow_res_matr = np.zeros((2 * n_seg + 1, 2 * n_seg + 1))
        infl_matr = np.zeros((N_seg_total + 1, N_seg_total + 1))
        main_matr = np.zeros((N_seg_total + 2, N_seg_total + 2))
        flow_res_matr = self.flow_resistance_asymmetrical(flow_res_matr, n_seg, Fcd)
        rght = self.rhs(N_seg_total, Qd)
        infl_matr = self.influence_matrix_multifrac_rect(S, infl_matr, n_seg, xwd, ywd, xed, yed,
                                                         xbound, ybound)
        main_matr = self.main_matrix_multifrac(main_matr, flow_res_matr, infl_matr, n_seg)
        rght = self.LU_(main_matr, rght, N_seg_total + 1)

        cinco_multifrac_rect = rght

        return cinco_multifrac_rect

    def influence_matrix_multifrac_rect(self, S, infl_matr, n_seg, xwd, ywd, xed, yed,
                                        xbound, ybound):
        """
        :param S: переменная пространства Лапласа
        :param infl_matr: матрица влияния
        :param n_seg: количество сегментов
        :param xwd: безразмерная координата скважины
        :param ywd: безразмерная координата скважины
        :param xed: безразмерная координата границы x
        :param yed: безразмерная координата границы y
        :param xbound: тип границы по x
        :param ybound: тип границы по y
        :return:
        """
        dx = 1 / n_seg
        unit_length = dx / 2
        Sd = S * unit_length ** 2
        xedd = xed / unit_length
        yedd = yed / unit_length
        S_frac_dd = self.sf / unit_length
        for i in range(1, self.n_frac + 1):
            ydd = ywd[i] / unit_length
            for j in range(1, 2 * n_seg + 1):
                xdd = (xwd[i] + (-1 + dx / 2 + dx * (j - 1))) / unit_length
                point_index = (i - 1) * 2 * n_seg + j
                for l in range(1, self.n_frac + 1):
                    ywdd = ywd[l] / unit_length
                    for m in range(1, 2 * n_seg + 1):
                        xwdd = (xwd[l] + (-1 + dx / 2 + dx * (m - 1))) / unit_length
                        source_index = (l - 1) * 2 * n_seg + m
                        if Sd > 0:
                            infl_matr[point_index, source_index] = self.pd_lapl_rect(Sd, xdd, xwdd, xedd, ydd, ywdd,
                                                                                     yedd,
                                                                                     xbound, ybound, "frac")
                            if l == i:
                                infl_matr[point_index, source_index] = infl_matr[point_index,
                                                                                 source_index] + self.S_choke / (2 * pi)
                        else:
                            infl_matr[point_index, source_index] = self.pd_rect_BD(xdd, xwdd, xedd, ydd, ywdd, yedd,
                                                                                   xbound,
                                                                                   ybound, "frac")
                        if point_index == source_index:
                            infl_matr[point_index, source_index] = infl_matr[point_index,
                                                                             source_index] + S_frac_dd / (2 * pi)

        return infl_matr

    def main_matrix_multifrac(self, main_matr, flow_res_matr, infl_matr, n_seg):
        """
        :param main_matr: главная матрица
        :param flow_res_matr: матрица сопротивлений
        :param infl_matr: матрица влияния
        :param n_seg: количество сегментов
        :return:
        """
        N_seg_total = self.n_frac * (2 * n_seg)
        for i in range(1, N_seg_total + 1):
            for j in range(1, N_seg_total + 1):
                main_matr[i, j] = infl_matr[i, j]
        for i in range(1, self.n_frac + 1):
            modulus = (i - 1) * 2 * n_seg
            for l in range(1, 2 * n_seg + 1):
                i_maped = modulus + l
                for m in range(1, 2 * n_seg + 1):
                    j_maped = modulus + m
                    main_matr[i_maped, j_maped] = main_matr[i_maped, j_maped] + flow_res_matr[l, m]
        for i in range(1, N_seg_total + 1):
            main_matr[i, N_seg_total + 1] = -1
        for j in range(1, N_seg_total + 1):
            main_matr[N_seg_total + 1, j] = 1
        main_matr[N_seg_total + 1, N_seg_total + 1] = 0

        return main_matr

    def multifrac_response_lapl_desuperposition(self, S, xwd, ywd, xed, yed, cfd, xbound,
                                                ybound,
                                                calc_type):
        """
        :param S: переменная пространства Лапласа
        :param xwd: безразмерная координата скважины
        :param ywd: безразмерная координата скважины
        :param xed: безразмерная координата границы x
        :param yed: безразмерная координата границы y
        :param cfd: безразмерный поправочный коэффициент на эффект ствола относительно полудлины трещины
        :param xbound: тип границы по x
        :param ybound: тип границы по y
        :param calc_type: тип расчёта
        :return:
        """
        N_seg_total = self.n_frac
        matrix = np.zeros((N_seg_total + 2, N_seg_total + 2))
        bb = np.zeros(N_seg_total + 2)

        for j in range(1, self.n_frac + 1):
            for i in range(1, self.n_frac + 1):
                if calc_type == "optimal":
                    matrix[i, j] = self.PD_Frac_Rect_Lapl(S, xwd[i], ywd[i], xwd[j], ywd[j], xed,
                                                          yed, xbound, ybound, cfd, 10000000
                                                          )
                    if j == i:
                        matrix[i, j] += self.S_choke / (2 * pi)
                else:
                    raise ValueError("Данный тип расчета не поддерживается в модели ГС с МГРП. Введите оптимальный "
                                     "или сегментацию.")
        for i in range(1, N_seg_total + 1):
            matrix[i, N_seg_total + 1] = -1
            matrix[N_seg_total + 1, i] = 1
            bb[i] = 0

        matrix[N_seg_total + 1, N_seg_total + 1] = 0
        bb[N_seg_total + 1] = 1

        bb = self.LU_(matrix, bb, N_seg_total + 1)

        sumq = 0
        for i in range(1, N_seg_total + 1):
            sumq += bb[i]

        multifrac_response_lapl_desuperposition = bb

        return multifrac_response_lapl_desuperposition, N_seg_total

    def multifrac_response_lapl_composite(self, S, xwd, ywd, xed, yed, cfd, R_in_d, rel_M, rel_D,
                                          xbound, ybound):
        """
        :param S: переменная пространства Лапласа
        :param xwd: безразмерная координата скважины
        :param ywd: безразмерная координата скважины
        :param xed: безразмерная координата границы x
        :param yed: безразмерная координата границы y
        :param cfd: безразмерный поправочный коэффициент на эффект ствола относительно полудлины трещины
        :param R_in_d: безразмерный внутренний диаметр
        :param rel_M: коэффициент подвижности
        :param rel_D: коэффициент пьезопроводности
        :param xbound: тип границы по x
        :param ybound: тип границы по y
        :return:
        """
        N_seg_total = self.n_frac
        matrix = np.zeros((N_seg_total + 2, N_seg_total + 2))
        bb = np.zeros(N_seg_total + 2)
        for j in range(1, self.n_frac + 1):
            for i in range(1, self.n_frac + 1):
                matrix[i, j] = self.PD_Frac_Rect_Lapl_composite(S, xwd[i], ywd[i], xwd[j], ywd[j], xed, yed,
                                                                xbound, ybound, cfd, R_in_d, rel_M, rel_D,
                                                                10000000)
                if j == i:
                    matrix[i, j] = matrix[i, j] + self.S_choke / (2 * pi)
        for i in range(1, N_seg_total + 1):
            matrix[i, N_seg_total + 1] = -1
            matrix[N_seg_total + 1, i] = 1
            bb[i] = 0
        matrix[N_seg_total + 1, N_seg_total + 1] = 0
        bb[N_seg_total + 1] = 1
        bb = self.LU_(matrix, bb, N_seg_total + 1)
        multifrac_response_lapl_composite = bb[N_seg_total + 1]

        return multifrac_response_lapl_composite

    def PD_Frac_Rect_Lapl(self, S, xd, yd, xwd, ywd, xed, yed, xbound, ybound, cfd, etad
                          ):
        """
        :param S: переменная пространства Лапласа
        :param xd: безразмерная координата в направлении x
        :param yd: безразмерная координата в направлении y
        :param xwd: безразмерная координата скважины
        :param ywd: безразмерная координата скважины
        :param xed: безразмерная координата границы x
        :param yed: безразмерная координата границы y
        :param xbound: тип границы по x
        :param ybound: тип границы по y
        :param cfd: безразмерный поправочный коэффициент на эффект ствола относительно полудлины трещины
        :param etad: безразмерный коэффициент диффузии трещины
        :return:
        """
        compl_type = "frac"
        if cfd < self.tiny_2:
            Uniform_Flux = True
        else:
            Uniform_Flux = False
        if Uniform_Flux:
            dxd = 0
        else:
            dxd = self.calcxd(cfd)
        if (abs(xd - xwd) < 1) and (abs(yd - ywd) < self.tiny_2):
            if S > 0:
                pwd = self.pd_lapl_rect(S, xwd + dxd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type) + \
                      self.dpd_fincond_sf(S, cfd, self.sf, etad, Uniform_Flux) + self.dpd_dxd(S, abs(xd - xwd))
            else:
                pwd = self.pd_rect_BD(xwd + dxd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type) + \
                      self.dpd_fincond_sf_BD(cfd, self.sf, Uniform_Flux) + self.dpd_dxd_BD(abs(xd - xwd))
            if abs(xd - xwd) < self.tiny_2:
                pwd = pwd + self.S_choke / (2 * pi)
        else:
            if S > 0:
                pwd = self.pd_lapl_rect(S, xd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type)
            else:
                pwd = self.pd_rect_BD(xd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type)

        PD_Frac_Rect_Lapl = pwd

        return PD_Frac_Rect_Lapl

    def dpd_dxd(self, S, dxd):
        """
        :param S: переменная пространства Лапласа
        :param dxd: значение x_D, при котором модели трещины бесконечной проводимости и однородного потока совпадают
        :return:
        """
        dpd = (self.unit_fracture_func(S, dxd, 0) - self.unit_fracture_func(S, 0, 0))

        return dpd

    def dpd_dxd_BD(self, dxd):
        """
        :param dxd: значение x_D, при котором модели трещины бесконечной проводимости и однородного потока совпадают
        :return:
        """
        dpd = -1 / 2 * 1 / (2 * pi) * ((dxd - 1) * log((dxd + 1) / (1 - dxd)) + log((dxd + 1) ** 2))

        return dpd
