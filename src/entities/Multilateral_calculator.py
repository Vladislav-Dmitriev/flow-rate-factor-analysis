import numpy as np
from math import cos, sin, pi, sqrt, exp, tan
from scipy import special
from src.entities.Vert_pp_calculator import Vert_pp_Calculator
from src.entities import series_fn


class Multilateral_Calculator(Vert_pp_Calculator):
    """
    Класс расчёта параметров модели многозабойной скважины
    """

    def __init__(self, extend_reflections_mode, tiny, tiny_2, large_s, small, part_sum_num, max_it, small_bessel_arg,
                 part_sum_num_b2_1, part_sum_num_f_4, n_lateral, l_lateral, psi_lateral):
        super().__init__(extend_reflections_mode, tiny, tiny_2, large_s, small, part_sum_num, max_it, small_bessel_arg,
                         part_sum_num_b2_1, part_sum_num_f_4)
        self.n_segments_main_bore_multilat_cinco = 20
        self.n_segments_lateral_bore_multilat_cinco = 2
        self.n_segments_main_bore_multilat = 5
        self.n_segments_lateral_bore_multilat = 2
        self.main_bore_control_point_offset = 0
        self.multilat_first_branch_coeff = 0
        self.n_lateral = n_lateral
        self.l_lateral = l_lateral
        self.lateral_control_point_offset = 1
        self.psi_lateral = psi_lateral * pi / 180

    def calc_multilateral_geometry(self, l_hor):
        """
        :param l_hor: длина горизонтального ствола
        :return: массив с размерами скважины
        """
        array_well_size = np.zeros(2)
        psi_lateral_rad = self.psi_lateral * pi / 180
        if self.n_lateral > 0:
            array_well_size[0] = l_hor / 2 - l_hor / (2 * self.n_lateral) + self.l_lateral * cos(psi_lateral_rad)
            array_well_size[1] = 2 * self.l_lateral * sin(psi_lateral_rad)
        else:
            array_well_size[0] = l_hor
            array_well_size[1] = 0

        return array_well_size

    def multilateral_rect_lapl(self, S, xwd, ywd, zwd, xed, yed, zed, rwd,
                               xbound, ybound, zbound_up, zbound_down, calc_type, num_segments_main, unit_length):
        """
        :param S: переменная пространства Лапласа
        :param xwd: безразмерная координата скважины в направлении x
        :param ywd: безразмерная координата скважины в направлении y
        :param zwd: безразмерная координата скважины в направлении z 
        :param xed: безразмерная координата границы x
        :param yed: безразмерная координата границы y 
        :param zed: безразмерная координата границы z
        :param rwd: безразмерный радиус скважины
        :param xbound: тип границы по x 
        :param ybound: тип границы по y
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z
        :param calc_type: тип расчёта
        :param num_segments_main: число сегментов
        :param unit_length: единица длины
        :return:
        """
        zd = zwd + rwd
        l_lateral_d = self.l_lateral / unit_length
        if calc_type == "segmentation":
            if num_segments_main is None:
                n_segments_main = self.n_segments_main_bore_multilat_cinco
                n_segments_lat = self.n_segments_lateral_bore_multilat_cinco
            else:
                n_segments_main = num_segments_main
                n_segments_lat = n_segments_main // 2
        else:
            n_segments_main = self.n_segments_main_bore_multilat
            n_segments_lat = self.n_segments_lateral_bore_multilat

        if 2 * self.n_lateral > n_segments_main:
            n_segments_main = 2 * self.n_lateral

        pwd = self.pd_cinco_multilateral_well_lapl_rect(S, self.n_lateral, n_segments_main, n_segments_lat, zd, xwd, ywd,
                                                        zwd,
                                                        xed, yed, zed, l_lateral_d, xbound, ybound, zbound_up,
                                                        zbound_down, calc_type)

        multilateral_rect_lapl = pwd

        return multilateral_rect_lapl

    def pd_cinco_multilateral_well_lapl_rect(self, S, N_laterals, Seg_num_hor, Seg_num_lat, zd, xwd, ywd, zwd,
                                             xed, yed, zed, l_lateral_d, xbound, ybound, zbound_up,
                                             zbound_down, calc_type):
        """
        :param S: переменная пространства Лапласа
        :param N_laterals: количество боковых стволов
        :param Seg_num_hor: количество сегментов основного ствола
        :param Seg_num_lat: количество сегментов латералей
        :param zd: безразмерная координата в направлении z
        :param xwd: безразмерная координата скважины в направлении x
        :param ywd: безразмерная координата скважины в направлении y
        :param zwd: безразмерная координата скважины в направлении z
        :param xed: безразмерная координата границы x
        :param yed: безразмерная координата границы y 
        :param zed: безразмерная координата границы z
        :param l_lateral_d: безразмерная длина бокового ствола
        :param xbound: тип границы по x
        :param ybound: тип границы по y
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z
        :param calc_type: тип расчёта
        :return:
        """
        ql = self.cinco_multilateral_well_lapl_rect(S, N_laterals, Seg_num_hor, Seg_num_lat, zd, xwd, ywd, zwd, xed,
                                                    yed,
                                                    zed, l_lateral_d, 0, xbound, ybound, zbound_up,
                                                    zbound_down, calc_type)

        pd_cinco_multilateral_well_lapl_rect = ql[N_laterals * Seg_num_lat + Seg_num_hor + 1]

        return pd_cinco_multilateral_well_lapl_rect

    def cinco_multilateral_well_lapl_rect(self, S, N_laterals, Seg_num_hor, Seg_num_lat, zd, xwd, ywd, zwd, xed, yed,
                                          zed, l_lateral_d, theta_main, xbound, ybound, zbound_up,
                                          zbound_down,
                                          calc_type):
        """
        :param S: переменная пространства Лапласа
        :param N_laterals: количество боковых стволов
        :param Seg_num_hor: количество сегментов основного ствола
        :param Seg_num_lat: количество сегментов латералей
        :param zd: безразмерная координата в направлении z
        :param xwd: безразмерная координата скважины в направлении x
        :param ywd: безразмерная координата скважины в направлении y
        :param zwd: безразмерная координата скважины в направлении z
        :param xed: безразмерная координата границы x
        :param yed: безразмерная координата границы y 
        :param zed: безразмерная координата границы z
        :param l_lateral_d: безразмерная длина бокового ствола
        :param theta_main: абсолютный угол основного ствола
        :param xbound: тип границы по x
        :param ybound: тип границы по y
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z
        :param calc_type: тип расчёта
        :return:
        """
        num_seg_total = N_laterals * Seg_num_lat + Seg_num_hor
        infl_matr = np.zeros((num_seg_total + 2, num_seg_total + 2))
        rght = self.rhs(num_seg_total, 1)
        infl_matr = self.influence_matrix_multilateral_well_rect(S, N_laterals, Seg_num_hor, Seg_num_lat, infl_matr, zd,
                                                                 xwd, ywd, zwd,
                                                                 xed, yed, zed, l_lateral_d, theta_main,
                                                                 xbound, ybound, zbound_up, zbound_down, calc_type)
        rght = self.LU_(infl_matr, rght, num_seg_total + 1)
        cinco_multilateral_well_lapl_rect = rght

        return cinco_multilateral_well_lapl_rect

    def influence_matrix_multilateral_well_rect(self, S, N_laterals, Seg_num_hor, Seg_num_lat, infl_matr, zd, xwd, ywd,
                                                zwd, xed, yed, zed, l_lateral_d, theta_main, xbound,
                                                ybound, zbound_up, zbound_down, calc_type):
        """ 
        :param S: переменная пространства Лапласа
        :param N_laterals: количество боковых стволов
        :param Seg_num_hor: количество сегментов основного ствола
        :param Seg_num_lat: количество сегментов латералей
        :param infl_matr: матрица влияния
        :param zd: безразмерная координата в направлении z
        :param xwd: безразмерная координата скважины в направлении x
        :param ywd: безразмерная координата скважины в направлении y
        :param zwd: безразмерная координата скважины в направлении z
        :param xed: безразмерная координата границы x
        :param yed: безразмерная координата границы y 
        :param zed: безразмерная координата границы z
        :param l_lateral_d: безразмерная длина бокового ствола
        :param theta_main: абсолютный угол основного ствола
        :param xbound: тип границы по x
        :param ybound: тип границы по y
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z 
        :param calc_type: тип расчёта
        :return:
        """
        N_seg_total = N_laterals * Seg_num_lat + Seg_num_hor

        if calc_type != "segmentation":
            if Seg_num_hor == 1:
                dldd_hor = self.calcxd(10000)
            else:
                dldd_hor = self.main_bore_control_point_offset * self.calcxd(10000)
            dldd_lateral = self.lateral_control_point_offset * self.calcxd(10000)
        else:
            dldd_hor = 0
            dldd_lateral = 0
        Zen_ang_main = pi / 2
        Zen_ang_lat = pi / 2
        l_main_d = 2

        x_start = np.zeros(N_laterals + 1)
        y_start = np.zeros(N_laterals + 1)
        z_start = np.zeros(N_laterals + 1)
        x_start[0] = xwd - l_main_d / 2 * cos(theta_main) * sin(Zen_ang_main)
        y_start[0] = ywd - l_main_d / 2 * sin(theta_main) * sin(Zen_ang_main)
        z_start[0] = zwd - l_main_d / 2 * cos(Zen_ang_main)
        if N_laterals > 0:
            dlat_lat_xd = (l_main_d / N_laterals) * cos(theta_main) * sin(Zen_ang_main)
            dlat_lat_yd = (l_main_d / N_laterals) * sin(theta_main) * sin(Zen_ang_main)
            dlat_lat_zd = (l_main_d / N_laterals) * cos(Zen_ang_main)
        else:
            dlat_lat_xd = l_main_d
            dlat_lat_yd = 0
            dlat_lat_zd = 0
        for i in range(1, N_laterals + 1):
            coeff = (self.multilat_first_branch_coeff / 2 + (i - 1))
            x_start[i] = x_start[0] + dlat_lat_xd * coeff
            y_start[i] = y_start[0] + dlat_lat_yd * coeff
            z_start[i] = z_start[0] + dlat_lat_zd * coeff

        src_type_i = [0] * (N_laterals + 1)
        n_seg_i = np.zeros(N_laterals + 1)
        l_wellbore_i = np.zeros(N_laterals + 1)
        azmth_i = np.zeros(N_laterals + 1)
        znth_i = np.zeros(N_laterals + 1)
        offset_i = np.zeros(N_laterals + 1)

        src_type_i[0] = "hor"
        n_seg_i[0] = Seg_num_hor
        l_wellbore_i[0] = l_main_d
        azmth_i[0] = theta_main
        znth_i[0] = Zen_ang_main
        offset_i[0] = dldd_hor
        for i in range(1, N_laterals + 1):
            src_type_i[i] = "hor"
            n_seg_i[i] = int(Seg_num_lat)
            l_wellbore_i[i] = l_lateral_d
            offset_i[i] = dldd_lateral
            azmth_i[i] = theta_main + self.psi_lateral * (-1) ** (i - 1)
            znth_i[i] = Zen_ang_lat

        dzwd = zd - zwd
        offs = 0
        arr_src_type = [0] * (N_seg_total + 1)
        arr_src_x = np.zeros(N_seg_total + 1)
        arr_src_y = np.zeros(N_seg_total + 1)
        arr_src_z = np.zeros(N_seg_total + 1)
        arr_cp_x = np.zeros(N_seg_total + 1)
        arr_cp_y = np.zeros(N_seg_total + 1)
        arr_cp_z = np.zeros(N_seg_total + 1)
        arr_offs_d = np.zeros(N_seg_total + 1)
        arr_dl = np.zeros(N_seg_total + 1)
        arr_az_ang = np.zeros(N_seg_total + 1)
        arr_zen_ang = np.zeros(N_seg_total + 1)
        for i in range(0, N_laterals + 1):
            dld = l_wellbore_i[i] / n_seg_i[i]
            dxd = dld * cos(azmth_i[i]) * sin(znth_i[i])
            dyd = dld * sin(azmth_i[i]) * sin(znth_i[i])
            dzd = dld * cos(znth_i[i])
            for j in range(1, int(n_seg_i[i]) + 1):
                coeff = (1 / 2 + (j - 1))
                ind = int(offs + j)
                arr_src_type[ind] = src_type_i[i]
                arr_src_x[ind] = x_start[i] + dxd * coeff
                arr_src_y[ind] = y_start[i] + dyd * coeff
                arr_src_z[ind] = z_start[i] + dzd * coeff
                arr_cp_x[ind] = arr_src_x[ind]
                arr_cp_y[ind] = arr_src_y[ind]
                arr_cp_z[ind] = arr_src_z[ind] + dzwd
                arr_offs_d[ind] = offset_i[i]
                arr_dl[ind] = dld
                arr_az_ang[ind] = azmth_i[i]
                arr_zen_ang[ind] = znth_i[i]
            offs = offs + n_seg_i[i]
        if calc_type == "optimal":
            speed_up_method = "hor_as_frac_on_neigbours"
        else:
            speed_up_method = "none"

        infl_matr = self.influence_matrix_U(S, N_seg_total, infl_matr, arr_src_x, arr_src_y, arr_src_z, arr_cp_x,
                                            arr_cp_y, arr_cp_z,
                                            arr_dl, arr_az_ang, arr_zen_ang,
                                            arr_offs_d, arr_src_type,
                                            xed, yed, zed, xbound, ybound, zbound_up, zbound_down, speed_up_method)

        return infl_matr

    def influence_matrix_U(self, S, N_seg_total, infl_matr, arr_src_x, arr_src_y, arr_src_z, arr_cp_x, arr_cp_y,
                           arr_cp_z,
                           arr_dl, arr_az_ang, arr_zen_ang, arr_offs_d, arr_src_type, xed, yed, zed, xbound, ybound,
                           zbound_up, zbound_down, speed_up_method):
        """
        :param S: переменная пространства Лапласа
        :param N_seg_total: общее число сегментов
        :param infl_matr: матрица влияния
        :param arr_src_x: координаты середины источников
        :param arr_src_y: координаты середины источников
        :param arr_src_z: координаты середины источников
        :param arr_cp_x: координаты точек контроля
        :param arr_cp_y: координаты точек контроля
        :param arr_cp_z: координаты точек контроля
        :param arr_dl: массив безразмерных длин всех источников
        :param arr_az_ang: азимутальные углы источников
        :param arr_zen_ang: зенитные углы источников
        :param arr_offs_d: безразмерный отступ от точек контроля
        :param arr_src_type: тип источника
        :param xed: безразмерная координата границы x
        :param yed: безразмерная координата границы y
        :param zed: безразмерная координата границы z
        :param xbound: тип границы по x
        :param ybound: тип границы по y
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z
        :param speed_up_method: метод ускорения
        :return:
        """
        for i in range(1, N_seg_total + 1):
            for j in range(1, N_seg_total + 1):
                unit_length = arr_dl[j] / 2
                xedd = xed / unit_length
                yedd = yed / unit_length
                zedd = zed / unit_length
                Sdd = S * unit_length ** 2
                Ang = arr_az_ang[j]
                if i == j:
                    dxdd = arr_offs_d[i] * cos(arr_az_ang[i]) * sin(arr_zen_ang[i])
                    dydd = arr_offs_d[i] * sin(arr_az_ang[i]) * sin(arr_zen_ang[i])
                    dzdd = arr_offs_d[i] * cos(arr_zen_ang[i])
                else:
                    dxdd = 0
                    dydd = 0
                    dzdd = 0
                xdd = arr_cp_x[i] / unit_length + dxdd
                ydd = arr_cp_y[i] / unit_length + dydd
                zdd = arr_cp_z[i] / unit_length + dzdd
                xw_dd = arr_src_x[j] / unit_length
                yw_dd = arr_src_y[j] / unit_length
                zw_dd = arr_src_z[j] / unit_length
                if i == j:
                    compl_type = arr_src_type[j]
                elif speed_up_method == "none" or zbound_up == "c" or zbound_down == "c":
                    compl_type = arr_src_type[j]
                elif speed_up_method == "hor_as_frac_on_neigbours":
                    compl_type = "frac"
                else:
                    compl_type = arr_src_type[j]

                infl_matr[i, j] = self.deviated_well_for_cinco(Sdd, xdd, ydd, zdd, xw_dd, yw_dd, zw_dd,
                                                               xedd, yedd, zedd, Ang, xbound, ybound, zbound_up,
                                                               zbound_down, compl_type)
        for i in range(1, N_seg_total + 1):
            infl_matr[i, N_seg_total + 1] = -1
        for j in range(1, N_seg_total + 1):
            infl_matr[N_seg_total + 1, j] = 1

        infl_matr[N_seg_total + 1, N_seg_total + 1] = 0

        return infl_matr

    def deviated_well_for_cinco(self, S, xd, yd, zd, xwd, ywd, zwd, xed, yed, zed, Psi, xbound, ybound, zbound_up,
                                zbound_down, compl_type):
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
        :param Psi: азимутальный угол
        :param xbound: тип границы по x
        :param ybound: тип границы по y
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z
        :param compl_type: тип скважины
        :return:
        """
        Theta = self.NormalizeAngle(Psi)
        Ld = 1 / zed
        if zbound_up == "n" and zbound_down == "n":
            if S > 0:
                pwd = self.pd_b1_DW(S, yd, ywd, yed, xed, Theta, xbound, ybound, True) + self.pd_b2_DW(S, xd, xwd, xed,
                                                                                                       yd, ywd, yed,
                                                                                                       Theta, xbound,
                                                                                                       ybound, True) + \
                      self.pd_b3_DW(S, xd, xwd, xed, yd, ywd, Theta, xbound, ybound)
                if compl_type == "hor":
                    pwd = pwd + self.F_b1_DW(S, xed, yd, ywd, yed, zd, zwd, zed, Ld, Theta, zbound_up, zbound_down,
                                             xbound,
                                             ybound, True) + self.F_b2_DW(S, yd, ywd, yed, xd, xwd, xed, zd, zwd, zed,
                                                                          Ld, Theta, zbound_up, zbound_down, xbound,
                                                                          ybound, True) + self.F_b3_DW(S, xd, xwd, xed,
                                                                                                       zd, zwd, zed, yd,
                                                                                                       ywd, Ld, Theta,
                                                                                                       zbound_up,
                                                                                                       zbound_down,
                                                                                                       xbound)
            else:
                pwd = self.pd_b1_DW(-1, yd, ywd, yed, xed, Theta, xbound, ybound, True) + self.pd_b2_DW(0, xd, xwd, xed,
                                                                                                        yd, ywd, yed,
                                                                                                        Theta, xbound,
                                                                                                        ybound, False)
                if compl_type == "hor":
                    pwd = pwd + self.F_b1_DW(0, xed, yd, ywd, yed, zd, zwd, zed, Ld, Theta, zbound_up, zbound_down,
                                             xbound,
                                             ybound, True) + self.F_b2_DW(0, yd, ywd, yed, xd, xwd, xed, zd, zwd, zed,
                                                                          Ld, Theta, zbound_up, zbound_down, xbound,
                                                                          ybound, True) + self.F_b3_DW(0, xd, xwd, xed,
                                                                                                       zd, zwd, zed, yd,
                                                                                                       ywd, Ld, Theta,
                                                                                                       zbound_up,
                                                                                                       zbound_down,
                                                                                                       xbound)
        else:
            if S > 0:
                pwd = self.F_b1_DW(S, xed, yd, ywd, yed, zd, zwd, zed, Ld, Theta,
                                   zbound_up, zbound_down, xbound, ybound, True) + self.F_b2_DW(S, yd, ywd, yed, xd,
                                                                                                xwd, xed, zd, zwd, zed,
                                                                                                Ld, Theta,
                                                                                                zbound_up, zbound_down,
                                                                                                xbound, ybound,
                                                                                                True) + self.F_b3_DW(S,
                                                                                                                     xd,
                                                                                                                     xwd,
                                                                                                                     xed,
                                                                                                                     zd,
                                                                                                                     zwd,
                                                                                                                     zed,
                                                                                                                     yd,
                                                                                                                     ywd,
                                                                                                                     Ld,
                                                                                                                     Theta,
                                                                                                                     zbound_up,
                                                                                                                     zbound_down,
                                                                                                                     xbound)
            else:
                pwd = self.F_b1_DW(0, xed, yd, ywd, yed, zd, zwd, zed, Ld, Theta,
                                   zbound_up, zbound_down, xbound, ybound, True) + self.F_b2_DW(0, yd, ywd, yed, xd,
                                                                                                xwd, xed, zd, zwd, zed,
                                                                                                Ld, Theta,
                                                                                                zbound_up, zbound_down,
                                                                                                xbound, ybound,
                                                                                                True) + self.F_b3_DW(0,
                                                                                                                     xd,
                                                                                                                     xwd,
                                                                                                                     xed,
                                                                                                                     zd,
                                                                                                                     zwd,
                                                                                                                     zed,
                                                                                                                     yd,
                                                                                                                     ywd,
                                                                                                                     Ld,
                                                                                                                     Theta,
                                                                                                                     zbound_up,
                                                                                                                     zbound_down,
                                                                                                                     xbound)

        deviated_well_for_cinco = pwd

        return deviated_well_for_cinco

    def pd_b1_DW(self, S, yd, ywd, yed, xed, Theta, xbound, ybound, subtract_inf):
        """
        :param S: переменная пространства Лапласа
        :param yd: безразмерная координата в направлении y
        :param ywd: безразмерная координата скважины в направлении y
        :param yed: безразмерная координата границы y
        :param xed: безразмерная координата границы x
        :param Theta: азимутальный угол в диапазоне от 0 до Пи
        :param xbound: тип границы по x
        :param ybound: тип границы по y
        :param subtract_inf: параметр, определяющий, нужно ли считать часть уравнения
        :return:
        """
        if xbound == "c":
            return 0
        signum = 1
        if ybound == "c":
            signum = -1
        if S <= 0:
            if Theta < 0:
                Theta = Theta + 2 * pi
            if ybound == "n":
                if Theta <= pi:
                    if ywd + sin(Theta) > yd > ywd - sin(Theta):
                        pd = (6 * yd ** 2 - 6 * yd * yed + 4 * yed ** 2 - 6 * ywd * yed + 6 * ywd ** 2 -
                              3 * yed * (yd - ywd) ** 2 * 1 / sin(Theta) - 3 * yed * sin(Theta) + 2 * sin(
                                    Theta) ** 2) / (
                                     12 * xed * yed)
                    elif yd >= ywd + sin(Theta):
                        pd = (3 * yd ** 2 - 6 * yd * yed + 2 * yed ** 2 + 3 * ywd ** 2 + sin(Theta) ** 2) / (
                                6 * xed * yed)
                    else:
                        pd = (3 * yd ** 2 - 6 * ywd * yed + 2 * yed ** 2 + 3 * ywd ** 2 + sin(Theta) ** 2) / (
                                6 * xed * yed)
                else:
                    if ywd + sin(Theta) < yd < ywd - sin(Theta):
                        pd = (6 * yd ** 2 - 6 * yd * yed + 4 * yed ** 2 - 6 * ywd * yed + 6 * ywd ** 2 +
                              3 * yed * (yd - ywd) ** 2 * 1 / sin(Theta) + 3 * yed * sin(Theta) + 2 * sin(
                                    Theta) ** 2) / (
                                     12 * xed * yed)
                    elif yd >= ywd - sin(Theta):
                        pd = (3 * yd ** 2 - 6 * yd * yed + 2 * yed ** 2 + 3 * ywd ** 2 + sin(Theta) ** 2) / (
                                6 * xed * yed)
                    else:
                        pd = (3 * yd ** 2 - 6 * ywd * yed + 2 * yed ** 2 + 3 * ywd ** 2 + sin(Theta) ** 2) / (
                                6 * xed * yed)
            else:
                if Theta <= pi:
                    if ywd + sin(Theta) > yd > ywd - sin(Theta):
                        pd = -(-2 * yd * yed + 4 * yd * ywd - 2 * ywd * yed + yed * (yd - ywd) ** 2 / sin(
                            Theta) + yed * sin(Theta)) / (4 * xed * yed)
                    elif yd >= ywd + sin(Theta):
                        pd = (ywd * (yed - yd)) / (xed * yed)
                    else:
                        pd = (yd * (yed - ywd)) / (xed * yed)
                else:
                    if ywd + sin(Theta) < yd < ywd - sin(Theta):
                        pd = -(-2 * yd * yed + 4 * yd * ywd - 2 * ywd * yed - yed * (yd - ywd) ** 2 / sin(
                            Theta) - yed * sin(Theta)) / (4 * xed * yed)
                    elif yd >= ywd - sin(Theta):
                        pd = (ywd * (yed - yd)) / (xed * yed)
                    else:
                        pd = (yd * (yed - ywd)) / (xed * yed)
        else:
            u = sqrt(S)
            summ = self.sumexp(u, yed)
            pd = 1 / 4 / xed / sin(Theta) / S * (signum * self.int_exp_11(u, yd, ywd, yed, xed, Theta) +
                                                 self.int_exp_21(u, yd, ywd, yed, xed, Theta) +
                                                 signum * self.int_exp_31(u, yd, ywd, yed, xed, Theta) +
                                                 self.int_exp_41(u, yd, ywd, yed, xed, Theta)) * (1 + summ)

        pd_b1_DW = pd

        return pd_b1_DW

    def pd_b2_DW(self, S, xd, xwd, xed, yd, ywd, yed, Theta, xbound, ybound, subtract_inf):
        """
        :param S: переменная пространства Лапласа
        :param xd: безразмерная координата в направлении x
        :param xwd: безразмерная координата скважины в направлении x
        :param xed: безразмерная координата границы x
        :param yd: безразмерная координата в направлении y
        :param ywd: безразмерная координата скважины в направлении y
        :param yed: безразмерная координата границы y
        :param Theta: азимутальный угол в диапазоне от 0 до Пи
        :param xbound: тип границы по x
        :param ybound: тип границы по y
        :param subtract_inf: параметр, определяющий, нужно ли считать часть уравнения
        :return:
        """
        if Theta < 0:
            Theta = Theta + 2 * pi
        summa = series_fn.nsum(
            lambda k: self.pd_b2_DW_k(k, S, xd, xwd, xed, yd, ywd, yed, Theta, xbound, ybound, subtract_inf),
            1,
            self.part_sum_num,
        )

        return summa

    def pd_b3_DW(self, S, xd, xwd, xed, yd, ywd, Theta, xbound, ybound, r_wd=0.001):
        """
        :param S: переменная пространства Лапласа
        :param xd: безразмерная координата в направлении x
        :param xwd: безразмерная координата скважины в направлении x
        :param xed: безразмерная координата границы x
        :param yd: безразмерная координата в направлении y
        :param ywd: безразмерная координата скважины в направлении y
        :param Theta: азимутальный угол в диапазоне от 0 до Пи
        :param xbound: тип границы по x
        :param ybound: тип границы по y
        :param r_wd: безразмерный радиус скважины
        :return:
        """
        u = sqrt(S)
        signum = 1
        if xbound == "c":
            signum = -1
        int_m = 5

        summa = 1 / 4 / pi * (
                signum * self.int_k0_gauss_DW(int_m, u, (xd - (-xwd)) * cos(Theta) - (yd - ywd) * sin(Theta),
                                              0, (xd - (-xwd)) * sin(Theta) + (yd - ywd) * cos(Theta),
                                              0) + signum * self.int_k0_gauss_DW(int_m, u,
                                                                                 (xd - (2 * xed - xwd)) * cos(
                                                                                     Theta) - (yd - ywd) * sin(
                                                                                     Theta),
                                                                                 0, (xd - (2 * xed - xwd)) * sin(
                Theta) + (yd - ywd) * cos(Theta), 0) +
                self.int_k0_gauss_DW(int_m, u, (xd - (xwd - 2 * xed)) * cos(Theta) + (yd - ywd) * sin(Theta), 0,
                                     - (xd - (xwd - 2 * xed)) * sin(Theta) + (yd - ywd) * cos(Theta), 0) +
                self.int_k0_gauss_DW(int_m, u, (xd - (xwd + 2 * xed)) * cos(Theta) + (yd - ywd) * sin(Theta), 0,
                                     - (xd - (xwd + 2 * xed)) * sin(Theta) + (yd - ywd) * cos(Theta), 0) +
                signum * self.int_k0_gauss_DW(int_m, u,
                                              (xd - (-xwd - 2 * xed)) * cos(Theta) - (yd - ywd) * sin(Theta), 0,
                                              (xd - (-xwd - 2 * xed)) * sin(Theta) + (yd - ywd) * cos(Theta), 0))

        if self.distance_XY_well_axis(xd, yd, xwd, ywd, Theta) < r_wd:
            xd1 = abs(sqrt((xd - xwd) ** 2 + (yd - ywd) ** 2))
            summa += self.unit_fracture_func(S, xd1, 0)
        else:
            summa += + 1 / 4 / pi * self.int_k0_gauss_DW(int_m, u,
                                                              (xd - xwd) * cos(Theta) + (yd - ywd) * sin(Theta),
                                                              0, -(xd - xwd) * sin(Theta) + (yd - ywd) * cos(Theta), 0)
        k = 2
        yd1 = (yd - ywd) ** 2
        summa = summa + series_fn.nsum(
            lambda k: (
                    signum
                    * self.unit_cylinder_source(S, sqrt((xd + xwd + 2 * xed * k) ** 2 + yd1))
                    + self.unit_cylinder_source(S, sqrt((xd - xwd + 2 * xed * k) ** 2 + yd1))
                    + signum
                    * self.unit_cylinder_source(S, sqrt((xd + xwd - 2 * xed * k) ** 2 + yd1))
                    + self.unit_cylinder_source(S, sqrt((xd - xwd - 2 * xed * k) ** 2 + yd1))
            ), k, np.inf,
        )

        pd_b3_DW = summa

        return pd_b3_DW

    def F_b1_DW(self, S, xed, yd, ywd, yed, zd, zwd, zed, Ld, Theta, zbound_up, zbound_down, xbound, ybound,
                subtract_inf):
        """
        :param S: переменная пространства Лапласа
        :param xed: безразмерная координата границы x
        :param yd: безразмерная координата в направлении y
        :param ywd: безразмерная координата скважины в направлении y
        :param yed: безразмерная координата границы y
        :param zd: безразмерная координата в направлении z
        :param zwd: безразмерная координата скважины в направлении z
        :param zed: безразмерная координата границы z
        :param Ld: безразмерная длина
        :param Theta: азимутальный угол в диапазоне от 0 до Пи
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z
        :param xbound: тип границы по x
        :param ybound: тип границы по y
        :param subtract_inf: параметр, определяющий, нужно ли считать часть уравнения
        :return:
        """
        if xbound == "c":
            return 0
        summa = series_fn.nsum(
            lambda k: self.F_b1_DW_k(k, S, xed, yd, ywd, yed, zd, zwd, zed, Ld, Theta, zbound_up, zbound_down, xbound,
                                     ybound, subtract_inf),
            1,
            self.part_sum_num,
        )

        return summa

    def F_b1_DW_k(self, k, S, xed, yd, ywd, yed, zd, zwd, zed, Ld, Theta, zbound_up, zbound_down, xbound, ybound,
                  subtract_inf):
        """
        :param k: номер слагаемого
        :param S: переменная пространства Лапласа
        :param xed: безразмерная координата границы x
        :param yd: безразмерная координата в направлении y
        :param ywd: безразмерная координата скважины в направлении y
        :param yed: безразмерная координата границы y
        :param zd: безразмерная координата в направлении z
        :param zwd: безразмерная координата скважины в направлении z
        :param zed: безразмерная координата границы z
        :param Ld: безразмерная длина
        :param Theta: азимутальный угол в диапазоне от 0 до Пи
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z
        :param xbound: тип границы по x
        :param ybound: тип границы по y
        :param subtract_inf: параметр, определяющий, нужно ли считать часть уравнения
        :return:
        """
        if (zbound_up == "c" and zbound_down == "n") or (zbound_up == "n" and zbound_down == "c"):
            N = (2 * k - 1) / 2
        else:
            N = k
        if zbound_up == "c" and zbound_down == "c":
            part__1 = sin(N * pi * zd / zed) * sin(N * pi * zwd / zed)
        else:
            part__1 = cos(N * pi * zd / zed) * cos(N * pi * zwd / zed)

        e_n = sqrt(S + (N * pi * Ld) ** 2)

        F_b1_DW_k = 2 * part__1 * self.pd_b1_DW(e_n, yd, ywd, yed, xed, Theta, xbound, ybound, subtract_inf)

        return F_b1_DW_k

    def F_b2_DW(self, S, yd, ywd, yed, xd, xwd, xed, zd, zwd, zed, Ld, Theta,
                zbound_up, zbound_down, xbound, ybound, subtract_inf):
        """
        :param S: переменная пространства Лапласа
        :param yd: безразмерная координата в направлении y
        :param ywd: безразмерная координата скважины в направлении y
        :param yed: безразмерная координата границы y
        :param xd: безразмерная координата в направлении x
        :param xwd: безразмерная координата скважины в направлении x
        :param xed: безразмерная координата границы x
        :param zd: безразмерная координата в направлении z
        :param zwd: безразмерная координата скважины в направлении z
        :param zed: безразмерная координата границы z
        :param Ld: безразмерная длина
        :param Theta: азимутальный угол в диапазоне от 0 до Пи
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z
        :param xbound: тип границы по x
        :param ybound: тип границы по y
        :param subtract_inf: параметр, определяющий, нужно ли считать часть уравнения
        :return:
        """
        summa = series_fn.nsum(
            lambda k: self.F_b2_DW_k_3(k, S, yd, ywd, yed, xd, xwd, xed, zd, zwd, zed, Ld, Theta,
                                       zbound_up, zbound_down, xbound, ybound, subtract_inf),
            1,
            self.part_sum_num,
        )

        return summa

    def pd_b2_DW_k(self, k, S, xd, xwd, xed, yd, ywd, yed, Theta, xbound, ybound, subtract_inf):
        """
        :param k: номер слагаемого
        :param S: переменная пространства Лапласа
        :param xd: безразмерная координата в направлении x
        :param xwd: безразмерная координата скважины в направлении x
        :param xed: безразмерная координата границы x
        :param yd: безразмерная координата в направлении y
        :param ywd: безразмерная координата скважины в направлении y
        :param yed: безразмерная координата границы y
        :param Theta: азимутальный угол в диапазоне от 0 до Пи
        :param xbound: тип границы по x
        :param ybound: тип границы по y
        :param subtract_inf: параметр, определяющий, нужно ли считать часть уравнения
        :return:
        """
        if subtract_inf:
            sbtr_inf = 0
        else:
            sbtr_inf = 1
        signum = 1
        if ybound == "c":
            signum = -1
        ek = sqrt(S + k ** 2 * pi ** 2 / xed ** 2)
        smexp = self.sumexp(ek, yed)

        if xbound == "n":
            part__1 = cos(k * pi * xd / xed)
            part__2 = (signum * self.int_exp_12(k, ek, yd, ywd, yed, xwd, xed, Theta) +
                       self.int_exp_22(k, ek, yd, ywd, yed, xwd, xed, Theta) +
            signum * self.int_exp_32(k, ek, yd, ywd, yed, xwd, xed, Theta) +
            self.int_exp_42(k, ek, yd, ywd, yed, xwd, xed, Theta)) * (1 + smexp)
            if sbtr_inf != 0:
                if Theta <= pi:
                    if ywd + sin(Theta) > yd > ywd - sin(Theta):
                        part__2 = part__2 + sbtr_inf * self.int_exp_02(k, ek, yd, ywd, yed, xwd, xed, Theta)
                    elif yd >= ywd + sin(abs(Theta)):
                        part__2 = part__2 + sbtr_inf * self.int_exp_02___1(k, ek, yd, ywd, yed, xwd, xed, Theta)
                    else:
                        part__2 = part__2 + sbtr_inf * self.int_exp_02___2(k, ek, yd, ywd, yed, xwd, xed, Theta)
                else:
                    if ywd + sin(Theta) < yd < ywd - sin(Theta):
                        part__2 = part__2 + sbtr_inf * self.int_exp_02_pi(k, ek, yd, ywd, yed, xwd, xed, Theta)
                    elif yd >= ywd - sin(abs(Theta)):
                        part__2 = part__2 + sbtr_inf * self.int_exp_02___1(k, ek, yd, ywd, yed, xwd, xed, Theta)
                    else:
                        part__2 = part__2 + sbtr_inf * self.int_exp_02___2(k, ek, yd, ywd, yed, xwd, xed, Theta)
        elif xbound == "c":
            part__1 = sin(k * pi * xd / xed)
            part__2 = (signum * self.int_exp_12_1(k, ek, yd, ywd, yed, xwd, xed, Theta) +
                       self.int_exp_22_1(k, ek, yd, ywd, yed, xwd, xed, Theta) +
            signum * self.int_exp_32_1(k, ek, yd, ywd, yed, xwd, xed, Theta) +
            self.int_exp_42_1(k, ek, yd, ywd, yed, xwd, xed, Theta)) * (1 + smexp)

            if sbtr_inf != 0:
                if Theta <= pi:
                    if ywd + sin(Theta) > yd > ywd - sin(Theta):
                        part__2 = part__2 + sbtr_inf * self.int_exp_02_1(k, ek, yd, ywd, yed, xwd, xed, Theta)
                    elif yd >= ywd + sin(abs(Theta)):
                        part__2 = part__2 + sbtr_inf * self.int_exp_02_1___1(k, ek, yd, ywd, yed, xwd, xed, Theta)
                    else:
                        part__2 = part__2 + sbtr_inf * self.int_exp_02_1___2(k, ek, yd, ywd, yed, xwd, xed, Theta)
                else:
                    if ywd + sin(Theta) < yd < ywd - sin(Theta):
                        part__2 = part__2 + sbtr_inf * self.int_exp_02_1_pi(k, ek, yd, ywd, yed, xwd, xed, Theta)
                    elif yd >= ywd - sin(abs(Theta)):
                        part__2 = part__2 + sbtr_inf * self.int_exp_02_1___1(k, ek, yd, ywd, yed, xwd, xed, Theta)
                    else:
                        part__2 = part__2 + sbtr_inf * self.int_exp_02_1___2(k, ek, yd, ywd, yed, xwd, xed, Theta)
        else:
            part__1 = 0

        part__2 = part__2 / (k ** 2 * pi ** 2 * (cos(Theta)) ** 2 + ek ** 2 * xed ** 2 * (sin(Theta)) ** 2)

        pd_b2_DW_k = part__1 * part__2 * 1 / 2 / xed / ek
            
        return pd_b2_DW_k

    def F_b2_DW_k_3(self, p, S, yd, ywd, yed, xd, xwd, xed, zd, zwd, zed, Ld, Theta,
                    zbound_up, zbound_down, xbound, ybound, subtract_inf):
        """
        :param p: номер слагаемого
        :param S: переменная пространства Лапласа
        :param yd: безразмерная координата в направлении y
        :param ywd: безразмерная координата скважины в направлении y
        :param yed: безразмерная координата границы y
        :param xd: безразмерная координата в направлении x
        :param xwd: безразмерная координата скважины в направлении x
        :param xed: безразмерная координата границы x
        :param zd: безразмерная координата в направлении z
        :param zwd: безразмерная координата скважины в направлении z
        :param zed: безразмерная координата границы z
        :param Ld: безразмерная длина
        :param Theta: азимутальный угол в диапазоне от 0 до Пи
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z
        :param xbound: тип границы по x
        :param ybound: тип границы по y
        :param subtract_inf: параметр, определяющий, нужно ли считать часть уравнения
        :return:
        """
        if (zbound_up == "c" and zbound_down == "n") or (zbound_up == "n" and zbound_down == "c"):
            N = (2 * p - 1) / 2
        else:
            N = p
        if zbound_up == "c" and zbound_down == "c":
            part__1 = sin(N * pi * zd / zed) * sin(N * pi * zwd / zed)
        else:
            part__1 = cos(N * pi * zd / zed) * cos(N * pi * zwd / zed)

        F_b2_DW_k_3 = part__1 * self.F_b2_DW_k_2(p, S, yd, ywd, yed, xd, xwd, xed, zd, zwd, zed, Ld, Theta,
                                                 zbound_up, zbound_down, xbound, ybound, subtract_inf)

        return F_b2_DW_k_3

    def F_b2_DW_k_2(self, p, S, yd, ywd, yed, xd, xwd, xed, zd, zwd, zed, Ld, Theta,
                    zbound_up, zbound_down, xbound, ybound, subtract_inf):
        """
        :param p: номер слагаемого
        :param S: переменная пространства Лапласа
        :param yd: безразмерная координата в направлении y
        :param ywd: безразмерная координата скважины в направлении y
        :param yed: безразмерная координата границы y
        :param xd: безразмерная координата в направлении x
        :param xwd: безразмерная координата скважины в направлении x
        :param xed: безразмерная координата границы x
        :param zd: безразмерная координата в направлении z
        :param zwd: безразмерная координата скважины в направлении z
        :param zed: безразмерная координата границы z
        :param Ld: безразмерная длина
        :param Theta: азимутальный угол в диапазоне от 0 до Пи
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z
        :param xbound: тип границы по x
        :param ybound: тип границы по y
        :param subtract_inf: параметр, определяющий, нужно ли считать часть уравнения
        :return:
        """
        summa = series_fn.nsum(
            lambda k: self.F_b2_DW_k_1(p, k, S, yd, ywd, yed, xd, xwd, xed, zd, zwd, zed, Ld, Theta,
                                       zbound_up, zbound_down, xbound, ybound, subtract_inf),
            1,
            self.part_sum_num,
        )

        return summa

    def F_b2_DW_k_1(self, p, m, S, yd, ywd, yed, xd, xwd, xed, zd, zwd, zed, Ld, Theta,
                    zbound_up, zbound_down, xbound, ybound, subtract_inf):
        """
        :param p: номер слагаемого
        :param m: номер слагаемого
        :param S: переменная пространства Лапласа
        :param yd: безразмерная координата в направлении y
        :param ywd: безразмерная координата скважины в направлении y
        :param yed: безразмерная координата границы y
        :param xd: безразмерная координата в направлении x
        :param xwd: безразмерная координата скважины в направлении x
        :param xed: безразмерная координата границы x
        :param zd: безразмерная координата в направлении z
        :param zwd: безразмерная координата скважины в направлении z
        :param zed: безразмерная координата границы z
        :param Ld: безразмерная длина
        :param Theta: азимутальный угол в диапазоне от 0 до Пи
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z
        :param xbound: тип границы по x
        :param ybound: тип границы по y
        :param subtract_inf: параметр, определяющий, нужно ли считать часть уравнения
        :return:
        """
        if subtract_inf:
            sbtr_inf = 0
        else:
            sbtr_inf = 1
        k = m

        if (zbound_up == "c" and zbound_down == "n") or (zbound_up == "n" and zbound_down == "c"):
            N = (2 * p - 1) / 2
        else:
            N = p
        sign = 1
        if ybound == "c":
            sign = -1
        e_n_k = sqrt(S + (N * pi * Ld) ** 2 + k ** 2 * pi ** 2 / xed ** 2)
        smexp = self.sumexp(e_n_k, yed)

        if xbound == "n":
            part__1 = cos(k * pi * xd / xed) * cos(k * pi * xwd / xed)
            part__2 = (self.int_exp_12(k, e_n_k, yd, ywd, yed, xwd, xed, Theta) +
                       self.int_exp_22(k, e_n_k, yd, ywd, yed, xwd, xed, Theta) +
                       self.int_exp_32(k, e_n_k, yd, ywd, yed, xwd, xed, Theta) +
                       self.int_exp_42(k, e_n_k, yd, ywd, yed, xwd, xed, Theta)) * (1 + smexp)
            if sbtr_inf != 0:
                if ywd + sin(abs(Theta)) > yd > ywd - sin(abs(Theta)):
                    part__2 = part__2 + sbtr_inf * self.int_exp_02(k, e_n_k, yd, ywd, yed, xwd, xed, Theta)
                elif yd > ywd + sin(abs(Theta)):
                    part__2 = part__2 + sbtr_inf * self.int_exp_02___1(k, e_n_k, yd, ywd, yed, xwd, xed, Theta)
                else:
                    part__2 = part__2 + sbtr_inf * self.int_exp_02___2(k, e_n_k, yd, ywd, yed, xwd, xed, Theta)
        elif xbound == "c":
            part__1 = sin(k * pi * xd / xed) * sin(k * pi * xwd / xed)
            part__2 = (self.int_exp_12_1(k, e_n_k, yd, ywd, yed, xwd, xed, Theta) +
                       self.int_exp_22_1(k, e_n_k, yd, ywd, yed, xwd, xed, Theta) +
                       self.int_exp_32_1(k, e_n_k, yd, ywd, yed, xwd, xed, Theta) +
                       self.int_exp_42_1(k, e_n_k, yd, ywd, yed, xwd, xed, Theta)) * (1 + smexp)
            if sbtr_inf != 0:
                if ywd + sin(abs(Theta)) > yd > ywd - sin(abs(Theta)):
                    part__2 = part__2 + sbtr_inf * self.int_exp_02_1(k, e_n_k, yd, ywd, yed, xwd, xed, Theta)
                elif yd > ywd + sin(abs(Theta)):
                    part__2 = part__2 + sbtr_inf * self.int_exp_02_1___1(k, e_n_k, yd, ywd, yed, xwd, xed, Theta)
                else:
                    part__2 = part__2 + sbtr_inf * self.int_exp_02_1___2(k, e_n_k, yd, ywd, yed, xwd, xed, Theta)
        else:
            part__1 = 0

        part__2 = 1 / (k ** 2 * pi ** 2 * (cos(Theta)) ** 2 + e_n_k ** 2 * xed ** 2 * (sin(Theta)) ** 2) * part__2

        F_b2_DW_k_1 = part__1 * part__2 / xed / e_n_k

        return F_b2_DW_k_1

    def F_b3_DW(self, S, xd, xwd, xed, zd, zwd, zed, yd, ywd, Ld, Theta, zbound_up, zbound_down, xbound, r_wd=0.001):
        """
        :param S: переменная пространства Лапласа
        :param xd: безразмерная координата в направлении x
        :param xwd: безразмерная координата скважины в направлении x
        :param xed: безразмерная координата границы x
        :param zd: безразмерная координата в направлении z
        :param zwd: безразмерная координата скважины в направлении z
        :param zed: безразмерная координата границы z
        :param yd: безразмерная координата в направлении y
        :param ywd: безразмерная координата скважины в направлении y
        :param Ld: безразмерная длина
        :param Theta: азимутальный угол в диапазоне от 0 до Пи
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z
        :param xbound: тип границы по x
        :param r_wd: безразмерный радиус скважины
        :return:
        """
        summa = series_fn.nsum(
            lambda k: self.F_b3_DW_k(k, S, xd, xwd, xed, zd, zwd, zed, yd, ywd, Ld, Theta, zbound_up, zbound_down,
                                     xbound),
            1,
            self.part_sum_num,
        )
        if S == 0:
            S = -1
        if self.distance_XY_well_axis(xd, yd, xwd, ywd, Theta) < r_wd:
            summa += self.f(S, sqrt((xd - xwd) ** 2 + (yd - ywd) ** 2), 0, zd, zwd, zed, Ld, 0, 0, zbound_up,
                            zbound_down)

        return summa

    def F_b3_DW_k(self, k, S, xd, xwd, xed, zd, zwd, zed, yd, ywd, Ld, Theta, zbound_up, zbound_down, xbound,
                  r_wd=0.001):
        """
        :param k: номер слагаемого
        :param S: переменная пространства Лапласа
        :param xd: безразмерная координата в направлении x
        :param xwd: безразмерная координата скважины в направлении x
        :param xed: безразмерная координата границы x
        :param zd: безразмерная координата в направлении z
        :param zwd: безразмерная координата скважины в направлении z
        :param zed: безразмерная координата границы z
        :param yd: безразмерная координата в направлении y
        :param ywd: безразмерная координата скважины в направлении y
        :param Ld: безразмерная длина
        :param Theta: азимутальный угол в диапазоне от 0 до Пи
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z
        :param xbound: тип границы по x
        :param r_wd: безразмерный радиус скважины
        :return:
        """
        if (zbound_up == "c" and zbound_down == "n") or (zbound_up == "n" and zbound_down == "c"):
            N = (2 * k - 1) / 2
        else:
            N = k
        if zbound_up == "c" and zbound_down == "c":
            part__1 = sin(N * pi * zd / zed) * sin(N * pi * zwd / zed)
        else:
            part__1 = cos(N * pi * zd / zed) * cos(N * pi * zwd / zed)
        signum = 1
        if xbound == 'n':
            signum = -1
        int_m = 5

        e_n = sqrt(S + (N * pi * Ld) ** 2)

        summa = 1 / 4 / pi * (
                signum * self.int_k0_gauss_DW(int_m, e_n, (xd - (-xwd)) * cos(Theta) - (yd - ywd) * sin(Theta),
                                              0, (xd - (-xwd)) * sin(Theta) + (yd - ywd) * cos(Theta), 0) +
                signum * self.int_k0_gauss_DW(int_m, e_n,
                                              (xd - (2 * xed - xwd)) * cos(Theta) - (yd - ywd) * sin(Theta),
                                              0, (xd - (2 * xed - xwd)) * sin(Theta) + (yd - ywd) * cos(Theta), 0) +
                self.int_k0_gauss_DW(int_m, e_n, (xd - (xwd - 2 * xed)) * cos(Theta) + (yd - ywd) * sin(Theta), 0,
                                     -(xd - (xwd - 2 * xed)) * sin(Theta) + (yd - ywd) * cos(Theta), 0) +
                self.int_k0_gauss_DW(int_m, e_n, (xd - (xwd + 2 * xed)) * cos(Theta) + (yd - ywd) * sin(Theta), 0,
                                     -(xd - (xwd + 2 * xed)) * sin(Theta) + (yd - ywd) * cos(Theta), 0) +
                signum * self.int_k0_gauss_DW(int_m, e_n,
                                              (xd - (-xwd - 2 * xed)) * cos(Theta) - (yd - ywd) * sin(Theta), 0,
                                              (xd - (-xwd - 2 * xed)) * sin(Theta) + (yd - ywd) * cos(Theta), 0))

        if self.distance_XY_well_axis(xd, yd, xwd, ywd, Theta) >= r_wd:
            summa = summa + 1 / 4 / pi * self.int_k0_gauss_DW(int_m, e_n,
                                                              (xd - xwd) * cos(Theta) + (yd - ywd) * sin(Theta), 0,
                                                              - (xd - xwd) * sin(Theta) + (yd - ywd) * cos(Theta), 0)

        p = 2

        signum = 1
        if xbound == "c":
            signum = -1
        yd1 = (yd - ywd) ** 2
        summa = summa + series_fn.nsum(
            lambda k: (
                    signum
                    * self.unit_cylinder_source(e_n, sqrt((xd + xwd + 2 * xed * k) ** 2 + yd1))
                    + self.unit_cylinder_source(e_n, sqrt((xd - xwd + 2 * xed * k) ** 2 + yd1))
                    + signum
                    * self.unit_cylinder_source(e_n, sqrt((xd + xwd - 2 * xed * k) ** 2 + yd1))
                    + self.unit_cylinder_source(e_n, sqrt((xd - xwd - 2 * xed * k) ** 2 + yd1))
            ), p, np.inf,
        )

        F_b3_DW_k = 2 * part__1 * summa

        return F_b3_DW_k

    def NormalizeAngle(self, Theta):
        """
        :param Theta: азимутальный угол
        :return:
        """
        while Theta > pi:
            Theta = Theta - pi
        while Theta < 0:
            Theta = Theta + pi

        NormalizeAngle = Theta

        return NormalizeAngle

    def int_exp_11(self, u, yd, ywd, yed, xed, Theta):
        """
        :param u: квадратный корень из S
        :param yd: безразмерная координата в направлении y
        :param ywd: безразмерная координата скважины в направлении y
        :param yed: безразмерная координата границы y
        :param xed: безразмерная координата границы x
        :param Theta: азимутальный угол в диапазоне от 0 до Пи
        :return:
        """
        int_exp_11 = -exp(-u * (yd + ywd + sin(Theta))) + exp(-u * (yd + ywd - sin(Theta)))

        return int_exp_11

    def int_exp_21(self, u, yd, ywd, yed, xed, Theta):
        """
        :param u: квадратный корень из S
        :param yd: безразмерная координата в направлении y
        :param ywd: безразмерная координата скважины в направлении y
        :param yed: безразмерная координата границы y
        :param xed: безразмерная координата границы x
        :param Theta: азимутальный угол в диапазоне от 0 до Пи
        :return:
        """
        int_exp_21 = -exp(-u * (2 * yed + yd - ywd + sin(Theta))) + exp(-u * (2 * yed + yd - ywd - sin(Theta)))

        return int_exp_21

    def int_exp_31(self, u, yd, ywd, yed, xed, Theta):
        """
        :param u: квадратный корень из S
        :param yd: безразмерная координата в направлении y
        :param ywd: безразмерная координата скважины в направлении y
        :param yed: безразмерная координата границы y
        :param xed: безразмерная координата границы x
        :param Theta: азимутальный угол в диапазоне от 0 до Пи
        :return:
        """
        int_exp_31 = -exp(-u * (2 * yed - yd - ywd + sin(Theta))) + exp(-u * (2 * yed - yd - ywd - sin(Theta)))

        return int_exp_31

    def int_exp_41(self, u, yd, ywd, yed, xed, Theta):
        """
        :param u: квадратный корень из S
        :param yd: безразмерная координата в направлении y
        :param ywd: безразмерная координата скважины в направлении y
        :param yed: безразмерная координата границы y
        :param xed: безразмерная координата границы x
        :param Theta: азимутальный угол в диапазоне от 0 до Пи
        :return:
        """
        int_exp_41 = -exp(-u * (2 * yed - yd + ywd + sin(Theta))) + exp(-u * (2 * yed - yd + ywd - sin(Theta)))

        return int_exp_41

    def int_k0_gauss_DW(self, m, S, xd, xwd, yd, ywd):
        """
        :param m: индекс, который задает кусочно-линейную функцию, зависимую от x
        :param S: переменная пространства Лапласа
        :param xd: безразмерная координата в направлении x
        :param xwd: безразмерная координата скважины в направлении x
        :param yd: безразмерная координата в направлении y
        :param ywd: безразмерная координата скважины в направлении y
        :return:
        """
        N = 8
        a = -1
        b = 1
        if N == 6:
            Points_xi = np.array([-0.93247, -0.6612094, -0.2386142, 0.2386142, 0.6612094, 0.93247])
            Coeff_ci = np.array([0.1713245, 0.3607616, 0.467914, 0.467914, 0.3607616, 0.1713245])
        elif N == 8:
            Points_xi = np.array(
                [-0.96028986, -0.79666648, -0.52553242, -0.18343464, 0.18343464, 0.52553242, 0.79666648, 0.96028986])
            Coeff_ci = np.array(
                [0.10122854, 0.22238104, 0.31370664, 0.36268378, 0.36268378, 0.31370664, 0.22238104, 0.10122854])
        elif N == 5:
            Points_xi = np.array([-0.9061798, -0.5384693, 0, 0.5384693, 0.9061798])
            Coeff_ci = np.array([0.4786287, 0.2369269, 0.5688888, 0.2369269, 0.4786287])
        summa = 0
        for j in range(0, m):
            for i in range(0, N):
                xi = (a + (b - a) / m * j + a + (b - a) / m * (j + 1)) / 2 + (b - a) / m / 2 * Points_xi[i]
                summa += Coeff_ci[i] * special.kn(0, S * sqrt((xd - xwd - xi) ** 2 + (yd - ywd) ** 2))

        int_k0_gauss_DW = (b - a) * summa / 2 / m

        return int_k0_gauss_DW

    def distance_XY_well_axis(self, X_0, Y_0, X_w, Y_w, Theta):
        """
        :param X_0: безразмерная координата в направлении x
        :param Y_0: безразмерная координата в направлении y
        :param X_w: безразмерная координата скважины в направлении x
        :param Y_w: безразмерная координата скважины в направлении y
        :param Theta: азимутальный угол в диапазоне от 0 до Пи
        :return:
        """

        distance_XY_well_axis = abs(X_0 * sin(Theta) - Y_0 * cos(Theta) + (Y_w * cos(Theta) - X_w * sin(Theta)))

        return distance_XY_well_axis

    def int_exp_12(self, k, u, yd, ywd, yed, xwd, xed, Theta):
        """
        :param k: номер слагаемого
        :param u: квадратный корень из S
        :param yd: безразмерная координата в направлении y
        :param ywd: безразмерная координата скважины в направлении y
        :param yed: безразмерная координата границы y
        :param xwd: безразмерная координата скважины в направлении x
        :param xed: безразмерная координата границы x
        :param Theta: азимутальный угол в диапазоне от 0 до Пи
        :return:
        """
        int_exp_12 = -exp(-u * (yd + ywd + sin(Theta))) * \
                     xed * (u * xed * cos(k * pi * (xwd + cos(Theta)) / xed) * sin(Theta) -
                            k * pi * cos(Theta) * sin(k * pi * (xwd + cos(Theta)) / xed)) + \
                     exp(-u * (yd + ywd - sin(Theta))) * \
                     xed * (u * xed * cos(k * pi * (xwd - cos(Theta)) / xed) * sin(Theta) -
                            k * pi * cos(Theta) * sin(k * pi * (xwd - cos(Theta)) / xed))

        return int_exp_12

    def int_exp_22(self, k, u, yd, ywd, yed, xwd, xed, Theta):
        """
        :param k: номер слагаемого
        :param u: квадратный корень из S
        :param yd: безразмерная координата в направлении y
        :param ywd: безразмерная координата скважины в направлении y
        :param yed: безразмерная координата границы y
        :param xwd: безразмерная координата скважины в направлении x
        :param xed: безразмерная координата границы x
        :param Theta: азимутальный угол в диапазоне от 0 до Пи
        :return:
        """
        int_exp_22 = -exp(-u * (2 * yed + yd - ywd + sin(Theta))) * \
                     xed * (u * xed * cos(k * pi * (xwd - cos(Theta)) / xed) * sin(Theta) +
                            k * pi * cos(Theta) * sin(k * pi * (xwd - cos(Theta)) / xed)) + \
                     xed * exp(-u * (2 * yed + yd - ywd - sin(Theta))) * \
                     (u * xed * cos(k * pi * (xwd + cos(Theta)) / xed) * sin(Theta) +
                      k * pi * cos(Theta) * sin(k * pi * (xwd + cos(Theta)) / xed))

        return int_exp_22

    def int_exp_32(self, k, u, yd, ywd, yed, xwd, xed, Theta):
        """
        :param k: номер слагаемого
        :param u: квадратный корень из S
        :param yd: безразмерная координата в направлении y
        :param ywd: безразмерная координата скважины в направлении y
        :param yed: безразмерная координата границы y
        :param xwd: безразмерная координата скважины в направлении x
        :param xed: безразмерная координата границы x
        :param Theta: азимутальный угол в диапазоне от 0 до Пи
        :return:
        """
        int_exp_32 = -exp(-u * (2 * yed - yd - ywd + sin(Theta))) * \
                     xed * (u * xed * cos(k * pi * (xwd - cos(Theta)) / xed) * sin(Theta) +
                            k * pi * cos(Theta) * sin(k * pi * (xwd - cos(Theta)) / xed)) + \
                     xed * exp(-u * (2 * yed - yd - ywd - sin(Theta))) * \
                     (u * xed * cos(k * pi * (xwd + cos(Theta)) / xed) * sin(Theta) +
                      k * pi * cos(Theta) * sin(k * pi * (xwd + cos(Theta)) / xed))

        return int_exp_32

    def int_exp_42(self, k, u, yd, ywd, yed, xwd, xed, Theta):
        """
        :param k: номер слагаемого
        :param u: квадратный корень из S
        :param yd: безразмерная координата в направлении y
        :param ywd: безразмерная координата скважины в направлении y
        :param yed: безразмерная координата границы y
        :param xwd: безразмерная координата скважины в направлении x
        :param xed: безразмерная координата границы x
        :param Theta: азимутальный угол в диапазоне от 0 до Пи
        :return:
        """
        int_exp_42 = -exp(-u * (2 * yed - yd + ywd + sin(Theta))) * \
                     xed * (u * xed * cos(k * pi * (xwd + cos(Theta)) / xed) * sin(Theta) -
                            k * pi * cos(Theta) * sin(k * pi * (xwd + cos(Theta)) / xed)) + \
                     xed * exp(-u * (2 * yed - yd + ywd - sin(Theta))) * \
                     (u * xed * cos(k * pi * (xwd - cos(Theta)) / xed) * sin(Theta) -
                      k * pi * cos(Theta) * sin(k * pi * (xwd - cos(Theta)) / xed))

        return int_exp_42

    def int_exp_02(self, k, u, yd, ywd, yed, xwd, xed, Theta):
        """
        :param k: номер слагаемого
        :param u: квадратный корень из S
        :param yd: безразмерная координата в направлении y
        :param ywd: безразмерная координата скважины в направлении y
        :param yed: безразмерная координата границы y
        :param xwd: безразмерная координата скважины в направлении x
        :param xed: безразмерная координата границы x
        :param Theta: азимутальный угол в диапазоне от 0 до Пи
        :return:
        """
        int_exp_02 = exp(-u * (-yd + ywd + sin(Theta))) * \
                     xed * (-u * xed * cos(k * pi * (xwd + cos(Theta)) / xed) * sin(Theta) +
                            k * pi * cos(Theta) * sin(k * pi * (xwd + cos(Theta)) / xed)) - \
                     exp(-u * (yd - ywd + sin(Theta))) * \
                     xed * (u * xed * cos(k * pi * (xwd - cos(Theta)) / xed) * sin(Theta) +
                            k * pi * cos(Theta) * sin(k * pi * (xwd - cos(Theta)) / xed)) - \
                     xed * (-u * xed * cos(k * pi * (xwd + (yd - ywd) / tan(Theta)) / xed) * sin(Theta) +
                            k * pi * cos(Theta) * sin(k * pi * (xwd + (yd - ywd) / tan(Theta)) / xed)) + \
                     xed * (u * xed * cos(k * pi * (xwd + (yd - ywd) / tan(Theta)) / xed) * sin(Theta) +
                            k * pi * cos(Theta) * sin(k * pi * (xwd + (yd - ywd) / tan(Theta)) / xed))

        return int_exp_02

    def int_exp_02___1(self, k, u, yd, ywd, yed, xwd, xed, Theta):
        """
        :param k: номер слагаемого
        :param u: квадратный корень из S
        :param yd: безразмерная координата в направлении y
        :param ywd: безразмерная координата скважины в направлении y
        :param yed: безразмерная координата границы y
        :param xwd: безразмерная координата скважины в направлении x
        :param xed: безразмерная координата границы x
        :param Theta: азимутальный угол в диапазоне от 0 до Пи
        :return:
        """
        int_exp_02___1 = -exp(-u * (yd - ywd + sin(Theta))) * xed * (
                u * xed * cos(k * pi * (xwd - cos(Theta)) / xed) * sin(Theta) +
                k * pi * cos(Theta) * sin(k * pi * (xwd - cos(Theta)) / xed)) + \
                         exp(-u * (yd - ywd - sin(Theta))) * xed * (
                                 u * xed * cos(k * pi * (xwd + cos(Theta)) / xed) * sin(Theta) +
                                 k * pi * cos(Theta) * sin(k * pi * (xwd + cos(Theta)) / xed))

        return int_exp_02___1

    def int_exp_02___2(self, k, u, yd, ywd, yed, xwd, xed, Theta):
        """
        :param k: номер слагаемого
        :param u: квадратный корень из S
        :param yd: безразмерная координата в направлении y
        :param ywd: безразмерная координата скважины в направлении y
        :param yed: безразмерная координата границы y
        :param xwd: безразмерная координата скважины в направлении x
        :param xed: безразмерная координата границы x
        :param Theta: азимутальный угол в диапазоне от 0 до Пи
        :return:
        """
        int_exp_02___2 = exp(-u * (-yd + ywd + sin(Theta))) * xed * (
                -u * xed * cos(k * pi * (xwd + cos(Theta)) / xed) * sin(Theta) +
                k * pi * cos(Theta) * sin(k * pi * (xwd + cos(Theta)) / xed)) + \
                         exp(-u * (-yd + ywd - sin(Theta))) * xed * (
                                 u * xed * cos(k * pi * (xwd - cos(Theta)) / xed) * sin(Theta) -
                                 k * pi * cos(Theta) * sin(k * pi * (xwd - cos(Theta)) / xed))

        return int_exp_02___2

    def int_exp_12_1(self, k, u, yd, ywd, yed, xwd, xed, Theta):
        """
        :param k: номер слагаемого
        :param u: квадратный корень из S
        :param yd: безразмерная координата в направлении y
        :param ywd: безразмерная координата скважины в направлении y
        :param yed: безразмерная координата границы y
        :param xwd: безразмерная координата скважины в направлении x
        :param xed: безразмерная координата границы x
        :param Theta: азимутальный угол в диапазоне от 0 до Пи
        :return:
        """
        int_exp_12_1 = -exp(-u * (yd + ywd + sin(Theta))) * xed * (
                u * xed * sin(k * pi * (xwd + cos(Theta)) / xed) * sin(Theta) +
                k * pi * cos(Theta) * cos(k * pi * (xwd + cos(Theta)) / xed)) + exp(
            -u * (yd + ywd - sin(Theta))) * xed * (u * xed * sin(k * pi * (xwd - cos(Theta)) / xed) * sin(Theta) +
                                                   k * pi * cos(Theta) * cos(k * pi * (xwd - cos(Theta)) / xed))

        return int_exp_12_1

    def int_exp_22_1(self, k, u, yd, ywd, yed, xwd, xed, Theta):
        """
        :param k: номер слагаемого
        :param u: квадратный корень из S
        :param yd: безразмерная координата в направлении y
        :param ywd: безразмерная координата скважины в направлении y
        :param yed: безразмерная координата границы y
        :param xwd: безразмерная координата скважины в направлении x
        :param xed: безразмерная координата границы x
        :param Theta: азимутальный угол в диапазоне от 0 до Пи
        :return:
        """
        int_exp_22_1 = -exp(-u * (2 * yed + yd - ywd + sin(Theta))) * xed * (
                u * xed * sin(k * pi * (xwd - cos(Theta)) / xed) * sin(Theta) - k * pi * cos(Theta) * cos(
            k * pi * (xwd - cos(Theta)) / xed)) \
                       + xed * exp(-u * (2 * yed + yd - ywd - sin(Theta))) * (
                               u * xed * sin(k * pi * (xwd + cos(Theta)) / xed) * sin(Theta) -
                               k * pi * cos(Theta) * cos(k * pi * (xwd + cos(Theta)) / xed))

        return int_exp_22_1

    def int_exp_32_1(self, k, u, yd, ywd, yed, xwd, xed, Theta):
        """
        :param k: номер слагаемого
        :param u: квадратный корень из S
        :param yd: безразмерная координата в направлении y
        :param ywd: безразмерная координата скважины в направлении y
        :param yed: безразмерная координата границы y
        :param xwd: безразмерная координата скважины в направлении x
        :param xed: безразмерная координата границы x
        :param Theta: азимутальный угол в диапазоне от 0 до Пи
        :return:
        """
        int_exp_32_1 = -exp(-u * (2 * yed - yd - ywd + sin(Theta))) * xed * (
                u * xed * sin(k * pi * (xwd - cos(Theta)) / xed) * sin(Theta) -
                k * pi * cos(Theta) * cos(k * pi * (xwd - cos(Theta)) / xed)) + xed * exp(
            -u * (2 * yed - yd - ywd - sin(Theta))) * (u * xed * sin(k * pi * (xwd + cos(Theta)) / xed) * sin(Theta) -
                                                       k * pi * cos(Theta) * cos(k * pi * (xwd + cos(Theta)) / xed))

        return int_exp_32_1

    def int_exp_42_1(self, k, u, yd, ywd, yed, xwd, xed, Theta):
        """
        :param k: номер слагаемого
        :param u: квадратный корень из S
        :param yd: безразмерная координата в направлении y
        :param ywd: безразмерная координата скважины в направлении y
        :param yed: безразмерная координата границы y
        :param xwd: безразмерная координата скважины в направлении x
        :param xed: безразмерная координата границы x
        :param Theta: азимутальный угол в диапазоне от 0 до Пи
        :return:
        """
        int_exp_42_1 = -exp(-u * (2 * yed - yd + ywd + sin(Theta))) * xed * (
                u * xed * sin(k * pi * (xwd + cos(Theta)) / xed) * sin(Theta) +
                k * pi * cos(Theta) * cos(k * pi * (xwd + cos(Theta)) / xed)) + xed * exp(
            -u * (2 * yed - yd + ywd - sin(Theta))) * (u * xed * sin(k * pi * (xwd - cos(Theta)) / xed) * sin(Theta) +
                                                       k * pi * cos(Theta) * cos(k * pi * (xwd - cos(Theta)) / xed))

        return int_exp_42_1

    def int_exp_02_1(self, k, u, yd, ywd, yed, xwd, xed, Theta):
        """
        :param k: номер слагаемого
        :param u: квадратный корень из S
        :param yd: безразмерная координата в направлении y
        :param ywd: безразмерная координата скважины в направлении y
        :param yed: безразмерная координата границы y
        :param xwd: безразмерная координата скважины в направлении x
        :param xed: безразмерная координата границы x
        :param Theta: азимутальный угол в диапазоне от 0 до Пи
        :return:
        """
        int_exp_02_1 = -exp(-u * (yd - ywd + sin(Theta))) * xed * (
                u * xed * sin(Theta) * sin(k * pi * (xwd - cos(Theta)) / xed) - k * pi * cos(Theta) * cos(
            k * pi * (xwd - cos(Theta)) / xed)) - \
                       exp(-u * (-yd + ywd + sin(Theta))) * xed * (
                               u * xed * sin(Theta) * sin(k * pi * (xwd + cos(Theta)) / xed) +
                               k * pi * cos(Theta) * cos(k * pi * (xwd + cos(Theta)) / xed)) + xed * (
                               u * xed * sin(k * pi * (xwd + (yd - ywd) / tan(Theta)) / xed) * sin(Theta) -
                               k * pi * cos(Theta) * cos(k * pi * (xwd + (yd - ywd) / tan(Theta)) / xed)) + xed * (
                               u * xed * sin(k * pi * (xwd + (yd - ywd) / tan(Theta)) / xed) * sin(Theta) +
                               k * pi * cos(Theta) * cos(k * pi * (xwd + (yd - ywd) / tan(Theta)) / xed))

        return int_exp_02_1

    def int_exp_02_1___1(self, k, u, yd, ywd, yed, xwd, xed, Theta):
        """
        :param k: номер слагаемого
        :param u: квадратный корень из S
        :param yd: безразмерная координата в направлении y
        :param ywd: безразмерная координата скважины в направлении y
        :param yed: безразмерная координата границы y
        :param xwd: безразмерная координата скважины в направлении x
        :param xed: безразмерная координата границы x
        :param Theta: азимутальный угол в диапазоне от 0 до Пи
        :return:
        """
        int_exp_02_1___1 = -exp(-u * (yd - ywd + sin(Theta))) * xed * (
                u * xed * sin(k * pi * (xwd - cos(Theta)) / xed) * sin(Theta) -
                k * pi * cos(Theta) * cos(k * pi * (xwd - cos(Theta)) / xed)) + exp(
            -u * (yd - ywd - sin(Theta))) * xed * (u * xed * sin(k * pi * (xwd + cos(Theta)) / xed) * sin(Theta) -
                                                   k * pi * cos(Theta) * cos(k * pi * (xwd + cos(Theta)) / xed))

        return int_exp_02_1___1

    def int_exp_02_1___2(self, k, u, yd, ywd, yed, xwd, xed, Theta):
        """
        :param k: номер слагаемого
        :param u: квадратный корень из S
        :param yd: безразмерная координата в направлении y
        :param ywd: безразмерная координата скважины в направлении y
        :param yed: безразмерная координата границы y
        :param xwd: безразмерная координата скважины в направлении x
        :param xed: безразмерная координата границы x
        :param Theta: азимутальный угол в диапазоне от 0 до Пи
        :return:
        """
        int_exp_02_1___2 = exp(-u * (-yd + ywd - sin(Theta))) * xed * (
                u * xed * sin(k * pi * (xwd - cos(Theta)) / xed) * sin(Theta) +
                k * pi * cos(Theta) * cos(k * pi * (xwd - cos(Theta)) / xed)) - exp(
            -u * (-yd + ywd + sin(Theta))) * xed * (u * xed * sin(k * pi * (xwd + cos(Theta)) / xed) * sin(Theta) +
                                                    k * pi * cos(Theta) * cos(k * pi * (xwd + cos(Theta)) / xed))

        return int_exp_02_1___2

    def int_exp_02_pi(self, k, u, yd, ywd, yed, xwd, xed, Theta):
        """
        :param k: номер слагаемого
        :param u: квадратный корень из S
        :param yd: безразмерная координата в направлении y
        :param ywd: безразмерная координата скважины в направлении y
        :param yed: безразмерная координата границы y
        :param xwd: безразмерная координата скважины в направлении x
        :param xed: безразмерная координата границы x
        :param Theta: азимутальный угол в диапазоне от 0 до Пи
        :return:
        """
        int_exp_02_pi = exp(u * (-yd + ywd + sin(Theta))) * \
        xed * (u * xed * cos(k * pi * (xwd + cos(Theta)) / xed) * sin(Theta) +
               k * pi * cos(Theta) * sin(k * pi * (xwd + cos(Theta)) / xed)) + \
        exp(u * (yd - ywd + sin(Theta))) * \
        xed * (u * xed * cos(k * pi * (xwd - cos(Theta)) / xed) * sin(Theta) -
               k * pi * cos(Theta) * sin(k * pi * (xwd - cos(Theta)) / xed)) + \
        xed * (-u * xed * cos(k * pi * (xwd + (yd - ywd) / tan(Theta)) / xed) * sin(Theta) +
               k * pi * cos(Theta) * sin(k * pi * (xwd + (yd - ywd) / tan(Theta)) / xed)) + \
        xed * (-u * xed * cos(k * pi * (xwd + (yd - ywd) / tan(Theta)) / xed) * sin(Theta) -
               k * pi * cos(Theta) * sin(k * pi * (xwd + (yd - ywd) / tan(Theta)) / xed))

        return int_exp_02_pi

    def int_exp_02_1_pi(self, k, u, yd, ywd, yed, xwd, xed, Theta):
        """
        :param k: номер слагаемого
        :param u: квадратный корень из S
        :param yd: безразмерная координата в направлении y
        :param ywd: безразмерная координата скважины в направлении y
        :param yed: безразмерная координата границы y
        :param xwd: безразмерная координата скважины в направлении x
        :param xed: безразмерная координата границы x
        :param Theta: азимутальный угол в диапазоне от 0 до Пи
        :return:
        """
        int_exp_02_1_pi = exp(u * (yd - ywd + sin(Theta))) * xed * \
        (u * xed * sin(Theta) * sin(k * pi * (xwd - cos(Theta)) / xed) +
         k * pi * cos(Theta) * cos(k * pi * (xwd - cos(Theta)) / xed)) + \
        exp(u * (-yd + ywd + sin(Theta))) * xed * \
        (u * xed * sin(Theta) * sin(k * pi * (xwd + cos(Theta)) / xed) -
         k * pi * cos(Theta) * cos(k * pi * (xwd + cos(Theta)) / xed)) + \
        xed * (-u * xed * sin(k * pi * (xwd + (yd - ywd) / tan(Theta)) / xed) * sin(Theta) +
               k * pi * cos(Theta) * cos(k * pi * (xwd + (yd - ywd) / tan(Theta)) / xed)) + \
        xed * (-u * xed * sin(k * pi * (xwd + (yd - ywd) / tan(Theta)) / xed) * sin(Theta) -
               k * pi * cos(Theta) * cos(k * pi * (xwd + (yd - ywd) / tan(Theta)) / xed))

        return int_exp_02_1_pi
