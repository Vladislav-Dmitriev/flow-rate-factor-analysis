from math import cos, sin, sqrt


class Layer:

    def __init__(self, k, h, ct, skin, phi, w, wf, l, lf, zw_h, xbound, ybound, zbound_up, zbound_down, kv_kh, rw, hw_f,
                 reservoir_model, unit_length, f_compressibility, f_porosity, lambda_, mgrp_flag, r=0):
        """
        :param k: проницаемость пласта
        :param h: эффективная толщина пласта
        :param ct: общая сжимаемость
        :param skin: скин-фактор
        :param phi: пористость пласта
        :param w: ширина прямоугольника
        :param wf: относительное расстояние от длинной стороны прямоугольника до скважины, доля ширины прямоугольника
        :param l: длина прямоугольника
        :param lf: относительное расстояние от короткой стороны прямоугольника до скважины, доля длины прямоугольника
        :param zw_h: вертикальное положение ствола в пласте (в долях от мощности) (для горизонтальных скважин)
        :param xbound: тип границы по x
        :param ybound: тип границы по y
        :param zbound_up: тип границы по z (верхняя)
        :param zbound_down: тип границы по z (нижняя)
        :param kv_kh: отношение вертикальной к горизонтальной проницаемости
        :param rw: радиус скважины
        :param hw_f: величина вскрытия пласта (в долях от мощности)
        :param reservoir_model: модель пласта
        :param unit_length: единица длины
        :param f_compressibility: сжимаемость трещин
        :param f_porosity: пористость трещин
        :param lambda_: параметр, характеризующий переток жидкости между матрицей и трещинами
        :param r: расстояние от центра скважины
        """
        self.k = k
        self.h = h
        self.ct = ct
        self.phi = phi
        self.skin = skin
        self.r = r
        self.w = w
        self.wf = wf
        self.l = l
        self.lf = lf
        self.zw_h = zw_h
        self.xbound = xbound
        self.ybound = ybound
        self.zbound_up = zbound_up
        self.zbound_down = zbound_down
        self.hw_f = hw_f
        self.reservoir_model = reservoir_model
        self.lambda_ = lambda_ * (unit_length / rw) ** 2
        if self.reservoir_model.startswith("DP"):
            self.omega = f_compressibility * f_porosity / (ct * phi)
        self.xd, self.yd, self.zd, self.hwd, self.xw, self.yw, self.zw, self.xe, self.ye, self.ze, self.s_kv_kh = self.calc_geometry(kv_kh, unit_length)
        if not mgrp_flag:
            self.xwd, self.ywd, self.zwd, self.xed, self.yed, self.zed, self.rwd = self.calc_dimensionless_geometry(self.xw, self.yw, self.zw, self.xe, self.ye, self.ze, rw, self.s_kv_kh,
                                             unit_length)

    def calc_geometry(self, kv_kh, unit_length, polar_angle=0):
        """
        :param kv_kh: отношение вертикальной к горизонтальной проницаемости
        :param unit_length: единица длины
        :param polar_angle: угол наклона
        :return: параметры геометрии пласта
        """
        xw = self.w * self.wf
        yw = self.l * self.lf
        xe = self.w
        ye = self.l
        ze = self.h
        zw = ze * self.zw_h
        X = self.r * cos(polar_angle) + xw
        Y = self.r * sin(polar_angle) + yw
        z = zw
        if kv_kh > 0:
            s_kv_kh = sqrt(kv_kh)
        else:
            s_kv_kh = 1
        hw = self.hw_f * self.h
        hwd = hw / unit_length / s_kv_kh
        xd = X / unit_length
        yd = Y / unit_length
        zd = z / unit_length / s_kv_kh
        return xd, yd, zd, hwd, xw, yw, zw, xe, ye, ze, s_kv_kh

    def calc_dimensionless_geometry(self, xw, yw, zw, xe, ye, ze, rw, s_kv_kh, unit_length):
        """
        :param xw: координата скважины
        :param yw: координата скважины
        :param zw: координата скважины
        :param xe: ширина прямоугольника
        :param ye: длина прямоугольника
        :param ze: координата границы по z
        :param rw: радиус скважины
        :param s_kv_kh: значение квадратного корня из kh_kv
        :param unit_length: единица длины
        :return: безразмерные параметры геометрии пласта
        """
        xwd = xw / unit_length
        ywd = yw / unit_length
        zwd = zw / unit_length / s_kv_kh
        xed = xe / unit_length
        yed = ye / unit_length
        zed = ze / unit_length / s_kv_kh
        rwd = rw * (1 + 1 / s_kv_kh) / 2 / unit_length

        return xwd, ywd, zwd, xed, yed, zed, rwd
