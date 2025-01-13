from math import pi


class Wellbore:

    def __init__(self, rw, mu, bl, l_hor, unit_length, h, phi, ct, wellbore_type, perf,
                 n_perf, well_permeability, layer_permeability):
        """
        :param rw: радиус скважины
        :param mu: вязкость
        :param bl: объемный коэффициент
        :param l_hor: длина горизонтального ствола
        :param unit_length: единица длины
        :param h: эффективная толщина пласта
        :param phi: пористость пласта
        :param ct: общая сжимаемость
        :param wellbore_type: тип скважины
        :param perf: доля перфораций горизонтальной скважины
        :param n_perf: количество перфорированных интервалов
        :param well_permeability: проницаемость скважины
        :param layer_permeability: проницаемость пласта
        """
        self.rw = rw
        self.mu = mu
        self.bl = bl
        self.l_hor = l_hor
        self.wellbore_type = wellbore_type
        self.perf = perf
        self.n_perf = n_perf
        self.well_permeability = well_permeability
        self.cd = self.calc_cd(unit_length, h, phi, ct)
        if self.wellbore_type == 'horizontal':
            self.fcd = (2 * self.well_permeability * 1000 * pi * self.rw ** 2) / (layer_permeability * self.l_hor * h)

    def calc_cd(self, unit_length, h, phi, ct, after_inflow_coef=0):
        """
        :param unit_length: единица длины
        :param h: эффективная толщина пласта
        :param phi: пористость пласта
        :param ct: общая сжимаемость
        :param after_inflow_coef: коэффициент послепритока
        :return: безразмерная константа эффекта ствола
        """
        cd = 1 / (2 * 3.14) / h / phi / ct / unit_length ** 2 * after_inflow_coef

        return cd

    def calc_pd_for_stable_model(self, S, pd, cd):
        """
        :param S: переменная пространства Лапласа
        :param pd: давление
        :param cd: безразмерная константа эффекта ствола
        :return: давление для модели послепритока Stable
        """
        pd = pd / (1 + (2 * pi) * cd * S * pd)

        return pd

    def well_fits_rectangle(self, xe, ye, array_well_size):
        """
        функция проверяет, помещается ли скважина в прямоугольную область дренирования
        :param xe: длина прямоугольника
        :param ye: ширина прямоугольника
        :param array_well_size: массив с размерами скважины
        :return: True или False
        """
        res = False

        if xe >= array_well_size[0] and ye >= array_well_size[1]:
            res = True
        well_fits_rectangle = res

        return well_fits_rectangle
