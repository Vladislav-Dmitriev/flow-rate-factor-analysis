import numpy as np
from math import pi, exp, sin, cos, sqrt, log
from src.entities.Laplace_calculator import LaplaceCalculator
from src.entities import series_fn


class Vert_pp_Calculator(LaplaceCalculator):

    def partial_penetration_vert_rect_lapl(self, S, xd, yd, zd, xwd, ywd, zwd, xed, yed, zed, hwd, xbound, ybound,
                                           zbound_up, zbound_down, rwd, calc_type):
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
        :param hwd: безразмерная глубина скважины безразмерная глубина скважины
        :param xbound: тип границы 
        :param ybound: тип границы 
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z 
        :param rwd: безразмерный радиус скважины 
        :param calc_type: тип расчета 
        :return: 
        """
        compl_type = "vert"
        Large_Cfd = 10000
        Ld = 1 / zed
        dzd = self.calcxd(Large_Cfd)
        zd = zwd + dzd * hwd / 2
        zd_1 = zwd + dzd * hwd / 2
        zd_2 = zwd - dzd * hwd / 2
        if calc_type == 'optimal':
            if zbound_up == "n" and zbound_down == "n":
                if S > 0:
                    if zwd / zed == 0.5 or dzd == 0:
                        pwd = self.pd_lapl_rect(S, xd, xwd, xed, yd, ywd, yed, xbound, ybound,
                                                compl_type) + self.sum_pp_1(S, xd, xwd, xed, yd, ywd, yed, zd, zwd, zed,
                                                                            hwd, Ld, xbound, ybound, zbound_up,
                                                                            zbound_down, compl_type)
                    elif zwd / zed != 0.5 and dzd != 0:
                        pwd = self.pd_lapl_rect(S, xd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type) + (
                                self.sum_pp_1(S, xd, xwd, xed, yd, ywd, yed, zd_1, zwd, zed, hwd, Ld, xbound, ybound,
                                              zbound_up, zbound_down, compl_type) + self.sum_pp_1(S, xd, xwd, xed, yd,
                                                                                                  ywd, yed, zd_2, zwd,
                                                                                                  zed, hwd, Ld, xbound,
                                                                                                  ybound, zbound_up,
                                                                                                  zbound_down,
                                                                                                  compl_type)) / 2
                else:
                    if zwd / zed == 0.5 or dzd == 0:
                        pwd = self.pd_rect_BD(xd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type) + self.sum_pp(0,
                                                                                                                    xd,
                                                                                                                    xwd,
                                                                                                                    xed,
                                                                                                                    yd,
                                                                                                                    ywd,
                                                                                                                    yed,
                                                                                                                    zd,
                                                                                                                    zwd,
                                                                                                                    zed,
                                                                                                                    hwd,
                                                                                                                    Ld,
                                                                                                                    xbound,
                                                                                                                    ybound,
                                                                                                                    zbound_up,
                                                                                                                    zbound_down,
                                                                                                                    compl_type)
                    elif zwd / zed != 0.5 and dzd != 0:
                        pwd = self.pd_rect_BD(xd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type) + (
                                self.sum_pp(0, xd, xwd, xed, yd, ywd, yed, zd_1, zwd, zed, hwd, Ld, xbound, ybound,
                                            zbound_up,
                                            zbound_down, compl_type) +
                                self.sum_pp(0, xd, xwd, xed, yd, ywd, yed, zd_2, zwd, zed, hwd, Ld, xbound, ybound,
                                            zbound_up, zbound_down, compl_type)) / 2
            elif zbound_up == "c" and zbound_down == "c":
                if S > 0:
                    if zwd / zed == 0.5 or dzd == 0:
                        pwd = self.sum_pp_1(S, xd, xwd, xed, yd, ywd, yed, zd, zwd, zed, hwd, Ld, xbound, ybound,
                                            zbound_up,
                                            zbound_down, compl_type)
                    elif zwd / zed != 0.5 and dzd != 0:
                        pwd = (self.sum_pp_1(S, xd, xwd, xed, yd, ywd, yed, zd_1, zwd, zed, hwd, Ld, xbound, ybound,
                                             zbound_up, zbound_down, compl_type) +
                               self.sum_pp_1(S, xd, xwd, xed, yd, ywd, yed, zd_2, zwd, zed, hwd, Ld, xbound, ybound,
                                             zbound_up, zbound_down, compl_type)) / 2
                else:
                    if zwd / zed == 0.5 or dzd == 0:
                        pwd = self.sum_pp(0, xd, xwd, xed, yd, ywd, yed, zd, zwd, zed, hwd, Ld, xbound, ybound,
                                          zbound_up,
                                          zbound_down, compl_type)
                    elif zwd / zed != 0.5 and dzd != 0:
                        pwd = (self.sum_pp(0, xd, xwd, xed, yd, ywd, yed, zd_1, zwd, zed, hwd, Ld, xbound, ybound,
                                           zbound_up,
                                           zbound_down, compl_type) + self.sum_pp(0, xd, xwd, xed, yd, ywd, yed, zd_2,
                                                                                  zwd, zed, hwd, Ld, xbound, ybound,
                                                                                  zbound_up, zbound_down,
                                                                                  compl_type)) / 2
            else:
                if S > 0:
                    pwd = (self.sum_pp_1(S, xd, xwd, xed, yd, ywd, yed, zd_1, zwd, zed, hwd, Ld, xbound, ybound,
                                         zbound_up,
                                         zbound_down, compl_type) + self.sum_pp_1(S, xd, xwd, xed, yd, ywd, yed, zd_2,
                                                                                  zwd, zed, hwd, Ld, xbound, ybound,
                                                                                  zbound_up, zbound_down,
                                                                                  compl_type)) / 2
                else:
                    pwd = (self.sum_pp(0, xd, xwd, xed, yd, ywd, yed, zd_1, zwd, zed, hwd, Ld, xbound, ybound,
                                       zbound_up,
                                       zbound_down, compl_type) +
                           self.sum_pp(0, xd, xwd, xed, yd, ywd, yed, zd_2, zwd, zed, hwd, Ld, xbound, ybound,
                                       zbound_up, zbound_down, compl_type)) / 2
        elif calc_type == 'segmentation':
            Ld = 1 / zed
            n_seg = 5
            pwd = self.pw_solve_pp(S, n_seg, xd, yd, xwd, ywd, zwd, xed, yed, zed, rwd, hwd, Ld, xbound, ybound,
                                   zbound_up,
                                   zbound_down)
        else:
            raise ValueError("Данный тип расчета не поддерживается в модели. Введите оптимальный "
                             "или сегментацию.")

        partial_penetration_vert_rect_lapl = pwd

        return partial_penetration_vert_rect_lapl

    def calcxd(self, cfdx):
        """
        :param cfdx: безразмерный поправочный коэффициент на эффект ствола относительно полудлины трещины
        :return: 
        """
        if cfdx < self.tiny_2:
            calcxd = 0
            return calcxd
        a0 = 0.759919
        a1 = 0.465301
        a2 = 0.562754
        a3 = 0.363093
        a4 = 0.0298881
        b0 = 1
        b1 = 0.99477
        b2 = 0.896679
        b3 = 0.430707
        b4 = 0.0467339
        var = log(cfdx)
        calcxd = (a0 + a1 * var + a2 * var ** 2 + a3 * var ** 3 + a4 * var ** 4) / (
                b0 + b1 * var + b2 * var ** 2 + b3 * var ** 3 + b4 * var ** 4)

        return calcxd

    def sum_pp_1(self, S, xd, xwd, xed, yd, ywd, yed, zd, zwd, zed, hwd, Ld, xbound, ybound,
                 zbound_up, zbound_down, compl_type):
        """
        :param S: переменная пространства Лапласа 
        :param xd: безразмерная координата точки, в которой рассчитывается давление
        :param xwd: безразмерная координата скважины
        :param xed: безразмерная координата границы x
        :param yd: безразмерная координата точки, в которой рассчитывается давление
        :param ywd: безразмерная координата скважины
        :param yed: безразмерная координата границы y
        :param zd: безразмерная координата точки, в которой рассчитывается давление
        :param zwd: безразмерная координата скважины
        :param zed: безразмерная координата границы z
        :param hwd: безразмерная глубина скважины
        :param Ld: безразмерная длина 
        :param xbound: тип границы
        :param ybound: тип границы
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z
        :param compl_type: тип скважины
        :return: 
        """
        sum_pp_1 = self.pd_b1_pp(S, yd, ywd, yed, xed, xbound, ybound, Ld, zd, zwd, zed, hwd, zbound_up,
                                 zbound_down) + self.pd_b2_pp(S, xd, xwd, xed, yd, ywd, yed, xbound, ybound, Ld, zd,
                                                              zwd, zed, hwd, zbound_up, zbound_down,
                                                              compl_type) + self.pd_b3_pp(S, xd, xwd, xed, yd, ywd,
                                                                                          xbound, ybound, Ld, zd, zwd,
                                                                                          zed, hwd, zbound_up,
                                                                                          zbound_down, compl_type)

        return sum_pp_1

    def sum_pp(self, S, xd, xwd, xed, yd, ywd, yed, zd, zwd, zed, hwd, Ld, xbound, ybound, zbound_up, zbound_down,
               compl_type):
        """
        :param S: переменная пространства Лапласа 
        :param xd: безразмерная координата точки, в которой рассчитывается давление
        :param xwd: безразмерная координата скважины 
        :param xed: безразмерная координата границы x 
        :param yd: безразмерная координата точки, в которой рассчитывается давление
        :param ywd: безразмерная координата скважины 
        :param yed: безразмерная координата границы y 
        :param zd: безразмерная координата точки, в которой рассчитывается давление
        :param zwd: безразмерная координата скважины 
        :param zed: безразмерная координата границы z 
        :param hwd: безразмерная глубина скважины 
        :param Ld: безразмерная длина 
        :param xbound: тип границы 
        :param ybound: тип границы 
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z 
        :param compl_type: тип скважины
        :return: 
        """
        PART_SUM_NUM = 1500
        MAXIT = 500
        summa = series_fn.nsum(
            lambda k: self.sum_pp_k(k, S, xd, xwd, xed, yd, ywd, yed, zd, zwd, zed, hwd, Ld, xbound, ybound, zbound_up,
                                    zbound_down, compl_type),
            1,
            PART_SUM_NUM,
        )
        return summa

    def sum_pp_k(self, N, S, xd, xwd, xed, yd, ywd, yed, zd, zwd, zed, hwd, Ld, xbound, ybound, zbound_up, zbound_down,
                 compl_type):
        """
        :param N: номер слагаемого
        :param S: переменная пространства Лапласа 
        :param xd: безразмерная координата точки, в которой рассчитывается давление
        :param xwd: безразмерная координата скважины 
        :param xed: безразмерная координата границы x 
        :param yd: безразмерная координата точки, в которой рассчитывается давление
        :param ywd: безразмерная координата скважины 
        :param yed: безразмерная координата границы y 
        :param zd: безразмерная координата точки, в которой рассчитывается давление
        :param zwd: безразмерная координата скважины 
        :param zed: безразмерная координата границы z 
        :param hwd: безразмерная глубина скважины 
        :param Ld: безразмерная длина 
        :param xbound: тип границы 
        :param ybound: тип границы 
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z 
        :param compl_type: тип скважины
        :return: 
        """
        if (zbound_up == "c" and zbound_down == "n") or (zbound_up == "n" and zbound_down == "c"):
            p = (2 * N - 1) / 2
        else:
            p = N
        mult1 = self.pd_lapl_rect(S + (p * pi * Ld) ** 2, xd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type)
        if zbound_up == "c" and zbound_down == "c":
            part__1 = sin(p * pi * zd / zed) * sin(p * pi * zwd / zed)
        else:
            part__1 = cos(p * pi * zd / zed) * cos(p * pi * zwd / zed)
        sum_pp_k = part__1 * sin(p * pi * hwd / zed / 2) * mult1 / p * zed / hwd * 4 / pi

        return sum_pp_k

    def pd_b1_pp(self, S, yd, ywd, yed, xed, xbound, ybound, Ld, zd, zwd, zed, hwd, zbound_up, zbound_down):
        """
        :param S: переменная пространства Лапласа 
        :param yd: безразмерная координата точки, в которой рассчитывается давление
        :param ywd: безразмерная координата скважины 
        :param yed: безразмерная координата границы y 
        :param xed: безразмерная координата границы x 
        :param xbound: тип границы 
        :param ybound: тип границы 
        :param Ld: безразмерная длина 
        :param zd: безразмерная координата точки, в которой рассчитывается давление
        :param zwd: безразмерная координата скважины 
        :param zed: безразмерная координата границы z 
        :param hwd: безразмерная глубина скважины 
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z 
        :return: 
        """
        MAXIT = 500
        summa = series_fn.nsum(
            lambda k: self.pd_b1_pp_n(k, S, yd, ywd, yed, xed, xbound, ybound, Ld, zd, zwd, zed, hwd, zbound_up,
                                      zbound_down),
            1,
            MAXIT - 1,
        )

        return summa

    def pd_b1_pp_n(self, N, S, yd, ywd, yed, xed, xbound, ybound, Ld, zd, zwd, zed, hwd, zbound_up, zbound_down):
        """
        :param N: номер слагаемого
        :param S: переменная пространства Лапласа 
        :param yd: безразмерная координата точки, в которой рассчитывается давление
        :param ywd: безразмерная координата скважины 
        :param yed: безразмерная координата границы y 
        :param xed: безразмерная координата границы x 
        :param xbound: тип границы 
        :param ybound: тип границы 
        :param Ld: безразмерная длина 
        :param zd: безразмерная координата точки, в которой рассчитывается давление
        :param zwd: безразмерная координата скважины 
        :param zed: безразмерная координата границы z 
        :param hwd: безразмерная глубина скважины 
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z 
        :return: 
        """
        small = 1E-300
        if xbound == "c":
            pd_b1_pp_n = 0
            return pd_b1_pp_n
        if ybound == "n":
            signum = 1
        elif ybound == "c":
            signum = -1
        else:
            signum = 1
        if (zbound_up == "c" and zbound_down == "n") or (zbound_up == "n" and zbound_down == "c"):
            p = (2 * N - 1) / 2
        else:
            p = N
        if S < small / yed ** 2:
            if ybound == "n":
                pd = yed / xed * (
                        1 / 3 - 1 / 2 * (abs(yd - ywd) + (yd + ywd)) / yed + (yd ** 2 + ywd ** 2) / (2 * yed ** 2))
                if S > 0:
                    pd = pd + 1 / (S * yed * xed)
            else:
                pd = yed / xed * (1 / 2 * ((yd + ywd) - abs(yd - ywd)) / yed - (yd * ywd) / yed ** 2)
        else:
            u = sqrt(S + (p * pi * Ld) ** 2)
            summa = self.sumexp(u, yed)
            yd1 = yed - abs(yd - ywd)
            yd2 = yed - (yd + ywd)
            pd = 1 / 2 / u / xed * (exp(-u * abs(yd - ywd)) * summa + (
                    signum * exp(-u * (yd + ywd)) + signum * exp(-u * (yed + yd2)) + exp(-u * (yed + yd1))) * (
                                            1 + summa))
        if zbound_up == "c" and zbound_down == "c":
            part__1 = sin(p * pi * zd / zed) * sin(p * pi * zwd / zed)
        else:
            part__1 = cos(p * pi * zd / zed) * cos(p * pi * zwd / zed)
        pd_b1_pp_n = part__1 * sin(p * pi * hwd / zed / 2) * pd / p * zed / hwd * 4 / pi

        return pd_b1_pp_n

    def pd_b2_pp(self, S, xd, xwd, xed, yd, ywd, yed, xbound, ybound, Ld, zd, zwd, zed, hwd, zbound_up, zbound_down,
                 compl_type):
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
        :param Ld: безразмерная длина 
        :param zd: безразмерная координата точки, в которой рассчитывается давление
        :param zwd: безразмерная координата скважины 
        :param zed: безразмерная координата границы z 
        :param hwd: безразмерная глубина скважины 
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z 
        :param compl_type: тип скважины
        :return: 
        """
        MAXIT = 500
        summa = series_fn.nsum(
            lambda k: self.pd_b2_pp_n(k, S, xd, xwd, xed, yd, ywd, yed, xbound, ybound, Ld, zd, zwd, zed, hwd,
                                      zbound_up,
                                      zbound_down, compl_type),
            1,
            MAXIT - 1,
        )

        return summa

    def pd_b2_pp_n(self, N, S, xd, xwd, xed, yd, ywd, yed, xbound, ybound, Ld, zd, zwd, zed, hwd, zbound_up,
                   zbound_down, compl_type):
        """
        :param N: номер слагаемого
        :param S: переменная пространства Лапласа 
        :param xd: безразмерная координата точки, в которой рассчитывается давление
        :param xwd: безразмерная координата скважины 
        :param xed: безразмерная координата границы x 
        :param yd: безразмерная координата точки, в которой рассчитывается давление
        :param ywd: безразмерная координата скважины 
        :param yed: безразмерная координата границы y 
        :param xbound: тип границы 
        :param ybound: тип границы 
        :param Ld: безразмерная длина 
        :param zd: безразмерная координата точки, в которой рассчитывается давление
        :param zwd: безразмерная координата скважины 
        :param zed: безразмерная координата границы z 
        :param hwd: безразмерная глубина скважины 
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z 
        :param compl_type: тип скважины
        :return: 
        """
        MAXIT = 500
        summa = series_fn.nsum(
            lambda k: self.pd_b2_pp_n_k(N, k, S, xd, xwd, xed, yd, ywd, yed, xbound, ybound, Ld, zd, zwd, zed, hwd,
                                        zbound_up, zbound_down, compl_type),
            1,
            MAXIT - 1,
        )
        return summa

    def pd_b2_pp_n_k(self, N, k, S, xd, xwd, xed, yd, ywd, yed, xbound, ybound, Ld, zd, zwd, zed, hwd, zbound_up,
                     zbound_down, compl_type):
        """
        :param N: номер слагаемого
        :param k: номер слагаемого
        :param S: переменная пространства Лапласа 
        :param xd: безразмерная координата точки, в которой рассчитывается давление
        :param xwd: безразмерная координата скважины 
        :param xed: безразмерная координата границы x 
        :param yd: безразмерная координата точки, в которой рассчитывается давление
        :param ywd: безразмерная координата скважины 
        :param yed: безразмерная координата границы y 
        :param xbound: тип границы 
        :param ybound: тип границы 
        :param Ld: безразмерная длина 
        :param zd: безразмерная координата точки, в которой рассчитывается давление
        :param zwd: безразмерная координата скважины 
        :param zed: безразмерная координата границы z 
        :param hwd: безразмерная глубина скважины 
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z 
        :param compl_type: тип скважины
        :return: 
        """
        if xbound == "n":
            part__1 = 2 / xed * cos(k * pi * xd / xed) * cos(k * pi * xwd / xed)
        elif xbound == "c":
            part__1 = 2 / xed * sin(k * pi * xd / xed) * sin(k * pi * xwd / xed)
        else:
            part__1 = 0
        if compl_type == "frac":
            part__1 = 2 / 2 * xed / pi / k * sin(k * pi / xed) * part__1
        elif compl_type == "vert":
            part__1 = part__1
        else:
            pd_b2_pp_n_k = 0
            return pd_b2_pp_n_k
        if ybound == "n":
            signum = 1
        elif ybound == "c":
            signum = -1
        else:
            signum = 1
        if (zbound_up == "c" and zbound_down == "n") or (zbound_up == "n" and zbound_down == "c"):
            p = (2 * N - 1) / 2
        else:
            p = N
        ek = sqrt(S + k ** 2 * pi ** 2 / xed ** 2 + (p * pi * Ld) ** 2)
        smexp = self.sumexp(ek, yed)
        part__2 = exp(-ek * abs(yd - ywd)) * smexp + (
                signum * exp(-ek * (yd + ywd)) + exp(-ek * (2 * yed - abs(yd - ywd))) + signum * exp(
            -ek * (2 * yed - (yd + ywd)))) * (1 + smexp)
        part__2 = 1 / (2 * ek) * part__2
        if zbound_up == "c" and zbound_down == "c":
            part__1_1 = sin(p * pi * zd / zed) * sin(p * pi * zwd / zed)
        else:
            part__1_1 = cos(p * pi * zd / zed) * cos(p * pi * zwd / zed)
        pd_b2_pp_n_k = part__1_1 * sin(p * pi * hwd / zed / 2) * part__1 * part__2 / p * zed / hwd * 4 / pi

        return pd_b2_pp_n_k

    def pw_solve_pp(self, S, n_seg, xd, yd, xwd, ywd, zwd, xed, yed, zed, rwd, hwd, Ld, xbound, ybound, zbound_up,
                    zbound_down):
        """
        :param S: переменная пространства Лапласа 
        :param n_seg: число сегментов 
        :param xd: безразмерная координата точки, в которой рассчитывается давление
        :param yd: безразмерная координата точки, в которой рассчитывается давление
        :param xwd: безразмерная координата скважины 
        :param ywd: безразмерная координата скважины 
        :param zwd: безразмерная координата скважины 
        :param xed: безразмерная координата границы x 
        :param yed: безразмерная координата границы y 
        :param zed: безразмерная координата границы z 
        :param rwd: безразмерный радиус скважины 
        :param hwd: безразмерная глубина скважины 
        :param Ld: безразмерная длина 
        :param xbound: тип границы 
        :param ybound: тип границы 
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z 
        :return: 
        """
        matrix, rh = self.fillmatrix_pp(S, n_seg, xd, yd, xwd, ywd, zwd, xed, yed, zed, rwd, hwd, Ld, xbound, ybound,
                                        zbound_up,
                                        zbound_down)
        rh = self.LU_(matrix, rh, 2 * n_seg + 1)
        pw_solve_pp = rh[2 * n_seg + 1]

        return pw_solve_pp

    def fillmatrix_pp(self, S, n_seg, xd, yd, xwd, ywd, zwd, xed, yed, zed, rwd, hwd, Ld, xbound, ybound, zbound_up,
                      zbound_down):
        """
        :param S: переменная пространства Лапласа 
        :param n_seg: число сегментов 
        :param xd: безразмерная координата точки, в которой рассчитывается давление
        :param yd: безразмерная координата точки, в которой рассчитывается давление
        :param xwd: безразмерная координата скважины 
        :param ywd: безразмерная координата скважины 
        :param zwd: безразмерная координата скважины 
        :param xed: безразмерная координата границы x 
        :param yed: безразмерная координата границы y 
        :param zed: безразмерная координата границы z 
        :param rwd: безразмерный радиус скважины 
        :param hwd: безразмерная глубина скважины 
        :param Ld: безразмерная длина 
        :param xbound: тип границы 
        :param ybound: тип границы 
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z 
        :return: 
        """
        matrix = np.zeros(shape=(2 * n_seg + 2, 2 * n_seg + 2))
        rh = np.zeros(2 * n_seg + 2)
        dz = hwd / 2 / n_seg
        for i in range(1, 2 * n_seg + 2):
            if i <= 2 * n_seg:
                rh[i] = 0
                zd = zwd + (-hwd / 2 + dz / 2 + dz * (i - 1))
                for j in range(1, 2 * n_seg + 1):
                    zwd1 = zwd + (-hwd / 2 + dz / 2 + dz * (j - 1))
                    if S > 0:
                        if zbound_up == "n" and zbound_down == "n":
                            matrix[i, j] = self.pd_lapl_rect(S, xd, xwd, xed, yd, ywd, yed, xbound, ybound,
                                                             "vert") + self.sum_pp_1(S, xd, xwd, xed, yd, ywd, yed, zd,
                                                                                     zwd, zed, hwd, Ld, xbound, ybound,
                                                                                     zbound_up,
                                                                                     zbound_down, "vert")
                        else:
                            matrix[i, j] = self.sum_pp_1(S, xd, xwd, xed, yd, ywd, yed, zd, zwd, zed, hwd, Ld, xbound,
                                                         ybound, zbound_up, zbound_down, "vert")
                    else:
                        if zbound_up == "n" and zbound_down == "n":
                            matrix[i, j] = self.pd_rect_BD(xd, xwd, xed, yd, ywd, yed, xbound, ybound,
                                                           "vert") + self.sum_pp(0, xd, xwd, xed, yd, ywd, yed, zd, zwd,
                                                                                 zed, hwd, Ld, xbound, ybound,
                                                                                 zbound_up,
                                                                                 zbound_down, "vert")
                        else:
                            matrix[i, j] = self.sum_pp(0, xd, xwd, xed, yd, ywd, yed, zd, zwd, zed, hwd, Ld, xbound,
                                                       ybound,
                                                       zbound_up, zbound_down, "vert")
                matrix[i, 2 * n_seg + 1] = -1
            else:
                rh[i] = 1
                for j in range(2 * n_seg + 1):
                    matrix[i, j] = 1
                matrix[i, 2 * n_seg + 1] = 0

        return matrix, rh

    def pd_b3_pp(self, S, xd, xwd, xed, yd, ywd, xbound, ybound, Ld, zd, zwd, zed, hwd, zbound_up, zbound_down,
                 compl_type):
        """
        :param S: переменная пространства Лапласа 
        :param xd: безразмерная координата точки, в которой рассчитывается давление
        :param xwd: безразмерная координата скважины 
        :param xed: безразмерная координата границы x 
        :param yd: безразмерная координата точки, в которой рассчитывается давление
        :param ywd: безразмерная координата скважины 
        :param xbound: тип границы 
        :param ybound: тип границы 
        :param Ld: безразмерная длина 
        :param zd: безразмерная координата точки, в которой рассчитывается давление
        :param zwd: безразмерная координата скважины 
        :param zed: безразмерная координата границы z 
        :param hwd: безразмерная глубина скважины 
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z 
        :param compl_type: тип скважины
        :return: 
        """
        summa = 0
        k = 1
        MAXIT = 500
        while True:
            psum = 0
            psumabs = 0
            for i in range(1, self.part_sum_num):
                add = self.pd_b3_pp_n(k, S, xd, xwd, xed, yd, ywd, xbound, ybound, Ld, zd, zwd, zed, hwd, zbound_up,
                                      zbound_down, compl_type)
                psum = psum + add
                psumabs = psumabs + abs(add)
                k = k + 1
            summa += psum
            if not (abs(psumabs) / self.part_sum_num >= self.tiny * (abs(summa) + self.tiny) and k < MAXIT):
                break
        pd_b3_pp = summa

        return pd_b3_pp

    def pd_b3_pp_n(self, N, S, xd, xwd, xed, yd, ywd, xbound, ybound, Ld, zd, zwd, zed, hwd, zbound_up, zbound_down,
                   compl_type):
        """
        :param N: номер слагаемого
        :param S: переменная пространства Лапласа 
        :param xd: безразмерная координата точки, в которой рассчитывается давление
        :param xwd: безразмерная координата скважины 
        :param xed: безразмерная координата границы x 
        :param yd: безразмерная координата точки, в которой рассчитывается давление
        :param ywd: безразмерная координата скважины 
        :param xbound: тип границы 
        :param ybound: тип границы 
        :param Ld: безразмерная длина 
        :param zd: безразмерная координата точки, в которой рассчитывается давление
        :param zwd: безразмерная координата скважины 
        :param zed: безразмерная координата границы z 
        :param hwd: безразмерная глубина скважины 
        :param zbound_up: верхняя граница по z
        :param zbound_down: нижняя граница по z 
        :param compl_type: тип скважины
        :return: 
        """
        signum = 1
        if xbound == "c":
            signum = -1
        if (zbound_up == "c" and zbound_down == "n") or (zbound_up == "n" and zbound_down == "c"):
            p = (2 * N - 1) / 2
        else:
            p = N
        u = S + (p * pi * Ld) ** 2
        if compl_type == 'vert':
            xd1m00 = (xd - xwd) ** 2
            xd1p00 = (xd + xwd) ** 2
            yd1 = (yd - ywd) ** 2
            sum_ = (self.unit_cylinder_source(u, sqrt(xd1m00 + yd1)) +
                    signum * self.unit_cylinder_source(u, sqrt(xd1p00 + yd1)))
            sum_ += series_fn.nsum(
                lambda k: (
                        signum
                        * self.unit_cylinder_source(
                    u, sqrt((xd + xwd + 2 * xed * k) ** 2 + yd1)
                )
                        + self.unit_cylinder_source(
                    u, sqrt((xd - xwd + 2 * xed * k) ** 2 + yd1)
                )
                        + signum
                        * self.unit_cylinder_source(
                    u, sqrt((xd + xwd - 2 * xed * k) ** 2 + yd1)
                )
                        + self.unit_cylinder_source(
                    u, sqrt((xd - xwd - 2 * xed * k) ** 2 + yd1)
                )
                ),
                1,
                np.inf,
            )
        elif compl_type == "frac":
            yd1 = abs(yd - ywd)
            sum_ = (self.unit_fracture_func(u, abs(xd - xwd), yd1) +
                    signum * self.unit_fracture_func(u, abs(xd + xwd), yd1))
            sum_ += series_fn.nsum(
                lambda k: (
                        signum
                        * self.unit_fracture_func(
                    u, abs(xd + xwd + 2 * k * xed), yd1)
                        + self.unit_fracture_func(
                    u, abs(xd - xwd + 2 * k * xed), yd1)
                        + signum * self.unit_fracture_func(
                    u, abs(xd + xwd - 2 * k * xed), yd1)
                        + self.unit_fracture_func(
                    u, abs(xd - xwd - 2 * k * xed), yd1)
                ),
                1,
                np.inf,
            )
        else:
            sum_ = 0
        if zbound_up == "c" and zbound_down == "c":
            part__1 = sin(p * pi * zd / zed) * sin(p * pi * zwd / zed)
        else:
            part__1 = cos(p * pi * zd / zed) * cos(p * pi * zwd / zed)

        pd_b3_pp_n = part__1 * sin(p * pi * hwd / zed / 2) * sum_ / p * zed / hwd * 4 / pi

        return pd_b3_pp_n
