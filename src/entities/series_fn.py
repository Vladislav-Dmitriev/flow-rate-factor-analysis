import math


def nsum(fn, a, b, precision=7):
    """Суммирования функции от a по b,
    работает с монотонно убывающими функциями

    :param fn: функция
    :type fn: Callable[[int|float], int|float]
    :param a: левая граница суммирования
    :type a: int|float
    :param b: правая граница суммирования (может быть np.inf)
    :type b: int|float
    :param precision: точность (до скольки знаков после запятой), по умолчанию 1
    :type precision: Optional[int]
    :return: сумма
    :rtype: int|float
    """
    eps = 10 ** - precision
    inc = fn(a)
    sum_ = inc
    while abs(inc) > eps and a <= b:
        a += 1
        inc = fn(a)
        sum_ = math.fsum((sum_, inc))
    return sum_
