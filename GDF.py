import math
import numpy as np

# поиск обратных функций выполняется с заданой точностью\

"""
----1. переделать функции обратные температуры, давления, плотности на прямые (без обращения к базовым)
----2. взаимосвязь лямбда - мах ( + обратная)
3. проверить правильность функций на таблицах (к = 1.4) (циклом)
----4. добавить описание к каждой функции через 
----5. начальное приближение для обратной газодинамической функции расхода
"""


def q(la, k):
    """
    ГДФ расхода по приведенной скорости и показателю адиабаты
    :param la: приведенная скорость
    :param k: показатель адиабаты
    :return: ГДФ расхода
    """
    x = la * math.pow(((k + 1.0) / 2.0) * (1.0 - (k - 1.0) / (k + 1.0) * la * la), (1.0 / (k - 1.0)))
    return x


def pi(la, k):
    """
    ГДФ давления по приведенной скорости и показателю адиабаты
    :param la: приведенная скорость
    :param k: показатель адиабаты
    :return: ГДФ давления
    """
    x = 1.0 - (k - 1.0) / (k + 1.0) * la * la
    x = math.pow(x, k / (k - 1.0))
    return x


def tau(la, k):
    """
    ГДФ температуры по приведенной скорости и показателю адиабаты
    :param la: приведенная скорость
    :param k: показатель адиабаты
    :return: ГДФ температуры
    """
    x = 1.0 - (k - 1.0) / (k + 1.0) * la * la
    return x


def eps(la, k):
    """
    ГДФ плотности по приведенной скорости и показателю адиабаты
    :param la: приведенная скорость
    :param k: показатель адиабаты
    :return: ГДФ плотности
    """
    x = 1.0 - (k - 1.0) / (k + 1.0) * la * la
    x = math.pow(x, 1.0 / (k - 1.0))
    return x


def y(la, k):
    """
    Произведение ГДФ q, pi по приведенной скорости и показателю адиабаты для уравнения расхода
    (используется в уравнении расхода при статических параметрах)
    :param la: приведенная скорость
    :param k: показатель адиабаты
    :return: ГДФ y
    """
    x = pi(la, k) * q(la, k)
    return x


def z(la):
    """
    ???
    :param la: приведенная скорость
    :return: ???
    """
    x = (la + 1.0/la) / 2.0
    return x


def M(la, k):
    """
    Число Маха по приведенной скорости и показателю адиабаты
    :param la: приведенная скорость
    :param k: показатель адиабаты
    :return: Число Маха
    """
    x = (2.0 / (k + 1.0)) / tau(la, k)
    x = la * math.sqrt(x)
    return x


def laM(Mach,k):
    """
    Приведенная скорость по числу Маха и показателю адиабаты
    :param Mach: число Маха
    :param k: показатель адиабаты
    :return: приведенная скорость
    """
    x = (k + 1.0) / 2.0 * Mach * Mach
    x = x / (1.0 + (k - 1.0) / 2.0 * Mach * Mach)
    return np.sqrt(x)


def qrsub(qq, k):
    """
    Приведенная скорость по значению ГДФ расхода (дозвук) и показателю адиабаты
    :param qq: значение ГДФ расхода
    :param k: показатель адиабаты
    :return: приведенная скорость
    """
    la = qq
    num = 0.0000001
    while math.fabs((q(la, k) - qq) / qq) > num:
        dydx = (q(la, k) - q(la + num, k)) / num
        la = la + (q(la, k) - qq) / dydx
    if la > 1.0:
        print('la(q) > 1')
    return la


def qrsup(qq, k):
    """
    Приведенная скорость по значению ГДФ расхода (сверхзвук) и показателю адиабаты
    :param qq: значение ГДФ расхода
    :param k: показатель адиабаты
    :return: приведенная скорость
    """
    la = (1.0 - qq) * (np.power((k + 1.0) / (k - 1.0), 0.5) - 1.0) + 1.0
    num = 0.0000001
    while math.fabs((q(la, k) - qq) / qq) > num:
        dydx = (q(la, k) - q(la + num, k)) / num
        la = la + (q(la, k) - qq) / dydx
    if la < 1.0:
        print('la(q) < 1')
    return la


def pir(pipi, k):
    """
    Приведенная скорость по значению ГДФ давления и показателю адиабаты
    :param pipi: значение ГДФ давления
    :param k: показатель адиабаты
    :return: приведенная скорость
    """
    la = (k + 1.0) / (k - 1.0) * (1.0 - np.power(pipi, (k - 1.0) / k))
    return np.sqrt(la)


def taur(tautau, k):
    """
    Приведенная скорость по значению ГДФ температуры и показателю адиабаты
    :param tautau: значение ГДФ температуры
    :param k: показатель адиабаты
    :return: приведенная скорость
    """
    la = (k + 1.0) / (k - 1.0) * (1.0 - tautau)
    return np.sqrt(la)


def epsr(epseps, k):
    """
    Приведенная скорость по значению ГДФ плотности и показателю адиабаты
    :param epseps: значение ГДФ плотности
    :param k: показатель адиабаты
    :return: приведенная скорость
    """
    la = (k + 1.0) / (k - 1.0) * (1.0 - np.power(epseps, (k - 1.0)))
    return np.sqrt(la)


def mk(k, r):
    """
    Коэффициент в уравнении расхода по газовой постоянной и показателю адиабаты
    :param r: газовая постоянная, Дж/(кг*К)
    :param k: показатель адиабаты
    :return: коэффициент в уравнении расхода, (Дж/(кг*К))^(-0.5)
    """
    x = math.sqrt(k / r * math.pow(2.0 / (k + 1.0), (k + 1.0) / (k - 1.0)))
    return x

breakpoint()
