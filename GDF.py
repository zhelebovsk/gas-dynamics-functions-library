import math
import numpy as np

# перед использованием требуется проверить точность (числа Маха)
# поиск обратных функций выполняется с заданой точностью


"""
1. переделать функции обратные температуры, давления, плотности на прямые (без обращения к базовым) (DONE)
2. взаимосвязь лямбда - мах (+ обратная)
3. проверить правильность функций на таблицах (к = 1.4) (циклом)
4. добавить описание к каждой функции через 
5. начальное приближение для обратной газодинамической функции расхода (DONE)
"""


def q(la, k):
    x = la*math.pow(((k+1.0)/2.0)*(1.0-(k-1.0)/(k+1.0)*la*la), (1.0/(k-1.0)))
    return x


def pi(la, k):
    x = 1.0 - (k - 1.0)/(k + 1.0)*la*la
    x = math.pow(x,k/(k-1))
    return x


def tau(la, k):
    x = 1.0-(k-1.0)/(k+1.0)*la*la
    return x


def eps(la, k):
    x = 1.0-(k-1.0)/(k+1.0)*la*la
    x = math.pow(x, 1.0/(k-1.0))
    return x


def y(la,k):
    x = pi(la, k) * q(la, k)
    return x


def z(la):
    x = (la + 1.0/la)/2.0
    return x


def M(la, k):
    x = (2.0/(k+1.0))/tau(la, k)
    x = la*math.sqrt(x)
    return x


def qrsub(qq, k):
    la = qq
    num = 0.0000001
    while math.fabs((q(la, k) - qq) / qq) > num:
        dydx = (q(la, k) - q(la + num, k)) / num
        la = la + (q(la, k) - qq) / dydx
    if la > 1.0:
        print('la(q) > 1')
    return la


def qrsup(qq, k):
    la = (1 - qq) * (np.power((k+1)/(k-1), 0.5) - 1) + 1
    num = 0.0000001
    while math.fabs((q(la, k) - qq) / qq) > num:
        dydx = (q(la, k) - q(la + num, k)) / num
        la = la + (q(la, k) - qq) / dydx
    if la < 1.0:
        print('la(q) < 1')
    return la


def pir(pipi, k):
    la = (k + 1) / (k - 1) * (1 - np.power(pipi, (k-1)/k))
    return np.sqrt(la)


def taur(tautau, k):
    la = (k+1) / (k-1) * (1 - tautau)
    return np.sqrt(la)


def epsr(epseps, k):
    la = (k + 1) / (k - 1) * (1 - np.power(epseps, (k-1)))
    return np.sqrt(la)


def mk(k, R):
    x = math.sqrt(k/R*math.pow(2.0/(k+1.0),(k+1.0)/(k-1.0)))
    return x

breakpoint()
