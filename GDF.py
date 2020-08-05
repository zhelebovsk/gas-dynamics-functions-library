import math
#перед использованием требуется проверить точность (числа Маха)
#поиск обратных функций выполняется с заданой точностью


"""
1. переделать функции обратные температуры, давления, плотности на прямые (без обращения к базовым)
2. взаимосвязь лямбда - мах (+ обратная)
3. проверить правильность функций на таблицах (к = 1.4)
4. добавить описание к каждой функции через breakpoint() (или циклами)
5. начальное приближение для обратной газодинамической функции расхода
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
    x = math.pow(x,1.0/(k-1.0))
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


def qr(qq, k, la):
    num = 0.0000001
    while math.fabs(q(la, k) - qq) > num:
        dydx = (q(la, k) - q(la + num, k)) / num
        la = la + (q(la, k) - qq) / dydx
    if la > 1.0:
        print('la(q) > 1')
    return la


def pir(pipi, k, la):
    num = 0.0000001
    while math.fabs(pi(la, k) - pipi) > num:
        dydx = (pi(la, k) - pi(la + num, k)) / num
        la = la + (pi(la, k) - pipi) / dydx
    if la > 1.0:
        print('la(pi) > 1')
    return la


def taur(tautau, k, la):
    num = 0.0000001
    while math.fabs(tau(la, k) - tautau) > num:
        dydx = (tau(la, k) - tau(la + num, k)) / num
        la = la + (tau(la, k) - tautau) / dydx
    if la > 1.0:
        print('la(tau) > 1')
    return la


def epsr(epseps, k, la):
    num = 0.0000001
    while math.fabs(eps(la, k) - epseps) > num:
        dydx = (eps(la, k) - eps(la + num, k)) / num
        la = la + (eps(la, k) - epseps) / dydx
    if la > 1.0:
        print('la(eps) > 1')
    return la


def mk(k, R):
    x = math.sqrt(k/R*math.pow(2.0/(k+1.0),(k+1.0)/(k-1.0)))
    return x

