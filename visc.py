import numpy as np


def ro(p, T, R):
    return p/(R*T)


def sutherland(T, b, S):
    return b * np.power(T, 3/2)/(T + S)


def airvisc(T):
    b = 1.458 * np.power(10.0, -6)
    S = 110.4
    return sutherland(T, b, S)


def airdens(p, T):
    R = 287
    return ro(p, T, R)


print(airdens(101325, 288))
print(airvisc(288 + 1000))
