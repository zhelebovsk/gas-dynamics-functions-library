import math
# упрощенный вариант определения давления и температуры (не по ГОСТ)
# до 20 км
# возможно изменение температуры и давления на заданную величину
# изменение температуры и давления раздельно (например, МСА + 10 град -- 101325 Па сохраняется)


def t(h,dt):
    if h < 11000.0:
        x = 288.15 - 6.5 * h / 1000.0 + dt
    else:
        x = 288.15 - 6.5 * 11000.0 / 1000.0 + dt
    return x


def p(h,dp):
    if h < 11000.0:
        x = math.pow(1-0.0065*h/288.15, 5.2561)
        x = 101325.0 * x + dp
    else:
        x = math.pow(1 - 0.0065 * 11000.0 / 288.15, 5.2561)
        x = 101325.0 * x * math.exp(-9.81/(287.0*t(11000.0))*(h - 11000.0)) + dp
    return x