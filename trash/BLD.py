import math
import GDF


def bld(k, R, b2, d2, c2r, alf2, t2s, t1, p1 ,eff, p2s,c2):
    if b2/d2 > 0.06:
        h3h2 = 0.8
    elif b2/d2 > 0.04:
        h3h2 = 1.0
    else:
        h3h2 = 1.2
    b3 = b2*h3h2
    # диаметр изменения канала 1..1.05
    d2c = d2 * 1.00
    # d3/d2 = 1.1..1.5 протяженность блд
    d3 = d2 * 1.2
    f3r = 3.1415 * d3 * b3
    c3r = c2r * (d2*b2)/(d3*b3)
    alf3 = math.atan(math.tan(alf2)/(b3/b2*1.02))
    c3u = c3r/math.tan(alf3)
    c3 = math.sqrt(c3r * c3r + c3u * c3u)
    t3s = t2s
    t3 = t3s - c3*c3/(2*k*R/(k-1))
    lac3 = c3/math.sqrt(2*k/(k+1)*R*t3s)
    p3 = p1*math.pow(t3/t1 , k/(k-1)*eff)
    p3s = p3/GDF.pi(lac3,k)
    sigmabld = p3s/p2s
    ro3 = p3/(R*t3)
    nuebld = 2*math.sqrt(b3/d3)/(1+ math.sqrt(d3/d2))*math.sin(alf3)
    nuebld = math.atan(nuebld)*2
    nuebld = math.degrees(nuebld)
    ebld = 0.147 + 0.0046* (nuebld - 12) * (nuebld - 12)
    dh = ebld * c2*c2/2
