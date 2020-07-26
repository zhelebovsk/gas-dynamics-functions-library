import math
import GDF
def compr(k,R,tinl,pik,eff,Hz,pinl,d_1k,d_1v,n,G,beta2):
    Lks = (k/(k-1.0))*R*tinl*(math.pow(pik,(k-1.0)/k)-1.0)
    Lz = Lks/eff
    pk = pinl*pik
    tk = tinl + Lz/(k*R/(k-1.0))
    u2 = math.sqrt(Lz/Hz)
    #print(u2,tk,Lz,Lks)
    d2 = 60.0*u2/(3.1415*n)
    d1k = d_1k*d2
    d1v = d_1v *d2
    d1sr = (d1k+d1v)/2.0
    h1 = (d1k-d1v)/2.0
    #ширина рабочего колеса s = 0.15..0.35
    S = 0.22*d2
    u1sr = 3.1415*d1sr*n/60.0
    #условный коэффициент расхода ступени (0.05..0.12)
    ro = pinl/(R*tinl)
    fip = 4.0*G/(3.14*ro*d2*d2*u2)
    print('fip = ',fip)
    print('check?')
    #zrk = beta2/4.0 + (105.0-beta2)*(beta2-10)/(200.0)
    #print(zrk)
    zrk = 16.0
    F1a = 3.1415/4.0*(d1k*d1k-d1v*d1v)
    #F1a = F1a - zrk*0.5*b1*(tk+tvt)/math.sin(beta1)
    mk = GDF.mk(k,R)
    q = G*math.sqrt(tinl)/(mk*pinl*F1a)

    print(GDF.qr(q, k, 0.2))
    la0 = GDF.qr(q, k, 0.2)

    c1a = la0*math.sqrt(2.0*k/(k+1.0)*R*tinl)
    print(c1a/u2)
    print('c1a = ',c1a)
    #c1u = c_1u*u1sr
    c1u = 0.0
    c1 = math.sqrt(c1a*c1a+c1u*c1u)
    eps = 0.000001
    if math.fabs(c1u) < eps:
        alf1 = 3.1415/2
    else:
        alf1 = math.atan(c1a/c1u)
    lac1 = c1/math.sqrt(2.0*k/(k+1.0)*R*tinl)
    p1 = pinl*GDF.pi(lac1,k)
    t1 = tinl * GDF.tau(lac1,k)
    ro1 = p1/(R*t1)
    #если закрутка в сторону вращения
    w1u = u1sr-c1u
    w1 = math.sqrt(c1a * c1a + w1u * w1u)

    if math.fabs(c1u) < eps:
        beta1 = 3.1415/2.0
    else:
        beta1 = math.atan(c1a/w1u)

    t1w = t1 + w1*w1/((2.0*k/(k+1.0))*R)

    law1 = w1/math.sqrt(2.0*k/(k+1.0)*R*t1w)
    print(law1)

# параметры на выхо де из рабочего колеса
    betatr = 1.05
    c2u = (Lz/(1+betatr) - c1u * u1sr) / u2
    #формула Виснера
    mu = 1 - math.sqrt(math.sin(beta2))/math.pow(zrk, 0.7)
    c2uinf = c2u/mu
    w2uinf = u2 - c2uinf
    c2r = w2uinf * math.tan(beta2)
    w2r = c2r
    alf2 = math.atan(c2r/c2u)
    w2u = u2 - c2u
    beta2 = math.atan(c2r/w2u)
    w2 = math.sqrt(w2u*w2u + w2r*w2r)
    c2 = math.sqrt(c2u * c2u + c2r * c2r)
    t2 = t1 + Lz/(k/(k-1)*R) + (c1*c1 - c2*c2)/(2*k/(k-1)*R)
    p2 = p1 * math.pow(t2/t1, k/(k-1)*eff)
    ro2 = p2/(R*t2)
    t2s = t2 + c2*c2/(2*k*R/(k-1))
    tw2s = t2 + w2 * w2 / (2 * k * R / (k - 1))
    lac2 = c2/math.sqrt(2*k/(k+1)*R*t2s)
    law2 = w2 / math.sqrt(2 * k / (k + 1) * R * tw2s)
    p2s = p2 / GDF.pi(lac2, k)
    pw2s = p2 / GDF.pi(law2, k)
    picheck = p2s/pinl
    # площадь на выходе из колеса
    F2a = G/(c2r*ro2)
    h2 = F2a / (3.1415*d2)

