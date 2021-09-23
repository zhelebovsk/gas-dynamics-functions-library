import GDF
import numpy as np
from matplotlib import pyplot as plt


def calc(d, Asupsonic, sonic, k):
    A = np.pi * d * d / 4
    q = Asupsonic / A
    if sonic:
        la = GDF.qrsup(q, k)
    else:
        la = GDF.qrsub(q, k)
    tau = GDF.tau(la, k)
    pi = GDF.pi(la, k)
    p = p0 * pi
    t = t0 * tau
    return tau, pi, la, p, t


if __name__ == "__main__":
    G = 12.5
    p0 = 23*100000
    ph = 0.8*100000
    psk = 1 * 100000
    t0 = 1120
    alpha = np.deg2rad(26)
    beta = np.deg2rad(11)
    k = 1.4
    R = 287.1
    r = 2
    m = GDF.mk(k, R)

    la = np.zeros([6])
    q = np.zeros([6])
    eps = np.zeros([6])
    tau = np.zeros([6])
    pi = np.zeros([6])
    t = np.zeros([6])
    p = np.zeros([6])
    d = np.zeros([6])
    A = np.zeros([6])
    v = np.zeros([6])
    ro = np.zeros([6])
    x = np.zeros([6])
    la[0] = 0
    tau[0] = GDF.tau(la[0], k)
    pi[0] = GDF.pi(la[0], k)
    q[0] = GDF.q(la[0], k)
    t[0] = tau[0] * t0
    p[0] = pi[0] * p0

    # Рассмотрим критическое сечение (lambda = 1)

    la[3] = 1
    q[3] = GDF.q(la[3], k)
    A[3] = G * np.sqrt(t[0]) / m / p[0] / q[3]
    d[3] = np.sqrt(4*A[3]/np.pi)
    pi[3] = GDF.pi(la[3], k)
    tau[3] = GDF.tau(la[3], k)
    p[3] = p[0] * pi[3]
    t[3] = t[0] * tau[3]
    # Входное сечение
    x[1] = 0
    l1 = r * d[3]
    x[3] = l1
    x[2] = l1/2
    d[1] = d[3] + 2 * l1 * np.tan(alpha/2)
    A[1] = np.pi * d[1]*d[1] / 4
    q[1] = A[3]/A[1]# ???
    la[1] = GDF.qrsub(q[1], k)
    tau[1] = GDF.tau(la[1], k)
    pi[1] = GDF.pi(la[1], k)
    p[1] = p[0] * pi[1]
    t[1] = t[0] * tau[1]
    # Промежуточная точка
    #x2 = l1/2
    d[2] = (d[1]+d[3]) / 2
    A[2] = np.pi * d[2]*d[2] / 4
    q[2] = A[3] / A[2]
    la[2] = GDF.qrsub(q[2], k)
    tau[2] = GDF.tau(la[2], k)
    pi[2] = GDF.pi(la[2], k)
    p[2] = p0 * pi[2]
    t[2] = t0 * tau[2]
    # Выходное сечение
    p[5] = ph
    pi[5] = p[5]/p[0]
    la[5] = GDF.pir(pi[5], k)
    tau[5] = GDF.tau(la[5], k)
    q[5] = GDF.q(la[5], k)
    p[5] = p[0] * pi[5]
    t[5] = t[0] * tau[5]
    A[5] = A[3]/q[5]
    d[5] = np.sqrt(4*A[5]/np.pi)
    l2 = (d[5] - d[3]) / (2*np.tan(beta/2))
    x[5] = l1 + l2
    x[4] = l1 + l2/2
    # Промежуточное сечение
    d[4] = (d[3] + d[5]) / 2
    A[4] = np.pi * d[4]*d[4] / 4
    q[4] = A[3] / A[4]
    la[4] = GDF.qrsup(q[4], k)
    tau[4] = GDF.tau(la[4], k)
    pi[4] = GDF.pi(la[4], k)
    p[4] = p[0] * pi[4]
    t[4] = t[0] * tau[4]

    eps[0] = GDF.eps(la[0], k)
    ro[0] = p[0] / R / t[0]
    eps[3] = GDF.eps(la[3], k)
    ro[3] = ro[0] * eps[3]
    eps[1] = GDF.eps(la[1], k)
    ro[1] = ro[0] * eps[1]
    eps[2] = GDF.eps(la[2], k)
    ro[2] = ro[0] * eps[2]
    eps[5] = GDF.eps(la[5], k)
    ro[5] = ro[0] * eps[5]
    eps[4] = GDF.eps(la[4], k)
    ro[4] = ro[0] * eps[4]
    #Определить положения x, построить графики

    #plt.plot(x[1:], la[1:])
    #plt.plot(x[1:], eps[1:])
    #plt.plot(x[1:], pi[1:])
    #plt.plot(x[1:], tau[1:])
    #plt.plot(x[1:], d[1:]/2 * 10)


    la = np.zeros([1000])
    pi = np.zeros([1000])
    tau = np.zeros([1000])
    p = np.zeros([1000])
    t = np.zeros([1000])

    d1 = np.linspace(d[1], d[3], 500)
    d2 = np.linspace(d[3], d[5], 500)
    x1 = np.linspace(x[1], x[3], 500)
    x2 = np.linspace(x[3], x[5], 500)
    x = np.hstack((x1,x2))
    d = np.hstack((d1,d2))
    for i, point in enumerate(d):
        if i < np.size(la)/2:
            tau[i], pi[i], la[i], p[i], t[i] = calc(point, A[3], False, k)
        else:
            tau[i], pi[i], la[i], p[i], t[i] = calc(point, A[3], True, k)
    #plt.plot(x[1:], la[1:])
    #plt.plot(x[1:], pi[1:])
    plt.plot(x[1:], p[1:])
    plt.show()