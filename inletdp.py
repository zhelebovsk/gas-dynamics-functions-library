import math
import GDF
import numpy as np
import matplotlib.pyplot as plt

"""
1. добавить комментарии по коду
2. пробелы перед запятой и прочая красота
"""

# Задаются параметры воздуха на входе (торможения), расход, требуемые скорости (по приведенной)
# Давление по каналу изменяется линейно (по длине)

G = 10.0
k = 1.4
Cp = 1005.0
T0s = 288.15
P0s = 101325.0
R = 287.0
mk = GDF.mk(k,R)

# длина канала
l = 3.0
# скорость на вхлде и выходе из канала
la0 =0.05
lak = 0.4

# площадь на входе и на выходе из канала
# по уравнению расхода S*q(la) = const
S0 = G*math.sqrt(R*T0s)/(mk*GDF.q(la0,k)*P0s)
Sk = S0*GDF.q(la0,k)/GDF.q(lak,k)
print('S0 = ',S0)
print('Sk = ',Sk)
print('-----------')
# рейт падения давления по длине
dpidx = (GDF.pi(lak,k)-GDF.pi(la0,k))/l
print('dpi(la)/dx = ', dpidx)
print('-----------')
#Давления и температура на входе и выходе из канала
T0 = T0s * GDF.tau(la0,k)
Tk = T0s * GDF.tau(lak,k)
print('T0s = ', T0s)
print('T0 = ', T0)
print('Tk = ', Tk)
print('-----------')
P0 = P0s * GDF.pi(la0,k)
Pk = P0s * GDF.pi(lak,k)
print('P0s = ', P0s)
print('P0 = ', P0)
print('Pk = ', Pk)


print('-----------')
# определение параметров потока в n сечений
n = 100
la = np.zeros((n+1))
q = np.zeros((n+1))
S = np.zeros((n+1))
d = np.zeros((n+1))
# заполнение массивов данными
x = np.arange(0, l+l/n, l/n)
pi = GDF.pi(la0, k) + dpidx * x
#la = GDF.pir(pi, k, 0.2)
for i in range(n+1):
    la[i] = GDF.pir(pi[i], k, lak)
    q[i] = GDF.q(la[i], k)
    S[i] = G*math.sqrt(R*T0s)/(mk*q[i]*P0s)
    d[i] = math.sqrt(4*S[i]/3.1415)
#print(x)
#print(pi)
print(la)
#print(q)
#print(S)
#print(d)
#plt.plot(x, pi)
plt.plot(x, d/2)
plt.plot(x, -d/2)
plt.show()
# оценить погрешность расчета относительно начального