import math
from matplotlib import pyplot as plt
import numpy as np
#
x0 = 1.0
r0 = 0.3
x = np.linspace(0, x0, 100)
R = 1/(2*r0)*(x0*x0+r0*r0)
print(R)
ycir = np.sqrt(R*R - (x-x0)*(x-x0)) + r0 - R

n = 0.5 # 1..0.3 ньютоновские тела?
yn = np.power(x/x0, n) * r0

k = 1.0 # парабола 0..1
yp = ((2*(x/x0) - k*(x/x0)*(x/x0))/(2-k)) * r0

plt.plot(x,yn)
plt.plot(x,ycir)
plt.plot(x,yp)
plt.show()