import math

x0 = 1.0
r0 = 0.3

R = 1/(2*r0)*(x0*x0+r0*r0)
print(R)
x = 0.0
ycir = math.sqrt(R*R - (x-x0)*(x-x0)) + r0 - R

n = 1 # 1..0.3 ньютоновские тела?
yn = math.pow(x/x0,n)

k= 1 # парабола
yp = (2*(x/x0) - k*(x/x0))/(2-k)


