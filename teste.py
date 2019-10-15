import numpy as np 
import matplotlib.pyplot as plt 
from statsmodels.graphics.tsaplots import plot_acf

n = 3000
x = np.zeros(n)
v = np.zeros(n)
dt = 0.00001

for i in range(1,n):
    q = np.random.rand()*2+1
    if np.random.rand() > 20:
        q = 0
        h = 0
    else:
        h = 2

    b = 0.1
    s = 0.0001
    v[i] = dt * ( -v[i-1]*h + np.sqrt(1/b*dt)*q)/s + v[i-1]
    x[i] = x[i] + v[i]*dt 

plt.figure()
plt.plot(x)
plt.title("espaço")

plt.figure()
plt.plot(v)
plt.title("velocidade")

plot_acf(v[2000:],lags=100)

plt.show()