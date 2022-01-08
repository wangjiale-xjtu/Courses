import numpy as np
import matplotlib.pyplot as plt

# 生成A矩阵
def set_a(scale):
    A = np.zeros((scale, scale))

    A[0, 0] = -2
    A[0, 1] = 1

    A[scale - 1, scale - 1] = -2
    A[scale - 1, scale - 2] = 1
    for i in range(1, scale - 1):
        A[i, i] = -2
        A[i, i - 1] = 1
        A[i, i + 1] = 1
    return A


# 生成b列向量
def set_b(scale):
    b = np.zeros(scale)
    b[0] = -1
    b[-1] = -1
    return b

kmax=10000000000
scale = 400
A = set_a(scale)
b = np.array([set_b(scale)]).T
Xplot=[]
Yplot=[]
x0 = np.zeros((scale, 1))
accuracy = 1e-6
rk = b - A @ x0

a = (rk.T @ rk) / ((A @ rk).T @ rk)
xk = x0 + a * rk
k=1
while np.linalg.norm(xk-x0) >= accuracy:
    x0 = xk
    rk = b - A @ x0
    a = (rk.T @ rk) / ((A @ rk).T @ rk)
    xk = x0 + a * rk
    Xplot.append(k)
    Yplot.append(np.log(np.linalg.norm(rk)))
    k=k+1
    if k >= kmax:
        print("超出最大迭代次数")
        break
print('迭代次数'+str(k))
print(xk)

plt.rcParams['font.sans-serif']=['SimHei']
plt.rcParams['axes.unicode_minus'] = False

plt.plot(Xplot, Yplot)
plt.xlabel("迭代次数")#x轴上的名字
plt.ylabel("最速下降法误差")#y轴上的名字
plt.grid()
plt.show()

