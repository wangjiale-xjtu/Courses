import numpy as np
import matplotlib.pyplot as plt
# 生成A矩阵
def set_a(scale):
    A = np.zeros((scale, scale))

    A[0, 0] = -2
    A[0, 1] = 1

    A[scale-1, scale-1] = -2
    A[scale-1, scale-2] = 1
    for i in range(1, scale-1):
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


scale = 400
A = set_a(scale)
b = set_b(scale)

# 给定初始向量及精度
Xplot=[]
Yplot=[]
xk = np.zeros(scale)
accuracy = 1e-6
rk = b-A@xk
dk = rk
ite = 0
while(1):
    ite += 1

    alph = rk.T@rk/(dk.T@A@dk)
    xk = xk + alph*dk
    rkk = b - A@xk
    Xplot.append(ite)
    Yplot.append(np.log(np.linalg.norm(rkk)))
    if np.linalg.norm(rkk) <= accuracy:
        break
    else:
        beta = (np.linalg.norm(rkk)*np.linalg.norm(rkk))/(np.linalg.norm(rk)*np.linalg.norm(rk))
        dk = rkk + beta*dk
        rk = rkk
print("求解完成,迭代次数"+str(ite))
print(xk)
print(Yplot)

plt.rcParams['font.sans-serif']=['SimHei']
plt.rcParams['axes.unicode_minus'] = False

plt.plot(Xplot, Yplot)
plt.xlabel("迭代次数")#x轴上的名字
plt.ylabel("共轭梯度法误差")#y轴上的名字
plt.grid()
plt.show()

