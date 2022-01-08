import numpy as np
import math
from sympy import *


def solve_NR(f, x, X, N_var):
    # define symbol variable


    # 求偏导
    dic={}
    grad=np.zeros((N_var,N_var))
    x_value=np.zeros(N_var)
    current_f = np.zeros((N_var,1))
    for i in range(N_var):
        x_value[i] = X[i, 0]
        dic[x[i]] = x_value[i]
    for i in range(N_var):
        current_f[i, 0] = float(f[i].subs(dic).evalf())
    for i in range(N_var):
        for j in range(N_var):
            grad[i, j] = diff(f[i], x[j]).subs(dic).evalf()

    iter_cnt = 0
    epsilon = 10 ** (-8)

    # print('current_F \n',F)
    while np.linalg.norm(current_f) >= epsilon:
        iter_cnt += 1
        print('迭代次数'+str(iter_cnt))
        print('误差'+str(np.linalg.norm(current_f)))
        # print('iter_cnt',#iter_cnt)
        # 矩阵的转置
        X = X - np.dot(np.linalg.inv(grad), current_f)
        for i in range(N_var):
            x_value[i] = X[i, 0]
            dic[x[i]] = x_value[i]
        for i in range(N_var):
            current_f[i, 0] = float(f[i].subs(dic).evalf())
        for i in range(N_var):
            for j in range(N_var):
                grad[i, j] = diff(f[i], x[j]).subs(dic).evalf()
    print('误差' + str(np.linalg.norm(current_f)))
    return X


N_var = 2
X = np.array([[1.04], [0.47]])
x=[]
for i in range(N_var):
    x.append(symbols('x'+str(i+1)))
#f=[x[0] ** 2 + x[1] ** 2 - x[2] ** 2 - 1,2 * x[0] ** 2 + x[1] ** 2 - 4 * x[2],3 * x[0] ** 2 - 4 * x[1] ** 2 + x[2] ** 2]
f=[cos(x[0]**2+0.4*x[1])+x[0]**2+x[1]**2-1.6,1.5*x[0]**2-1/0.36*x[1]**2-1]
X = solve_NR(f,x,X, N_var)
print(X)
