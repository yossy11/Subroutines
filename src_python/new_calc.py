import math
import csv
# import pdb
import numpy as np


def calc_S(c_params, stress):
    """
        Calculate S1,S2,S3 used for yld calculation

        Args:
            c_params: anisotropic params, orderd as [c12, c13, c21, c23, c31, c32, c44, c55, c66]
            stress:

        Returns:
            S1, S2, S3
    """
    T = np.array([[2, -1, -1, 0, 0, 0],
                  [-1, 2, -1, 0, 0, 0],
                  [-1, -1, 2, 0, 0, 0],
                  [0, 0, 0, 3, 0, 0],
                  [0, 0, 0, 0, 3, 0],
                  [0, 0, 0, 0, 0, 3],
                  ])/3

    C = np.array([[0, -c_params[0], -c_params[1], 0, 0, 0],
                  [-c_params[2], 0, -c_params[3], 0, 0, 0],
                  [-c_params[4], -c_params[5], 0, 0, 0, 0],
                  [0, 0, 0, c_params[6], 0, 0],
                  [0, 0, 0, 0, c_params[7], 0],
                  [0, 0, 0, 0, 0, c_params[8]],
                  ])

    s = C@T@stress

    H1 = (s[0] + s[1] + s[2])/3
    H2 = (s[5]**2 + s[4]**2 + s[3]**2 - s[1]*s[2] - s[2]*s[0] - s[0]*s[1])/3
    # H3 = (2*s[5]*s[4]*s[3] + s[0]*s[1]*s[2] - s[0]*s[5]**2 - s[1]*s[4]**2 - s[2]*s[3]**2)/2
    H3 = (2*s[5]*s[4]*s[3] + s[0]*s[1]*s[2] - s[0]*s[3]**2 - s[1]*s[4]**2 - s[2]*s[5]**2)/2

    p = (H1**2 + H2)
    q = (2*H1**3 + 3*H1*H2 + 2*H3)/2
    cos_theta = q/(p**1.5)
    if cos_theta > 1:
        # print(cos_theta)
        cos_theta = 1
    theta = math.acos(cos_theta)

    S1 = 2*math.sqrt(H1**2 + H2)*math.cos(theta/3) + H1
    S2 = 2*math.sqrt(H1**2 + H2)*math.cos(theta+4*math.pi/3) + H1
    S3 = 2*math.sqrt(H1**2 + H2)*math.cos(theta+2*math.pi/3) + H1

    return [S1, S2, S3]


def calc_phi(a, c_prime, c_Db_prime, stress):
    """
        Calculate eq_stress by yld

        Args:
            a: degree of a yld expression
            c_prime, c_Db_prime: anisotropic params,
                                 orderd as [c12, c13, c21, c23, c31, c32, c44, c55, c66]
            stress:

        Returns:
            eq_stress
    """
    S_prime = calc_S(c_prime, stress)
    S_Db_prime = calc_S(c_Db_prime, stress)
    phi = 0
    for i in range(3):
        for j in range(3):
            phi += (S_prime[i] - S_Db_prime[j])**a

    return phi


def calc_yld(a, c_prime, c_Db_prime, stress):
    phi = calc_phi(a, c_prime, c_Db_prime, stress)
    eq_stress = (phi/4)**(1/a)
    return eq_stress


def error_func(c_params):
    error = 0
    for data in datas:
        direction = (int(data["direction"])/10)*math.pi/180
        stress0 = float(data["yield_stress"])
        stress = [stress0*math.cos(direction)**2, stress0*math.sin(direction) ** 2,
                  0, 0, 0, stress0*math.sin(direction)*math.cos(direction)]
        pr_stress = calc_yld(6, c_params[:9], c_params[9:], np.array(stress))
        error += weight_s * (pr_stress/stress0 - 1)**2
    return error


def num_gradient(f, x):    # 勾配関数
    """
          Args:
              f: target func(error func)
              x: target params
          Returns:
              grad
    """
    h = 1e-4
    grad = np.zeros_like(x)    # xと同じ形状の配列で値がすべて 0
    for idx in range(x.size):    # x の次元分ループする。 (下記例は、5.0, 10.0 の 2 周ループ)
        idx_val = x[idx]
        x[idx] = idx_val + h    # f(x + h)の算出。
        fxh1 = f(x)

        x[idx] = idx_val - h    # f(x - h)の算出。
        fxh2 = f(x)

        grad[idx] = (fxh1 - fxh2) / (2 * h)
        x[idx] = idx_val    # 値をループ先頭の状態に戻す。
    return grad


with open('Datas/material_data.csv') as f:
    reader = csv.DictReader(f)
    datas = [row for row in reader]


datas.pop(-1)
weight_s = 1


c_values = np.array([0.5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

# 以下バッチ勾配降下法
tolerance = 0.00001
count = 0
learning_rate = 0.1
error = error_func(c_values)
print(error)

with open('Datas/error_func.csv', 'w') as f:
    writer = csv.writer(f)
    while error > tolerance:
        print(error)
        row = np.insert(c_values, 0, error)
        writer.writerow(row)
        grad = num_gradient(error_func, c_values)    # 上記 勾配関数「num_gradient」をコール
        c_values -= learning_rate * grad
        error = error_func(c_values)
        count += 1
