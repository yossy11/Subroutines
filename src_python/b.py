import math
import csv
import numpy as np
import numpy.linalg as LA


yld_T = np.array([[2, -1, -1, 0, 0, 0],
                  [-1, 2, -1, 0, 0, 0],
                  [-1, -1, 2, 0, 0, 0],
                  [0, 0, 0, 3, 0, 0],
                  [0, 0, 0, 0, 3, 0],
                  [0, 0, 0, 0, 0, 3]])
yld_T = yld_T/3.0
c_prime = np.array([-0.0698, 0.9364, 0.0791, 1.0030, 0.5247, 1.3631, 0.9543, 1.0690, 1.0237])
c_Db_prime = np.array([0.9811, 0.4767, 0.5753, 0.8668, 1.1450, -0.0792, 1.4046, 1.1471, 1.0516])
YLDM = 8


def make_C_matrix(c_params):
    C = np.array([[0, -1.0*c_params[0], -1.0*c_params[1], 0, 0, 0],
                  [-1.0*c_params[2], 0, -1.0*c_params[3], 0, 0, 0],
                  [-1.0*c_params[4], -1.0*c_params[5], 0, 0, 0, 0],
                  [0, 0, 0, c_params[6], 0, 0],
                  [0, 0, 0, 0, c_params[7], 0],
                  [0, 0, 0, 0, 0, c_params[8]]])
    return C


def calc_phi(stress, c_prime, c_Db_prime, YLDM):
    yld_C_prime = make_C_matrix(c_prime)
    yld_C_Db_prime = make_C_matrix(c_Db_prime)
    s_prime = yld_C_prime@yld_T@stress
    s_Db_prime = yld_C_Db_prime@yld_T@stress
    s_prime_mat = np.array([[s_prime[0], s_prime[3], s_prime[4]],
                            [s_prime[3], s_prime[1], s_prime[5]],
                            [s_prime[4], s_prime[5], s_prime[2]]])
    s_Db_prime_mat = np.array([[s_Db_prime[0], s_Db_prime[3], s_Db_prime[4]],
                               [s_Db_prime[3], s_Db_prime[1], s_Db_prime[5]],
                               [s_Db_prime[4], s_Db_prime[5], s_Db_prime[2]]])
    pri_stress, _v = LA.eig(s_prime_mat)
    pri_Db_stress, _v = LA.eig(s_Db_prime_mat)
    phi = 0.0
    for i in range(3):
        for j in range(3):
            phi += abs(pri_stress[i] - pri_Db_stress[j])**YLDM

    return phi


def calc_angled_eqStress(angle, c_prime, c_Db_prime, YLDM):
    angled_stress = np.array([math.cos(angle)**2, math.sin(angle)**2,
                              0, math.sin(angle)*math.cos(angle), 0, 0])
    phi = calc_phi(angled_stress, c_prime, c_Db_prime, YLDM)
    angled_eqStress = (4/phi)**(1/YLDM)
    return angled_eqStress


def error_func(exp_data, c_prime, c_Db_prime, YLDM, wp, wb):
    error = 0
    for data in exp_data:
        if not data["orientation"].isdecimal():
            stress = np.array([0, 0, 1, 0, 0, 0])
            phi = calc_phi(stress, c_prime, c_Db_prime, YLDM)
            predicted_stress = (4/phi)**(1/YLDM)
            error += wb*(predicted_stress/float(data["normalized_yield_stress"]) - 1.0)**2
            continue
        predicted_stress = calc_angled_eqStress(
            math.pi*float(data["orientation"])/180.0, c_prime, c_Db_prime, YLDM)
        print(predicted_stress)
        error += wp*(predicted_stress/float(data["normalized_yield_stress"])-1.0)**2
    return error


if __name__ == "__main__":
    with open("Datas/AA2090-T3.csv") as f:
        reader = csv.DictReader(f)
        exp_data = [row for row in reader]
    print(error_func(exp_data, c_prime, c_Db_prime, YLDM, 1.0, 0.01))
