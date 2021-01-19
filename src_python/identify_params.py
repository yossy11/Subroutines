import numpy as np
import numpy.linalg as LA


yld_C_prime = np.array([[0, 0.0698, -0.9364, 0, 0, 0],
                        [-0.0791, 0, -1.0030, 0, 0, 0],
                        [-0.5247, -1.3631, 0, 0, 0, 0],
                        [0, 0, 0, 0.9543, 0, 0],
                        [0, 0, 0, 0, 1.0690, 0],
                        [0, 0, 0, 0, 0, 1.0237]])


yld_C_Db_prime = np.array([[0, -0.9811, -0.4767, 0, 0, 0],
                           [-0.5753, 0, -0.8668, 0, 0, 0],
                           [-1.1450, 0.0792, 0, 0, 0, 0],
                           [0, 0, 0, 1.4046, 0, 0],
                           [0, 0, 0, 0, 1.1471, 0],
                           [0, 0, 0, 0, 0, 1.0516]])

yld_T = np.array([[2, -1, -1, 0, 0, 0],
                  [-1, 2, -1, 0, 0, 0],
                  [-1, -1, 2, 0, 0, 0],
                  [0, 0, 0, 3, 0, 0],
                  [0, 0, 0, 0, 3, 0],
                  [0, 0, 0, 0, 0, 3]])

yld_T = yld_T/3.0
yld_stress = 646*(0.025**0.227)


def calc_eqStress(YLDM, yld_C_prime, yld_C_Db_prime, stress):
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
            phi += abs(pri_stress[i] - pri_Db_stress[j])**9

    eq_stress = (phi/4.0)**(1/YLDM)
    return eq_stress
