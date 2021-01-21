# calculate yield stress and r_value using yld2000-2d
import math
import numpy as np
import numpy.linalg as LA

yld_T = np.array([[2, -1, -1, 0, 0, 0],
                  [-1, 2, -1, 0, 0, 0],
                  [-1, -1, 2, 0, 0, 0],
                  [0, 0, 0, 3, 0, 0],
                  [0, 0, 0, 0, 3, 0],
                  [0, 0, 0, 0, 0, 3]])
yld_T = yld_T/3.0

# AA2090-T3
# a_params = np.array([0.4865, 1.3783, 0.7536, 1.0246, 1.0363, 0.9036, 1.2321, 1.4858])
# YLDM = 8

# SPCD
a_params = np.array([0.9477, 1.1033, 0.7719, 0.9393, 0.9524, 0.8086, 0.9935, 1.1375])
YLDM = 6


def make_C_matrix(a_params):
    yld_C_prime = np.array([[a_params[0], 0, 0, 0, 0, 0],
                            [0, a_params[1], 0, 0, 0, 0],
                            [0, 0, 0, 0, 0, 0],
                            [0, 0, 0, a_params[6], 0, 0],
                            [0, 0, 0, 0, 0, 0],
                            [0, 0, 0, 0, 0, 0]])
    yld_C_Db_prime = np.array([[(4*a_params[4]-a_params[2])/3, (2*a_params[5]-2*a_params[3])/3, 0,
                                0, 0, 0],
                               [(2*a_params[2]-2*a_params[4])/3, (4*a_params[3]-a_params[5])/3, 0,
                                0, 0, 0],
                               [0, 0, 0, 0, 0, 0],
                               [0, 0, 0, a_params[7], 0, 0],
                               [0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0]])
    return yld_C_prime, yld_C_Db_prime


def calc_phi(stress, a_params, YLDM):
    yld_C_prime, yld_C_Db_prime = make_C_matrix(a_params)
    s_prime = yld_C_prime@yld_T@stress
    s_Db_prime = yld_C_Db_prime@yld_T@stress
    s_prime_mat = np.array([[s_prime[0], s_prime[3]], [s_prime[3], s_prime[1]]])
    s_Db_prime_mat = np.array([[s_Db_prime[0], s_Db_prime[3]], [s_Db_prime[3], s_Db_prime[1]]])
    pri_stress, _v = LA.eig(s_prime_mat)
    pri_Db_stress, _v = LA.eig(s_Db_prime_mat)
    phi1 = abs(pri_stress[0]-pri_stress[1])**YLDM
    phi2 = abs(2*pri_Db_stress[1]+pri_Db_stress[0])**YLDM + \
        abs(2*pri_Db_stress[0]+pri_Db_stress[1])**YLDM
    phi = phi1 + phi2
    return phi


def calc_angled_eqStress(angle, a_params, YLDM):
    if angle == "EB":
        angled_stress = np.array([0, 0, 1, 0, 0, 0])
    elif angle == "TDND":
        angled_stress = np.array([0, 0, 0, 0, 0, 1.0])
    elif angle == "NDRD":
        angled_stress = np.array([0, 0, 0, 0, 1.0, 0])
    elif angle == "TDND45":
        angled_stress = np.array([0, 0.5, 0.5, 0, 0, 0.5])
    elif angle == "NDRD45":
        angled_stress = np.array([0.5, 0, 0.5, 0, 0.5, 0])
    else:
        angle = math.pi*float(angle)/180.0
        angled_stress = np.array([math.cos(angle)**2, math.sin(angle)**2,
                                  0, math.sin(angle)*math.cos(angle), 0, 0])
    phi = calc_phi(angled_stress, a_params, YLDM)
    angled_eqStress = (2/phi)**(1/YLDM)
    return angled_eqStress


def calc_dphids(angle, a_params, YLDM):
    DELTAX = 1.0e-6
    dphids = np.zeros(6)
    if angle == "EB":
        angled_stress = np.array([-1/3.0, -1/3.0, 2/3.0, 0, 0, 0])
    elif angle == "TDND":
        angled_stress = np.array([0, 0, 0, 0, 0, 1.0])
    elif angle == "NDRD":
        angled_stress = np.array([0, 0, 0, 0, 1.0, 0])
    elif angle == "TDND45":
        angled_stress = np.array([0, 0.5, 0.5, 0, 0, 0.5])
    elif angle == "NDRD45":
        angled_stress = np.array([0.5, 0, 0.5, 0, 0.5, 0])
    else:
        angle = math.pi*float(angle)/180.0
        angled_stress = np.array([math.cos(angle)**2, math.sin(angle)**2,
                                  0, math.sin(angle)*math.cos(angle), 0, 0])
    phi = calc_phi(angled_stress, a_params, YLDM)
    for i in range(6):
        sub_angled_stress = angled_stress.copy()
        sub_angled_stress[i] += DELTAX
        sub_phi = calc_phi(sub_angled_stress, a_params, YLDM)
        dphids[i] = (sub_phi - phi)/DELTAX
    return dphids


def calc_angled_r(angle, a_params, YLDM):
    dphids = calc_dphids(angle, a_params, YLDM)
    if angle == "EB":
        r = dphids[1]/dphids[0]
    else:
        angle = math.pi*float(angle)/180.0
        r = (dphids[3]*math.cos(angle)*math.sin(angle)-dphids[0]*math.sin(angle)**2 -
             dphids[1]*math.cos(angle)**2)/(dphids[0]+dphids[1])
    return r


if __name__ == "__main__":
    angle = "90"
    print(calc_angled_eqStress(angle, a_params, YLDM))
    print(calc_angled_r(angle, a_params, YLDM))
