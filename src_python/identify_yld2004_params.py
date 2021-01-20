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
c_params = np.array([-0.0698, 0.9364, 0.0791, 1.0030, 0.5247, 1.3631, 0.9543, 1.0690, 1.0237,
                     0.9811, 0.4767, 0.5753, 0.8668, 1.1450, -0.0792, 1.4046, 1.1471, 1.0516])
YLDM = 8


def make_C_matrix(c_params):
    yld_C_prime = np.array([[0, -1.0*c_params[0], -1.0*c_params[1], 0, 0, 0],
                            [-1.0*c_params[2], 0, -1.0*c_params[3], 0, 0, 0],
                            [-1.0*c_params[4], -1.0*c_params[5], 0, 0, 0, 0],
                            [0, 0, 0, c_params[6], 0, 0],
                            [0, 0, 0, 0, c_params[7], 0],
                            [0, 0, 0, 0, 0, c_params[8]]])
    yld_C_Db_prime = np.array([[0, -1.0*c_params[9], -1.0*c_params[10], 0, 0, 0],
                               [-1.0*c_params[11], 0, -1.0*c_params[12], 0, 0, 0],
                               [-1.0*c_params[13], -1.0*c_params[14], 0, 0, 0, 0],
                               [0, 0, 0, c_params[15], 0, 0],
                               [0, 0, 0, 0, c_params[16], 0],
                               [0, 0, 0, 0, 0, c_params[17]]])
    return yld_C_prime, yld_C_Db_prime


def calc_phi(stress, c_params, YLDM):
    yld_C_prime, yld_C_Db_prime = make_C_matrix(c_params)
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


def calc_angled_eqStress(angle, c_params, YLDM):
    if angle == "EB":
        angled_stress = np.array([0, 0, 1, 0, 0, 0])
    else:
        angle = math.pi*float(angle)/180.0
        angled_stress = np.array([math.cos(angle)**2, math.sin(angle)**2,
                                  0, math.sin(angle)*math.cos(angle), 0, 0])
    phi = calc_phi(angled_stress, c_params, YLDM)
    angled_eqStress = (4/phi)**(1/YLDM)
    return angled_eqStress


def error_func(exp_data, c_params, YLDM, wp, wb):
    error = 0
    for data in exp_data:
        predicted_stress = calc_angled_eqStress(data["orientation"], c_params, YLDM)
        if data["orientation"] == "EB":
            error += wb*(predicted_stress/float(data["normalized_yield_stress"])-1.0)**2
        else:
            error += wp*(predicted_stress/float(data["normalized_yield_stress"])-1.0)**2
    return error


def calc_dsdc(angle, c_params, YLDM):
    DELTAX = 1.0e-6
    dsdc = np.ones(18)
    predicted_stress = calc_angled_eqStress(angle, c_params, YLDM)
    for i in range(18):
        sub_c_params = c_params.copy()
        sub_c_params[i] += DELTAX
        sub_predicted_stress = calc_angled_eqStress(angle, sub_c_params, YLDM)
        dsdc[i] = (sub_predicted_stress - predicted_stress)/DELTAX
    return dsdc


def calc_dphids(angle, c_params, YLDM):
    DELTAX = 1.0e-6
    dphids = np.zeros(6)
    if angle == "EB":
        angled_stress = np.array([-1/3.0, -1/3.0, 2/3.0, 0, 0, 0])
    else:
        angle = math.pi*float(angle)/180.0
        angled_stress = np.array([math.cos(angle)**2, math.sin(angle)**2,
                                  0, math.sin(angle)*math.cos(angle), 0, 0])
    phi = calc_phi(angled_stress, c_params, YLDM)
    for i in range(6):
        sub_angled_stress = angled_stress.copy()
        sub_angled_stress[i] += DELTAX
        sub_phi = calc_phi(sub_angled_stress, c_params, YLDM)
        dphids[i] = (sub_phi - phi)/DELTAX
    return dphids


def calc_dfds(angle, c_params, YLDM):
    dphids = calc_dphids(angle, c_params, YLDM)
    angled_eqStress = calc_angled_eqStress(angle, c_params, YLDM)
    dfds = (dphids*angled_eqStress**(1-YLDM))/(4*YLDM)
    return dfds


def calc_angled_r(angle, c_params, YLDM):
    dphids = calc_dphids(angle, c_params, YLDM)
    dfds = calc_dfds(angle, c_params, YLDM)
    if angle == "EB":
        r = dphids[1]/dphids[0]
    else:
        angle = math.pi*float(angle)/180.0
        r = (dfds[3]*math.cos(angle)*math.sin(angle)-dfds[0]*math.sin(angle)**2 -
             dfds[1]*math.cos(angle)**2)/(dfds[0]+dfds[1])
    return r


def calc_error_gradient(exp_data, c_params, YLDM, wp, wb):
    gradient = np.zeros(18)
    for data in exp_data:
        predicted_stress = calc_angled_eqStress(data["orientation"], c_params, YLDM)
        exp_stress = float(data["normalized_yield_stress"])
        dsdc = calc_dsdc(data["orientation"], c_params, YLDM)
        if data["orientation"] == "EB":
            weight = wb
        else:
            weight = wp
        gradient += (weight*2.0*(predicted_stress/exp_stress - 1.0)/exp_stress)*dsdc
    return gradient


def gradient_descent(exp_data, YLDM, wp, wb):
    learning_rate = 1.0e-1
    c_params = np.ones(18)
    with open('Datas/error.csv', 'w') as f:
        writer = csv.writer(f)
        header = ["iterationNum", "error", "c_params"]
        writer.writerow(header)
        for i in range(30000):
            error = error_func(exp_data, c_params, YLDM, wp, wb)
            writer.writerow([i, error, ' '.join(map(str, c_params))])
            print(i, error)
            gradient = calc_error_gradient(exp_data, c_params, YLDM, wp, wb)
            c_params -= learning_rate*gradient
    return c_params


if __name__ == "__main__":
    with open("Datas/AA2090-T3.csv") as f:
        reader = csv.DictReader(f)
        exp_data = [row for row in reader]
    wp = 1.0
    wb = 0.01
    # error = error_func(exp_data, c_params, YLDM, wp, wb)
    # print(error)
    c_params = gradient_descent(exp_data, YLDM, wp, wb)
    print(c_params)
