# calculate the appropriate scale for DELTAX

import math
import csv
import random
import glob
import numpy as np

HILL_F = 0.25216953733566727
HILL_G = 0.8254230293025175
HILL_H = 0.17457697069748246
HILL_L = 2.2380520016508463
HILL_M = 2.2380520016508463
HILL_N = 2.2380520016508463
# HILL_F = 1.0
# HILL_G = 1.0
# HILL_H = 1.0
# HILL_L = 3.0
# HILL_M = 3.0
# HILL_N = 3.0


def hill(sxx, syy, szz, sxy, sxz, syz):
    numerator = HILL_F*(syy - szz)**2 + HILL_G*(szz - sxx)**2 + \
        HILL_H*(sxx - syy)**2 + 2*HILL_L*syz**2 + 2*HILL_M*sxz**2 + 2*HILL_N*sxy**2
    return math.sqrt(numerator*1.5/(HILL_F+HILL_G+HILL_H))


def diff_hill(DELTAX, sxx, syy, szz, sxy, sxz, syz):
    # analytical method
    equivalent_stress = hill(sxx, syy, szz, sxy, sxz, syz)
    multiplier = 0.75/((HILL_F+HILL_G+HILL_H)*equivalent_stress)
    dsxx = multiplier * (-2*HILL_G*(szz-sxx)+2*HILL_H*(sxx-syy))
    dsyy = multiplier * (2*HILL_F*(syy-szz)-2*HILL_H*(sxx-syy))
    dszz = multiplier * (-2*HILL_F*(syy-szz)+2*HILL_G*(szz-sxx))
    dsxy = multiplier * (4*HILL_N*sxy)
    dsxz = multiplier * (4*HILL_M*sxz)
    dsyz = multiplier * (4*HILL_L*syz)
    result = [dsxx, dsyy, dszz, dsxy, dsxz, dsyz]

    # numerical method
    init_value = [sxx, syy, szz, sxy, sxz, syz]
    for i in range(6):
        value = init_value.copy()
        value[i] += DELTAX
        result.append((hill(*value) - equivalent_stress) / DELTAX)

    for i in range(6):
        result.append(result[i] - result[i+6])

    return result


if __name__ == "__main__":
    for i in range(8):
        DELTAX = 1.0*10**(-4-i)  # 1.0e-4 ~ 1.0e-11
        file_name = f'Datas/{DELTAX:.0e}.csv'
        with open(file_name, 'w') as f:
            writer = csv.writer(f)
            count = 0
            while count < 100:
                init_value = []
                for j in range(6):
                    init_value.append(450*random.random())
                writer.writerow(diff_hill(DELTAX, *init_value))
                count += 1

    with open("Datas/diff_of_diffs.csv", 'w') as f:
        writer = csv.writer(f)
        file_list = glob.glob('Datas/1e-*.csv')
        for file_name in file_list:
            print(file_name)
            with open(file_name) as f:
                reader = csv.reader(f)
                datas = [row for row in reader]

            new_datas = []
            for data in datas:
                data = [float(i) for i in data]
                new_datas.append(data)
            new_datas = np.array(new_datas)
            result = np.zeros(6)
            for data in new_datas:
                result += data[12:]

            result /= len(new_datas)
            print(result)
            writer.writerow([file_name, *result])
