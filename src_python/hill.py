import math
import csv
import random
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
    # print(diff_hill(1.0e-6, 300.0, 200.0, 100.0, 100.0, 100.0, 100.0))
    print(diff_hill(1.0e-6, 84.5349, 113.6041, -1.2214, -118.2640, 9.8132, -11.6963))

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
