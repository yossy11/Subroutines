import numpy.linalg as LA
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt


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


yld_T = np.array([[2, -1, -1, 0, 0, 0],
                  [-1, 2, -1, 0, 0, 0],
                  [-1, -1, 2, 0, 0, 0],
                  [0, 0, 0, 3, 0, 0],
                  [0, 0, 0, 0, 3, 0],
                  [0, 0, 0, 0, 0, 3]])
yld_T = yld_T/3.0
c_params = np.array([-0.0698, 0.9364, 0.0791, 1.0030, 0.5247, 1.3631, 0.9543, 1.0690, 1.0237,
                     0.9811, 0.4767, 0.5753, 0.8668, 1.1450, -0.0792, 1.4046, 1.1471, 1.0516])
c_params = np.array([-0.11714269625987313, 0.4677946760370396, -0.6329214118208218,
                     0.1958292433240264, 0.20155498293410648, -0.5137137661110223,
                     -1.1624636739743532, 1.082419111390934, 1.1164359346997557,
                     1.1246439061566151, -0.3140781585816851, 0.2885538299610835,
                     1.19435689224697, 1.3081518030712693, 0.5495801370099768,
                     -1.2006650067491322, 1.090309427765219, 1.1181902696635893])
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


def calc_eqStress(stress, c_params, YLDM):
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

    eq_stress = (phi/4.0)**(1/YLDM)
    return eq_stress


def make_stress_vector(x, y, shear):
    stress = np.array([x, y, 0, shear, 0, 0])
    return stress


fig = plt.figure(dpi=192)
ax = fig.add_subplot(111)

# 下軸と左軸をそれぞれ中央へもってくる
ax.spines['bottom'].set_position(('data', 0))
ax.spines['left'].set_position(('data', 0))

# 上軸と右軸を表示しない
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

# ax.spines['left'].set(position='zero')
# ax.spines['bottom'].set(position='zero')
delta = 0.025
xrange = np.arange(-2, 2, delta)
yrange = np.arange(-2, 2, delta)
X, Y = np.meshgrid(xrange, yrange)

# 軸の設定
plt.axis([-2, 2, -2, 2])
plt.gca().set_aspect('equal', adjustable='box')

# 描画
Z = np.empty_like(X)
Zs = np.array([Z]*10)
for i in range(len(X)):
    for j in range(len(X)):
        stress = make_stress_vector(X[i][j], Y[i][j], 0)
        Z[i][j] = calc_eqStress(stress, c_params, YLDM)
        # for k in range(10):
        #     stress = make_stress_vector(X[i][j], Y[i][j], k*0.05)
        #     Zs[k][i][j] = calc_eqStress(stress, c_params, YLDM)


# contours = plt.contour(X, Y, Z, [1], colors=["red"])

for i in range(10):
    blue = format(int(255*i/10), '02x')
    color = f'#9900{blue}'
    plt.contour(X, Y, Zs[i], [1], colors=[color])

# exp_x = [-1.01739, 0.01144, 1.00051, 1.04027, 0.0015, -0.98757]
# exp_y = [-1.02588, -0.90504, 0.00771, 1.03052, 0.9104, 0.00771, ]
# exps, = plt.plot(exp_x, exp_y, marker='o', markersize=4,
#                  color="black", linestyle='None', label="exp")


yld_patch = mpatches.Patch(color='#990000', linewidth=1, label='m=9')
# plt.legend(handles=[(exps), yld_patch])
plt.legend(handles=[yld_patch])
plt.grid(color='black', linestyle='dotted', linewidth=1)
# ax.clabel(contours, ["Yld2004-18P"], inline=True, fontsize=10)
# plt.legend(bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0, fontsize=18)
plt.show()
plt.savefig("figuretest.png")  # -----(2)
