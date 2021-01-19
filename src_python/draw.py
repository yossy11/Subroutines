# import pickle
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
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


def calc_eqStress(x, y, shear):
    stress = np.array([x, y, 0, shear, 0, 0])
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

    eq_stress = (phi/4.0)**(1/9.0)
    return eq_stress


fig = plt.figure(dpi=192)
ax = fig.add_subplot(111)
ax.spines['left'].set(position='zero')
ax.spines['bottom'].set(position='zero')
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
        # Z[i][j] = calc_eqStress(X[i][j], Y[i][j], 0)
        for k in range(10):
            Zs[k][i][j] = calc_eqStress(X[i][j], Y[i][j], k*0.05)

# with open("Datas/Zs.pkl", 'wb') as f:
#     pickle.dump(Zs, f)

# with open("Datas/Zs.pkl", "rb") as f:
#     Zs = pickle.load(f)

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
plt.savefig("figure9.png")  # -----(2)
