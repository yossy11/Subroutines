import numpy.linalg as LA
import numpy as np
# import matplotlib.patches as mpatches
import matplotlib.pyplot as plt


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


def draw_yield_surface(shear_stresses):
    fig = plt.figure(dpi=192)
    ax = fig.add_subplot(111)

    # 下軸と左軸をそれぞれ中央へもってくる
    ax.spines['bottom'].set_position(('data', 0))
    ax.spines['left'].set_position(('data', 0))

    # 上軸と右軸を表示しない
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    delta = 0.025
    xrange = np.arange(-2, 2, delta)
    yrange = np.arange(-2, 2, delta)
    X, Y = np.meshgrid(xrange, yrange)

    # 軸の設定
    ax.set_xlim(-2, 2)
    ax.set_ylim(-2, 2)

    # アスペクト比の設定
    ax.set_aspect('equal', adjustable='box')

    # 描画
    Z = np.empty_like(X)
    Zs = np.array([Z]*len(shear_stresses))
    for i in range(len(X)):
        for j in range(len(X)):
            for k, shear_stress in enumerate(shear_stresses):
                stress = make_stress_vector(X[i][j], Y[i][j], shear_stress)
                Zs[k][i][j] = calc_eqStress(stress, c_params, YLDM)

    for i in range(len(shear_stresses)):
        blue = format(int(255*i/len(shear_stresses)), '02x')
        color = f'#9900{blue}'
        ax.contour(X, Y, Zs[i], [1], colors=[color])

    # 実験値をplot
    # AA2090-T3
    exp_x = [1.0, 0.0, 1.035]
    exp_y = [0.0, 0.9102, 1.035]
    # DP980
    # exp_x = [1.0, 1.0608, 1.0785, 1.0513, 0.9873, 0.8302, 0.5775, 0.2794, 0.0]
    # exp_y = [0.0, 0.2652, 0.5402, 0.7885, 0.9859, 1.1054, 1.1524, 1.1159, 1.0311]
    exps, = plt.plot(exp_x, exp_y, marker='o', markersize=4,
                     color="black", linestyle='None', label="exp")

    # 凡例の設定
    # yld_patch = mpatches.Patch(color='#990000', linewidth=1, label=f'm={YLDM}')
    # ax.legend(handles=[(exps), yld_patch])
    # ax.legend(handles=[yld_patch])
    ax.legend(handles=[(exps)])
    ax.grid(color='black', linestyle='dotted', linewidth=1)
    fig.savefig("test.png")


if __name__ == "__main__":
    # shear_stresses = np.arange(0, 0.5, 0.05)
    shear_stresses = [0]
    draw_yield_surface(shear_stresses)
