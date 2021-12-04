"""
Основные уравнения динамики КА и их решение
"""

import numpy as np
import matplotlib.pyplot as plt
import flywheel_engine as dm


# Решение дифференциальных уравнений в конкретный момент времени
# На вход получает моменты, кинетические моменты ДМ и проекции
# угловой скорости на ССК
def f(t, w, M, H, HH, I):
    ww = np.array([[0.0], [0.0], [0.0], [0.0], [0.0], [0.0]])
    a = 1/(I[0,0]*I[1,1]*I[2,2] - I[0,0]*I[1,2]**2 - I[1,1]*I[0,2]**2 - I[2,2]*I[0,1]**2 - 2*I[0,1]*I[0,2]*I[1,2])
    sigma = [0, 0, 0]
    sigma[0] = M[0] - (HH[0] + H[2]*w[1] - H[1]*w[2])
    sigma[1] = M[1] - (HH[1] + H[0]*w[2] - H[2]*w[0])
    sigma[2] = M[2] - (HH[2] + H[1]*w[0] - H[0]*w[1])
    ww[0,0] = a*((I[1,1]*I[2,2]-I[1,2]**2)*sigma[0] + (I[2,2]*I[0,1]+I[0,2]*I[1,2])*sigma[1] +
               + (I[1,1]*I[0,2]+I[0,1]*I[1,2])*sigma[2])
    ww[1,0] = a*((I[2,2]*I[0,1]+I[0,2]*I[1,2])*sigma[0] + (I[0,0]*I[2,2]-I[0,2]**2)*sigma[1] +
               + (I[0,0]*I[1,2]+I[0,1]*I[0,2])*sigma[2])
    ww[2,0] = a*((I[1,1]*I[0,2]+I[0,1]*I[1,2])*sigma[0] + (I[0,0]*I[1,2]+I[0,1]*I[0,2])*sigma[1] +
               + (I[0,0]*I[1,1]-I[0,1]**2)*sigma[2])
    ww[3,0] = w[0]
    ww[4,0] = w[1]
    ww[5,0] = w[2]
    return ww


# Решение системы диффуров на определенном промежутке времени
# методом Рунге-Кутты
def runge_kutta(t_span, dt, x0, M0, I):

    # инициализация моделей ДМ
    dm_all = []
    for j in range(4):
        dm_all.append(dm.Flywheel(0,0,0,0))

    # начальные условия
    t_begin = t_span[0]
    t_end = t_span[1]
    k = 0
    x = x0

    alpha = 20.0 / 180 * np.pi
    t_curr = t_begin + dt
    h = dt
    t = [t_curr]

    H_out = np.array([[0.0], [0.0], [0.0], [0.0], [0.0], [0.0]])

    k1 = np.array([[0.0], [0.0], [0.0], [0.0], [0.0], [0.0]])
    k2 = np.array([[0.0], [0.0], [0.0], [0.0], [0.0], [0.0]])
    k3 = np.array([[0.0], [0.0], [0.0], [0.0], [0.0], [0.0]])
    k4 = np.array([[0.0], [0.0], [0.0], [0.0], [0.0], [0.0]])

    # интегрирование
    while t_curr <= t_end:

        t.append(t_curr)

        # текущие значения скоростей и углов
        right_column = x[:, k:]
        angles = right_column[3:,0]
        vel = right_column[0:3,0]
        # right_column = right_column.T

        # пока все Н считать одинаковыми
        # функция управления; пока для тестовых целей просто синусоида

        [Md, H, HH] = dm.get_all(angles, vel, dm_all)

        for j in range(3):
            M[j] = M0[j] + Md[j]

        k1 = f(t_curr, right_column, M, H, HH, I)
        k2 = f(t_curr, right_column + h * k1 / 2, M, H, HH, I)
        k3 = f(t_curr, right_column + h * k2 / 2, M, H, HH, I)
        k4 = f(t_curr, right_column + h * k3, M, H, HH, I)

        H_out_curr = np.array([[H[0]], [H[1]], [H[2]], [HH[0]], [HH[1]], [HH[2]]])
        right_column = right_column + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        x = np.append(x, right_column, axis=1)
        H_out = np.append(H_out, H_out_curr, axis=1)
        k += 1
        t_curr += h
        """for i in range(3, 6):
            if x[i, k] >= 2*np.pi:
                x[i, k] -= 2*np.pi
            if x[i, k] <= -2*np.pi:
                x[i, k] += 2*np.pi"""
        if t_curr % 1000 == 0:
            print('t_curr = ', t_curr)

    return [t, x, H_out]


if __name__ == '__main__':

    # время интегрирования
    t_span = [0, 300]

    # такт вычислений
    dt = 0.25

    # внешний возмущающий постоянный момент, угловые скорости и углы (в град/с и град)
    M = [0.0, 0.0, 0.0]
    w0 = np.array([[0.0], [0.0], [0.0], [-180.0], [-180.0], [-180.0]])
    w0 = w0 / 180 * np.pi

    # тензор инерции (данные из документации, не менять)
    I = np.array([[3337.94, 25.6, 3.2], [25.6, 9283.9, 19.6], [3.2, 19.6, 10578.5]])

    [t, x, H_out] = runge_kutta(t_span, dt, w0, M, I)
    # sol = solve_ivp(f, t_span, w0, 'RK45', t_eval, args=(M, H, HH, I), dense_output=True)

    # t = t_eval
    # z = sol.sol(t)
    # z = z.T

    # отображение графиков
    n = len(t)
    for i in range(n):
        x[0, i] = x[0, i] / np.pi * 180
        x[1, i] = x[1, i] / np.pi * 180
        x[2, i] = x[2, i] / np.pi * 180
        x[3, i] = x[3, i] / np.pi * 180
        x[4, i] = x[4, i] / np.pi * 180
        x[5, i] = x[5, i] / np.pi * 180

    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4,1)

    ax1.plot(t, x[3, :], label="gamma, град")
    ax1.plot(t, x[4, :], label="theta, град")
    ax1.plot(t, x[5, :], label="psi, град")
    ax1.set_title("Изменение углов")
    ax1.legend()
    ax1.grid(True)

    ax2.plot(t, x[0,:], label="wx, град/с")
    ax2.plot(t, x[1,:], label="wy, град/с")
    ax2.plot(t, x[2,:], label="wz, град/с")
    ax2.set_title("Изменение угловых скоростей")
    ax2.legend()
    ax2.grid(True)

    ax3.plot(t, H_out[0, :], label="Hx")
    ax3.plot(t, H_out[1, :], label="Hy")
    ax3.plot(t, H_out[2, :], label="Hz")
    ax3.set_title("Изменение кинетического момента ДМ")
    ax3.legend()
    ax3.grid(True)

    ax4.plot(t, H_out[3, :], label="HHx")
    ax4.plot(t, H_out[4, :], label="HHy")
    ax4.plot(t, H_out[5, :], label="HHz")
    ax4.set_title("Изменение производной кинетического момента ДМ")
    ax4.legend()
    ax4.grid(True)

    fig.tight_layout()
    plt.show()
