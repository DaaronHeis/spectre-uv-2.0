"""
    Основные уравнения динамики КА и их решение
"""

# TODO: перейти с quaternion на quaternionic

import numpy as np
import quaternion as qt
import matplotlib.pyplot as plt
import flywheel_engine as dm
import control_unit as ctrl


# Решение дифференциальных уравнений в конкретный момент времени
# На вход получает моменты, кинетические моменты ДМ и проекции
# угловой скорости на ССК
def f(t, w, M, H, HH, I):
    """
        Диффуры динамики КА
        Представляют из себя следующую систему:
            w' = A * [a] @ [sigma]
            phi' = w
            [a] и [sigma] - матрица 3х3 и 3х1 соответственно
            [a] - матрица моменттов инерции
            [sigma] - вектор-столбец действующих моментов
    """
    ww = np.array([[0.0], [0.0], [0.0], [0.0], [0.0], [0.0]])
    A = 1/(I[0,0]*I[1,1]*I[2,2] - I[0,0]*I[1,2]**2 - I[1,1]*I[0,2]**2 - I[2,2]*I[0,1]**2 - 2*I[0,1]*I[0,2]*I[1,2])
    sigma = np.array([[0.], [0.], [0.]])
    sigma[0,0] = M[0] - (HH[0] + H[2]*w[1] - H[1]*w[2])
    sigma[1,0] = M[1] - (HH[1] + H[0]*w[2] - H[2]*w[0])
    sigma[2,0] = M[2] - (HH[2] + H[1]*w[0] - H[0]*w[1])
    a = np.array([[I[1,1]*I[2,2]-I[1,2]**2, I[2,2]*I[0,1]+I[0,2]*I[1,2], I[1,1]*I[0,2]+I[0,1]*I[1,2]],
                  [I[2,2]*I[0,1]+I[0,2]*I[1,2], I[0,0]*I[2,2]-I[0,2]**2, I[0,0]*I[1,2]+I[0,1]*I[0,2]],
                  [I[1,1]*I[0,2]+I[0,1]*I[1,2], I[0,0]*I[1,2]+I[0,1]*I[0,2], I[0,0]*I[1,1]-I[0,1]**2]])
    ww[0:3, :] = A * (a @ sigma)
    """
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
    """
    ww[3,0] = w[0]
    ww[4,0] = w[1]
    ww[5,0] = w[2]
    return ww


# Решение системы диффуров на определенном промежутке времени
# методом Рунге-Кутты
def runge_kutta(t_span, dt, angles_0, angles_end, vel_0, vel_end, M0, I):

    # инициализация моделей ДМ
    dm_all = []
    # TODO: учесть ограничения на скорость из-за астродатчика (не более 0.2 град/с и 0.1 град/с^2 отн. ЭСК)
    sigma_max = 0.2
    w_bw = 5
    for j in range(4):
        dm_all.append(dm.Flywheel(0,0,0,0,w_bw))

    # инициализация модуля управления
    vel = vel_0
    angles = angles_0
    omega_pr = vel_end
    gamma_pr = angles_end
    # TODO: разобраться с начальными параметрами и их корректным представлением
    l_0 = ctrl.from_euler_to_quat(angles, 'YZXr')
    print('l_0 = ', l_0)
    l_k = ctrl.from_euler_to_quat(gamma_pr, 'YZXr')

    # пока итоговая ориентация задается здесь
    # l_pr = l_k.inverse() * l_0
    l_pr = qt.quaternion(1, 0, 0, 0)
    print('l_pr = ', l_pr)
    l_cur = l_0
    l_delta = l_pr.inverse() * l_cur
    print('l_delta = ', l_delta)
    l_delta_out = np.array([qt.as_float_array(l_delta)])
    l_cur_out = np.array([qt.as_float_array(l_cur)])
    angles = np.array([2*l_delta.w*l_delta.x, 2*l_delta.w*l_delta.y, 2*l_delta.w*l_delta.z])

    ctrl_unit = ctrl.ControlUnit(l_pr, l_cur, l_delta, vel, angles, omega_pr, I, w_bw, sigma_max, dt)

    # начальные условия
    t_begin = t_span[0]
    t_end = t_span[1]
    k = 0
    x = np.array([[vel[0]], [vel[1]], [vel[2]],
                  [angles[0]], [angles[1]], [angles[2]]])

    h = dt
    t = [t_begin]
    t_curr = t_begin + h

    H_out = np.array([[0.0], [0.0], [0.0], [0.0], [0.0], [0.0], [0.0]])
    sigma_out = np.array([[0.0], [0.0], [0.0]])

    # интегрирование
    while t_curr <= t_end:

        t.append(t_curr)

        # текущие значения скоростей и углов
        right_column = x[:, k:]
        angles = right_column[3:,0]
        vel = right_column[0:3,0]
        # right_column = right_column.T

        # расчет управляющего момента
        sigma = ctrl_unit.update(vel, t_curr)
        sigma = [sigma[0], sigma[1], sigma[2]]
        npsigma = np.array([[sigma[0]], [sigma[1]], [sigma[2]]])
        sigma_out = np.append(sigma_out, npsigma, axis=1)

        # получение текущих параметров модуля управления
        [a, b, c, d] = ctrl_unit.get_parameters()
        a = qt.as_float_array(a)
        l_cur = np.array([[a[0], a[1], a[2], a[3]]])
        b = qt.as_float_array(b)
        l_delta = np.array([[b[0], b[1], b[2], b[3]]])
        l_cur_out = np.append(l_cur_out, l_cur, axis=0)
        l_delta_out = np.append(l_delta_out, l_delta, axis=0)
        angles = c #np.array([c[0], c[2], c[1]])

        # расчет динамического момента комплекса двигателей-маховиков
        Md = dm.update_block(sigma, dm_all)

        # получение параметров комплекса двигателей-маховиков
        dm_param = dm.get_all(dm_all)
        H_dm = []
        HH_dm = []
        w_self_dm = []
        for param in dm_param:
            H_dm.append(param[1])
            HH_dm.append(param[2])
            w_self_dm.append(param[3])
        H = dm.from_dm_to_xyz(H_dm)
        HH = dm.from_dm_to_xyz(HH_dm)

        # расчет действующего момента
        for j in range(3):
            M[j] = M0[j] + Md[j]

        # коэффициенты для Рунге-Кутты
        # right_column[3:, 0] = angles
        k1 = f(t_curr, right_column, M, H, HH, I)
        k2 = f(t_curr, right_column + h * k1 / 2, M, H, HH, I)
        k3 = f(t_curr, right_column + h * k2 / 2, M, H, HH, I)
        k4 = f(t_curr, right_column + h * k3, M, H, HH, I)

        # TODO: написать отдельную функцию по формированию выходных массивов
        H_out_curr = np.array([[H[0]], [H[1]], [H[2]], [HH_dm[0]], [HH_dm[1]], [HH_dm[2]], [HH_dm[3]]])
        right_column = right_column + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        right_column[3:, 0] = angles
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

    # final_param = [t, x, H_out, l_cur_out, l_delta_out, sigma_out]
    return [t, x, H_out, l_cur_out, l_delta_out, sigma_out]


if __name__ == '__main__':

    # время интегрирования
    t_span = [0, 1820]

    # такт вычислений
    dt = 0.25

    # коэффициент перехода от градусов к радианам
    k = np.pi / 180

    """
        Небольшое пояснение по поводу физического смысла вектора состояния
        Вектор состояния состоит из трех угловых скоростей и трех углов (по каждой оси)
        Углы получаются путем интегрирования угловой скорости при решении 
        кинематических уравнений, а значит, они являются углами поворота вокруг каждой оси
        Скорости измеряются БИУС, которые, по идее, измеряют их относительно ИСК, а значит, и углы
        будут являться углами вращения вокруг каждой из осей ИСК
        Изначальные и конечные углы задаются поворотами в порядке 'YZXr', то есть сначала
        поворот по тангажу (ось Y), затем по рысканью (вокруг оси Z'), и затем по крену (вокруг оси X'')  
    
        Значения внутри векторов - в град/с и град, но сразу же после этого они
        пересчитываются в радианы    
    """
    angles_0 = np.array([30.0, 30.0, 30.0]) * k
    angles_end = np.array([0.0, 0.0, 0.0]) * k
    vel_0 = np.array([0.0, 0.0, 0.0]) * k
    vel_end = np.array([0.0, 0.0, 0.0]) * k

    # внешний возмущающий постоянный момент
    M = [0.0, 0.0, 0.0]

    # с этого момента все углы в радианах

    # тензор инерции (данные из документации, не менять)
    I = np.array([[3337.94, 25.6, 3.2], [25.6, 9283.9, 19.6], [3.2, 19.6, 10578.5]])

    """
            !!!ВАЖНО!!!
        Выходные углы - это углы, которые получаются из текущего кватерниона ориентации
        Они (по идее) отображают углы, на которые нужно повернутся, чтобы прийти к итоговой ориентации
    """
    [t, x, H_out, l_cur, l_delta, sigma_out] = runge_kutta(t_span, dt, angles_0, angles_end, vel_0, vel_end, M, I)

    # отображение графиков
    n = len(t)
    x[:, :] = x[:, :] / k

    angles = x[3:6, :]
    vel = x[0:3, :]

    # с этого момента все углы в градусах

    fig1, (ax1, ax2) = plt.subplots(2,1)

    ax1.plot(t, angles[0, :], label="gamma, град")
    ax1.plot(t, angles[1, :], label="theta, град")
    ax1.plot(t, angles[2, :], label="psi, град")
    ax1.set_title("Изменение углов")
    ax1.legend()
    ax1.grid(True)

    ax2.plot(t, vel[0,:], label="wx, град/с")
    ax2.plot(t, vel[1,:], label="wy, град/с")
    ax2.plot(t, vel[2,:], label="wz, град/с")
    ax2.set_title("Изменение угловых скоростей")
    ax2.legend()
    ax2.grid(True)

    fig2, (ax1, ax2) = plt.subplots(2,1)

    ax1.plot(t, H_out[0, :], label="Hx")
    ax1.plot(t, H_out[1, :], label="Hy")
    ax1.plot(t, H_out[2, :], label="Hz")
    ax1.set_title("Изменение кинетического момента в проекциях")
    ax1.legend()
    ax1.grid(True)

    ax2.plot(t, H_out[3, :], label="H1")
    ax2.plot(t, H_out[4, :], label="H2")
    ax2.plot(t, H_out[5, :], label="H3")
    ax2.plot(t, H_out[6, :], label="H4")
    ax2.set_title("Изменение кинетического момента каждого ДМ")
    ax2.legend()
    ax2.grid(True)

    fig3, (ax1, ax2) = plt.subplots(2,1)

    ax1.plot(t, sigma_out[0, :], label="sigma_x")
    ax1.plot(t, sigma_out[1, :], label="sigma_y")
    ax1.plot(t, sigma_out[2, :], label="sigma_z")
    ax1.set_title("Изменение управляющего момента")
    ax1.legend()
    ax1.grid(True)

    ax2.plot(t, l_delta[:, 0], label="L[0]")
    ax2.plot(t, l_delta[:, 1], label="L[1]")
    ax2.plot(t, l_delta[:, 2], label="L[2]")
    ax2.plot(t, l_delta[:, 3], label="L[3]")
    ax2.set_title("Изменение кватерниона рассогласования")
    ax2.legend()
    ax2.grid(True)

    """
    ax1.plot(t, l_cur[:, 0], label="L[0]")
    ax1.plot(t, l_cur[:, 1], label="L[1]")
    ax1.plot(t, l_cur[:, 2], label="L[2]")
    ax1.plot(t, l_cur[:, 3], label="L[3]")
    ax1.set_title("Изменение кватерниона текущей ориентации")
    ax1.legend()
    ax1.grid(True)
    fig4, (ax1, ax2) = plt.subplots(2,1)
    
    """

    fig1.tight_layout()
    fig2.tight_layout()
    fig3.tight_layout()
    # fig4.tight_layout()
    plt.show()
