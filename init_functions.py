"""
    Функции инициализации различных элементов модели
"""

import numpy as np
import quaternion as qt
import flywheel_engine as dm
import control_unit as ctrl
import kalman_filter as kf


def init_time(t_span, dt):
    """ Инициализация временного отрезка моделирования """
    t_begin = t_span[0]
    t_end = t_span[1]
    h = dt
    t = [t_begin]
    t_curr = t_begin
    return t, t_curr, h, t_end


def init_flywheels(n, w_bw, sigma_max):
    """ Инициализация двигателей-маховиков """
    dm_all = []
    for j in range(n):
        dm_all.append(dm.Flywheel(0,0,0,0,w_bw))
    return dm_all


def init_target_orientation(angles_end, vel_end):
    """ Выставление требуемой ориентации через заданные углы """
    gamma_pr = angles_end.copy()
    l_k = ctrl.from_euler_to_quat(gamma_pr, 'YZXr')
    # l_pr = l_k.inverse() * l_0
    # l_pr = qt.quaternion(1, 0, 0, 0)
    l_pr = qt.quaternion(0.86472620667614, 0.256908589358277, 0.334920502035886, 0.272166900113631)
    print('l_pr = ', l_pr)
    return l_pr


def init_start_orientation(angles_0):
    """ Выставление начальной ориентации (сразу в виде кватерниона) """
    angles = angles_0.copy()
    l_0 = ctrl.from_euler_to_quat(angles, 'YZXr')
    print('l_0 = ', l_0)
    return l_0


def init_control_unit(l_0, l_pr, vel_0, omega_pr, w_bw, sigma_max,
                      CORR_KEY: bool, ARTIF_ERR_KEY: bool):
    """ Инициализация модуля блока управления """
    vel = vel_0.copy()
    l_cur = l_0.copy()
    l_delta = l_pr.inverse() * l_cur
    print('l_delta = ', l_delta)
    angles = np.array([2 * l_delta.w * l_delta.x, 2 * l_delta.w * l_delta.y, 2 * l_delta.w * l_delta.z])
    ctrl_unit = ctrl.ControlUnit(l_pr, l_cur, l_delta, vel, angles, omega_pr,
                                 w_bw, sigma_max, CORR_KEY, ARTIF_ERR_KEY)
    return angles, vel, ctrl_unit


def init_GIVUS(GIVUS_ERR_KEY: bool):
    """ Инициализация модуля ГИВУСа """
    givus = ctrl.GIVUS(GIVUS_ERR_KEY)
    return givus


def init_astrosensor(L_0, L_pr, dt, A_S_ERR_KEY: bool):
    """ Инициализация модуля астродатчика """
    KF = init_kalman(L_0, L_pr, dt)
    astrosensor = ctrl.AstroSensor(L_0, np.array([0,0,0]), KF, A_S_ERR_KEY)
    return astrosensor


def init_kalman(L_KA, L_pr, dt):
    """ Инициализация фильтра Калмана """
    # инициализация начального состояния
    l_delta = qt.as_float_array(L_pr.inverse() * L_KA)
    x = [0., 0., 0.]
    x[0] = 2 * l_delta[0] * l_delta[1]
    x[1] = 2 * l_delta[0] * l_delta[2]
    x[2] = 2 * l_delta[0] * l_delta[3]
    x = np.array(x)
    P = np.diag([1e-10, 1e-10, 1e-10])
    Q = np.diag([1e-10, 1e-10, 1e-10])
    delta = 12.0 / 60/60/2/3.1415926
    R = np.diag([delta, delta, delta])
    F = np.zeros_like(P)
    H = np.eye(3)
    B = np.eye(3)
    KF = kf.LinearKalmanFilter(x, P, Q, R, dt, L_pr, F, H, B)
    return KF
    