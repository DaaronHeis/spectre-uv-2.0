"""
    Модуль, вырабатывающий закон управления
"""

import numpy as np
import quaternion as qt
import random


def from_euler_to_quat(angles, order):
    """
        Перевод из углов Эйлера (или их аналогов) в кватернион

        Принимает:
            angles: массив float углов
                gamma
                theta
                psi
            order: строка, определяющая порядок углов в формате XYZr
                r - оси вращаются вместе с телом
                s - оси неподвижны относительно тела
    """
    q = qt.quaternion(0, 0, 0, 0)
    a = angles[0]/2
    b = angles[1]/2
    g = angles[2]/2
    if order == 'ZYXr':
        q.x = np.cos(b) * np.sin(g) * np.cos(a) - np.sin(b) * np.cos(g) * np.sin(a)
        q.y = np.cos(b) * np.sin(g) * np.sin(a) + np.sin(b) * np.cos(g) * np.cos(a)
        q.z = np.cos(b) * np.cos(g) * np.sin(a) - np.sin(b) * np.sin(g) * np.cos(a)
        q.w = np.cos(b) * np.cos(g) * np.cos(a) + np.sin(b) * np.sin(g) * np.sin(a)
    elif order == 'XYXr':
        q.x = np.cos(b) * np.cos(g) * np.sin(a) + np.cos(b) * np.sin(g) * np.cos(a)
        q.y = np.sin(b) * np.cos(g) * np.cos(a) + np.sin(b) * np.sin(g) * np.sin(a)
        q.z = np.sin(b) * np.cos(g) * np.sin(a) - np.sin(b) * np.sin(g) * np.cos(a)
        q.w = np.cos(b) * np.cos(g) * np.cos(a) - np.cos(b) * np.sin(g) * np.sin(a)
    elif order == 'YZXr':
        q.x = np.cos(-b) * np.sin(g) * np.cos(a) - np.sin(-b) * np.cos(g) * np.sin(a)
        q.z = -np.cos(-b) * np.sin(g) * np.sin(a) - np.sin(-b) * np.cos(g) * np.cos(a)
        q.y = np.cos(-b) * np.cos(g) * np.sin(a) - np.sin(-b) * np.sin(g) * np.cos(a)
        q.w = np.cos(-b) * np.cos(g) * np.cos(a) + np.sin(-b) * np.sin(g) * np.sin(a)
    elif order == 'XZXr':
        q.x = np.cos(-b) * np.cos(g) * np.sin(a) + np.cos(-b) * np.sin(g) * np.cos(a)
        q.z = -np.sin(-b) * np.cos(g) * np.cos(a) - np.sin(-b) * np.sin(g) * np.sin(a)
        q.y = np.sin(-b) * np.cos(g) * np.sin(a) - np.sin(-b) * np.sin(g) * np.cos(a)
        q.w = np.cos(-b) * np.cos(g) * np.cos(a) - np.cos(-b) * np.sin(g) * np.sin(a)
    elif order == 'XZYr':
        q.y = np.cos(b) * np.sin(g) * np.cos(a) - np.sin(b) * np.cos(g) * np.sin(a)
        q.z = np.cos(b) * np.sin(g) * np.sin(a) + np.sin(b) * np.cos(g) * np.cos(a)
        q.x = np.cos(b) * np.cos(g) * np.sin(a) - np.sin(b) * np.sin(g) * np.cos(a)
        q.w = np.cos(b) * np.cos(g) * np.cos(a) + np.sin(b) * np.sin(g) * np.sin(a)
    elif order == 'YZYr':
        q.y = np.cos(b) * np.cos(g) * np.sin(a) + np.cos(b) * np.sin(g) * np.cos(a)
        q.z = np.sin(b) * np.cos(g) * np.cos(a) + np.sin(b) * np.sin(g) * np.sin(a)
        q.x = np.sin(b) * np.cos(g) * np.sin(a) - np.sin(b) * np.sin(g) * np.cos(a)
        q.w = np.cos(b) * np.cos(g) * np.cos(a) - np.cos(b) * np.sin(g) * np.sin(a)
    elif order == 'ZXYr':
        q.y = np.cos(-b) * np.sin(g) * np.cos(a) - np.sin(-b) * np.cos(g) * np.sin(a)
        q.x = -np.cos(-b) * np.sin(g) * np.sin(a) - np.sin(-b) * np.cos(g) * np.cos(a)
        q.z = np.cos(-b) * np.cos(g) * np.sin(a) - np.sin(-b) * np.sin(g) * np.cos(a)
        q.w = np.cos(-b) * np.cos(g) * np.cos(a) + np.sin(-b) * np.sin(g) * np.sin(a)
    elif order == 'YXYr':
        q.y = np.cos(-b) * np.cos(g) * np.sin(a) + np.cos(-b) * np.sin(g) * np.cos(a)
        q.x = -np.sin(-b) * np.cos(g) * np.cos(a) - np.sin(-b) * np.sin(g) * np.sin(a)
        q.z = np.sin(-b) * np.cos(g) * np.sin(a) - np.sin(-b) * np.sin(g) * np.cos(a)
        q.w = np.cos(-b) * np.cos(g) * np.cos(a) - np.cos(-b) * np.sin(g) * np.sin(a)
    elif order == 'YXZr':
        q.z = np.cos(b) * np.sin(g) * np.cos(a) - np.sin(b) * np.cos(g) * np.sin(a)
        q.x = np.cos(b) * np.sin(g) * np.sin(a) + np.sin(b) * np.cos(g) * np.cos(a)
        q.y = np.cos(b) * np.cos(g) * np.sin(a) - np.sin(b) * np.sin(g) * np.cos(a)
        q.w = np.cos(b) * np.cos(g) * np.cos(a) + np.sin(b) * np.sin(g) * np.sin(a)
    elif order == 'ZXZr':
        q.z = np.cos(b) * np.cos(g) * np.sin(a) + np.cos(b) * np.sin(g) * np.cos(a)
        q.x = np.sin(b) * np.cos(g) * np.cos(a) + np.sin(b) * np.sin(g) * np.sin(a)
        q.y = np.sin(b) * np.cos(g) * np.sin(a) - np.sin(b) * np.sin(g) * np.cos(a)
        q.w = np.cos(b) * np.cos(g) * np.cos(a) - np.cos(b) * np.sin(g) * np.sin(a)
    elif order == 'XYZr':
        q.z = np.cos(-b) * np.sin(g) * np.cos(a) - np.sin(-b) * np.cos(g) * np.sin(a)
        q.y = -np.cos(-b) * np.sin(g) * np.sin(a) - np.sin(-b) * np.cos(g) * np.cos(a)
        q.x = np.cos(-b) * np.cos(g) * np.sin(a) - np.sin(-b) * np.sin(g) * np.cos(a)
        q.w = np.cos(-b) * np.cos(g) * np.cos(a) + np.sin(-b) * np.sin(g) * np.sin(a)
    elif order == 'ZYZr':
        q.z = np.cos(-b) * np.cos(g) * np.sin(a) + np.cos(-b) * np.sin(g) * np.cos(a)
        q.y = -np.sin(-b) * np.cos(g) * np.cos(a) - np.sin(-b) * np.sin(g) * np.sin(a)
        q.x = np.sin(-b) * np.cos(g) * np.sin(a) - np.sin(-b) * np.sin(g) * np.cos(a)
        q.w = np.cos(-b) * np.cos(g) * np.cos(a) - np.cos(-b) * np.sin(g) * np.sin(a)
    elif order == 'ZYXs':
        q.z = np.cos(-b) * np.sin(a) * np.cos(g) - np.sin(-b) * np.cos(a) * np.sin(g)
        q.y = -np.cos(-b) * np.sin(a) * np.sin(g) - np.sin(-b) * np.cos(a) * np.cos(g)
        q.x = np.cos(-b) * np.cos(a) * np.sin(g) - np.sin(-b) * np.sin(a) * np.cos(g)
        q.w = np.cos(-b) * np.cos(a) * np.cos(g) + np.sin(-b) * np.sin(a) * np.sin(g)
    elif order == 'XYXs':
        q.x = np.cos(b) * np.cos(a) * np.sin(g) + np.cos(b) * np.sin(a) * np.cos(g)
        q.y = np.sin(b) * np.cos(a) * np.cos(g) + np.sin(b) * np.sin(a) * np.sin(g)
        q.z = np.sin(b) * np.cos(a) * np.sin(g) - np.sin(b) * np.sin(a) * np.cos(g)
        q.w = np.cos(b) * np.cos(a) * np.cos(g) - np.cos(b) * np.sin(a) * np.sin(g)
    elif order == 'YZXs':
        q.y = np.cos(b) * np.sin(a) * np.cos(g) - np.sin(b) * np.cos(a) * np.sin(g)
        q.z = np.cos(b) * np.sin(a) * np.sin(g) + np.sin(b) * np.cos(a) * np.cos(g)
        q.x = np.cos(b) * np.cos(a) * np.sin(g) - np.sin(b) * np.sin(a) * np.cos(g)
        q.w = np.cos(b) * np.cos(a) * np.cos(g) + np.sin(b) * np.sin(a) * np.sin(g)
    elif order == 'XZXs':
        q.x = np.cos(-b) * np.cos(a) * np.sin(g) + np.cos(-b) * np.sin(a) * np.cos(g)
        q.z = -np.sin(-b) * np.cos(a) * np.cos(g) - np.sin(-b) * np.sin(a) * np.sin(g)
        q.y = np.sin(-b) * np.cos(a) * np.sin(g) - np.sin(-b) * np.sin(a) * np.cos(g)
        q.w = np.cos(-b) * np.cos(a) * np.cos(g) - np.cos(-b) * np.sin(a) * np.sin(g)
    elif order == 'XZYs':
        q.x = np.cos(-b) * np.sin(a) * np.cos(g) - np.sin(-b) * np.cos(a) * np.sin(g)
        q.z = -np.cos(-b) * np.sin(a) * np.sin(g) - np.sin(-b) * np.cos(a) * np.cos(g)
        q.y = np.cos(-b) * np.cos(a) * np.sin(g) - np.sin(-b) * np.sin(a) * np.cos(g)
        q.w = np.cos(-b) * np.cos(a) * np.cos(g) + np.sin(-b) * np.sin(a) * np.sin(g)
    elif order == 'YZYs':
        q.y = np.cos(b) * np.cos(a) * np.sin(g) + np.cos(b) * np.sin(a) * np.cos(g)
        q.z = np.sin(b) * np.cos(a) * np.cos(g) + np.sin(b) * np.sin(a) * np.sin(g)
        q.x = np.sin(b) * np.cos(a) * np.sin(g) - np.sin(b) * np.sin(a) * np.cos(g)
        q.w = np.cos(b) * np.cos(a) * np.cos(g) - np.cos(b) * np.sin(a) * np.sin(g)
    elif order == 'ZXYs':
        q.z = np.cos(b) * np.sin(a) * np.cos(g) - np.sin(b) * np.cos(a) * np.sin(g)
        q.x = np.cos(b) * np.sin(a) * np.sin(g) + np.sin(b) * np.cos(a) * np.cos(g)
        q.y = np.cos(b) * np.cos(a) * np.sin(g) - np.sin(b) * np.sin(a) * np.cos(g)
        q.w = np.cos(b) * np.cos(a) * np.cos(g) + np.sin(b) * np.sin(a) * np.sin(g)
    elif order == 'YXYs':
        q.y = np.cos(-b) * np.cos(a) * np.sin(g) + np.cos(-b) * np.sin(a) * np.cos(g)
        q.x = -np.sin(-b) * np.cos(a) * np.cos(g) - np.sin(-b) * np.sin(a) * np.sin(g)
        q.z = np.sin(-b) * np.cos(a) * np.sin(g) - np.sin(-b) * np.sin(a) * np.cos(g)
        q.w = np.cos(-b) * np.cos(a) * np.cos(g) - np.cos(-b) * np.sin(a) * np.sin(g)
    elif order == 'YXZs':
        q.y = np.cos(-b) * np.sin(a) * np.cos(g) - np.sin(-b) * np.cos(a) * np.sin(g)
        q.x = -np.cos(-b) * np.sin(a) * np.sin(g) - np.sin(-b) * np.cos(a) * np.cos(g)
        q.z = np.cos(-b) * np.cos(a) * np.sin(g) - np.sin(-b) * np.sin(a) * np.cos(g)
        q.w = np.cos(-b) * np.cos(a) * np.cos(g) + np.sin(-b) * np.sin(a) * np.sin(g)
    elif order == 'ZXZs':
        q.z = np.cos(b) * np.cos(a) * np.sin(g) + np.cos(b) * np.sin(a) * np.cos(g)
        q.x = np.sin(b) * np.cos(a) * np.cos(g) + np.sin(b) * np.sin(a) * np.sin(g)
        q.y = np.sin(b) * np.cos(a) * np.sin(g) - np.sin(b) * np.sin(a) * np.cos(g)
        q.w = np.cos(b) * np.cos(a) * np.cos(g) - np.cos(b) * np.sin(a) * np.sin(g)
    elif order == 'XYZs':
        q.x = np.cos(b) * np.sin(a) * np.cos(g) - np.sin(b) * np.cos(a) * np.sin(g)
        q.y = np.cos(b) * np.sin(a) * np.sin(g) + np.sin(b) * np.cos(a) * np.cos(g)
        q.z = np.cos(b) * np.cos(a) * np.sin(g) - np.sin(b) * np.sin(a) * np.cos(g)
        q.w = np.cos(b) * np.cos(a) * np.cos(g) + np.sin(b) * np.sin(a) * np.sin(g)
    elif order == 'ZYZs':
        q.z = np.cos(-b) * np.cos(a) * np.sin(g) + np.cos(-b) * np.sin(a) * np.cos(g)
        q.y = -np.sin(-b) * np.cos(a) * np.cos(g) - np.sin(-b) * np.sin(a) * np.sin(g)
        q.x = np.sin(-b) * np.cos(a) * np.sin(g) - np.sin(-b) * np.sin(a) * np.cos(g)
        q.w = np.cos(-b) * np.cos(a) * np.cos(g) - np.cos(-b) * np.sin(a) * np.sin(g)
    else:
        print('{order} is not supported, choose another')
    return q


def integrate_angular_velocity(L: qt.quaternion, w, w_prev, dt):
    """
        Интегрирование уравнения

            dL
            --  =  0.5 * L * w
            dt

            L: quaternion - кватернион текущей ориентации
            w: list[float, float, float] -  текущее измерение угловой скорости
            w_prev: list[float, float, float] - предыдущее измерение угловой скорости
            dt: float - шаг интегрирования

        Возвращает:
            L_rotated: quaternion - кватернион ориентации после поворота
    """
    avg_w = [0.0, 0.0, 0.0]
    avg_w[0] = 0.5 * (w[0] + w_prev[0])
    avg_w[1] = 0.5 * (w[1] + w_prev[1])
    avg_w[2] = 0.5 * (w[2] + w_prev[2])
    len = np.sqrt(avg_w[0]**2 + avg_w[1]**2 + avg_w[2]**2)
    theta = 0.5 * len * dt
    if len > 1.0e-12:
        w = np.cos(theta)
        s = np.sin(theta) / len
        q = qt.quaternion(w, avg_w[0]*s, avg_w[1]*s, avg_w[2]*s)
        q = q.normalized()
    else:
        q = qt.quaternion(1,0,0,0)

    L_rotated = L * q
    return L_rotated


def add_angular_velocity_noise(w: np.array, t, dt):
    """
        Внесение шума в измерение угловой скорости от ГИВУСов по двум параметрам:
            случайный дрейф
            нестабильность масштабного коэффициента

        Параметры
        ----------
            w: np.float 3x1 - "чистая" текущая угловая скорость
            t: текущее время моделирования
            dt: шаг интегрирования

        Возвращает
        ----------
            w_noise: np.float 3x1 - "грязная" угловая скорость
    """
    w_noise = w.copy()

    # добавление случайного дрейфа
    step_size = 0.0005 / 3600 /2/3.1415926 * dt
    sigma = (step_size**2) * t/dt
    w_noise += random.choice([-1, 1]) * np.sqrt(sigma) + random.choice([-1, 1])*0.0001

    # добавление нестабильности масштабного коэффициента
    eps = 0.003 / 100
    dw = eps * w
    w_noise += random.choice([-1, 1]) * dw

    return w_noise


class AstroSensor:
    """
        Матмодель астродатчика
        Измеряет (или каким-то другим образом получает) независимую от полученной интегрированием ориентацию КА,
        по которой, используя фильтр Чебышева 2-го рода (2 пор.), получает угловую скорость астрокоррекции
    """
    def __init__(self, L_meas, w_cor, ERROR_KEY):

        self.name = '348K'                              # имя астродатчика (из документации)
        self.L_meas = L_meas                            # текущий кватернион ориентации, измеренный астродатчиком
        self.L_cor = L_meas                             # текущий кватернион коррекции
        self.w_cor = w_cor                              # текущая рассчитанная угловая скорость астрокоррекции
#       self.K = np.diag([0.1, 0.1, 0.1])               # матрица фильтра
        self.delta = 12.0 / 60/60/2/3.1415926           # погрешность измерения астродатчика
        self.ERROR_KEY = ERROR_KEY                      # флаг учёта погрешностей измерения
        self.fi_prev = None
        self.w_cor_prev = np.array([0.0, 0.0, 0.0])

    def add_noise(self):
        """
            Добавление погрешностей измерения астродатчика

            Параметры
            ------------
            L: quaternion - измеренный кватернион ориентации
        """
        self.L_meas *= qt.quaternion(1, self.delta, self.delta, self.delta)

    def set_current_orientation(self, L_old, w, w_prev, dt):
        """
            Так как звездное небо не моделируется, то измерения взять неоткуда, поэтому текущая ориентация
            КА получается интегрированием незашумленной угловой скорости таким же образом, что и в СУ

            Параметры
            ------------
            w: np.array 3х1 - текущая угловая скорость
            w_prev: np.array 3x1 - прошлая угловая скорость
            dt: float - шаг интегрироваия
        """
        self.L_meas = integrate_angular_velocity(L_old, w, w_prev, dt)

    def set_correction(self, L_old, w, w_prev, dt):
        """
            Определяет текущую ориентацию аппарата и
            соответствующую угловую скорость астрокоррекции

            Параметры
            ------------
            w: np.array 3х1 - текущая угловая скорость
            w_prev: np.array 3x1 - прошлая угловая скорость
            dt: float - шаг интегрироваия
        """
        """ Получение текущей ориентации """
        self.set_current_orientation(L_old, w, w_prev, dt)
        L_cur = self.L_meas
        if self.ERROR_KEY:
            self.add_noise()

        """ Вычисление угловой скорости астрокоррекции """
        self.L_cor = self.L_meas.inverse() * L_cur
        L_delta = qt.as_float_array(self.L_cor)
        fi = [0, 0, 0]
        fi[0] = 2 * L_delta[0] * L_delta[1]
        fi[1] = 2 * L_delta[0] * L_delta[2]
        fi[2] = 2 * L_delta[0] * L_delta[3]
        fi = np.array(fi)

        if self.fi_prev is None:
            self.fi_prev = fi
        else:
            # self.w_cor = -self.K @ fi
            self.w_cor = 0.6628*fi + 0.6628*self.fi_prev - 0.3255*self.w_cor_prev
            self.w_cor_prev = self.w_cor
            self.fi_prev = fi

    def get_correction(self, L_old, w, w_prev, dt):
        """
            Вычисляет угловую скорость астрокоррекции и возвращает её.
            Вызывается из модуля СУ
        """
        self.set_correction(L_old, w, w_prev, dt)
        return self.w_cor


class ControlUnit:
    """
        Матмодель расчёта управления
        Нужный управляемый момент рассчитывается с использованием кватернионов ориентации
    """
    def __init__(self, l_pr: qt.quaternion,
                       l_cur: qt.quaternion,
                       l_delta: qt.quaternion,
                       omega, gamma, omega_pr, I, w_bw, sigma_max, t, dt, CORR_KEY: bool, GIVUS_ERR_KEY: bool):
        """
            Кватернионы ориентации отсчитываются от ЭСК-2. То есть, высчитывается положение
            начальной точки поворота, конечной точки в ЭСК-2, и затем рассчитывается по этим известным
            углам кватернион поворота. Базовое положение (то, от которого отсчитываются повороты) - [1 0 0 0]
        """
        self.l_pr = l_pr                                       # кватернион программной ориентации
        self.l_cur = l_cur                                     # кватернион текущей ориентации
        self.l_delta = l_delta                                 # кватернион рассогласования

        """
            Диагональные матрицы для расчёта управляющего момента
            Являются управляющими коэффициентами, подбираются вручную
            
            --------------------------------------------------------
                            Основной вариант
                self.K1_pt = np.diag([0.03, 0.1175, 0.1175]) * 0.6
                self.K2_pt = np.diag([1.35, 4.55, 4.85]) * 1.1

                self.K1_st = np.diag([0.05, 0.0875, 0.1175]) * 0.5
                self.K2_st = np.diag([1.95, 3.75, 4.85]) * 0.9
            --------------------------------------------------------
            
        self.K1_pt = np.diag([0.1075, 0.0875, 0.1075]) * 0.5
        self.K2_pt = np.diag([4.85, 3.85, 4.85]) * 1.3

        self.K1_st = np.diag([0.0875, 0.0875, 0.1075]) * 0.5
        self.K2_st = np.diag([4.85, 3.85, 4.85]) * 1.1
        
        
        self.K1_pt = np.diag([0.1055, 0.1175, 0.1175]) * 0.6
        self.K2_pt = np.diag([2.985, 4.55, 4.85]) * 1.1

        self.K1_st = np.diag([0.1275, 0.0875, 0.1175]) * 0.5
        self.K2_st = np.diag([3.55, 3.85, 4.85]) * 0.9
        """

        # коэффициенты управления для перенацеливания
        # self.K1_pt = np.diag([0.0675, 0.0675, 0.0775]) * 0.6
        # self.K2_pt = np.diag([4.85, 6.85, 6.55]) * 0.8
        self.K1_pt = np.diag([0.03, 0.1175, 0.1175]) * 0.6
        self.K2_pt = np.diag([1.35, 4.55, 4.85]) * 1.1

        self.K1_st = np.diag([0.05, 0.0875, 0.1175]) * 0.5
        self.K2_st = np.diag([1.75, 3.85, 4.65]) * 0.9

        self.threshold = 3 / 180 * np.pi                        # порог переключения

        self.omega = omega                                      # текущий вектор угловой скорости
        self.gamma = gamma                                      # текущие вычисляемые углы поворота КА
        # угловая скорость берётся из показаний ГИВУСа (т.е. из основной матмодели)

        self.omega_pr = omega_pr                                # программный вектор угловой скорости
        self.omega_delta = np.array([0.0, 0.0, 0.0])            # разница между текущей и программной скоростью
        self.omega_cor = np.array([0.0, 0.0, 0.0])              # угловая скорость коррекции (от астродатчика)
        self.omega_prev = omega
        self.sigma_prev = np.array([0.0, 0.0, 0.0])

        self.sigma_max = sigma_max
        self.w_bw = w_bw
        self.flag = False                                       # флаг перехода в режим стабилизации
        self.dt = dt
        self.t = t
        self.corr_key = CORR_KEY                                # ключ, подключающий коррекцию по астродатчику
        self.vel_key = GIVUS_ERR_KEY                            # ключ, подключающий погрешности измерения ГИВУСа

    def set_gamma(self):
        """ Определение вектора углового положения КА """
        l_delta = qt.as_float_array(self.l_delta)
        self.gamma[0] = 2 * l_delta[0] * l_delta[1]
        self.gamma[1] = 2 * l_delta[0] * l_delta[2]
        self.gamma[2] = 2 * l_delta[0] * l_delta[3]
        self.gamma = np.array(self.gamma)

    def set_delta(self):
        """ Вычисление кватерниона рассогласования """
        self.l_delta = self.l_pr.inverse() * self.l_cur

    def get_correction(self, astrosensor):
        """ Получение угловой скорости коррекции от астродатчика """
        self.omega_cor = astrosensor.get_correction(self.l_cur, self.omega, self.omega_prev, self.dt)

    def get_control_moment(self):
        """ Получение управляющего момента """
        if (abs(self.gamma[0]) < self.threshold
            and abs(self.gamma[1]) < self.threshold
            and abs(self.gamma[2]) < self.threshold):

            # стабилизация
            if self.flag == False:
                self.flag = True
                print(self.t, self.flag)
            self.omega_delta = self.omega - self.omega_pr
            sigma = -(self.K1_st @ self.gamma + self.K2_st @ self.omega_delta)
            #sigma = sigma + np.sign(sigma)*0.0012
        else:

            # перенацеливание
            if self.flag == True:
                self.flag = False
                print(self.t, self.flag)
            # self.omega_delta = self.omega + self.K1_pt @ self.gamma
            # sigma = - self.K2_pt @ self.omega_delta
            sigma = -self.K1_pt @ self.gamma - self.K2_pt @ self.omega

        # насыщение управления
        for i in range(len(sigma)):
            if abs(sigma[i]) > self.sigma_max:
                sigma[i] = self.sigma_max * np.sign(sigma[i])
        for i in range(len(sigma)):
            delta_sigma = sigma[i] - self.sigma_prev[i]
            if abs(delta_sigma) > self.w_bw * 0.2:
                sigma[i] = self.sigma_prev[i] + np.sign(delta_sigma)*self.w_bw*0.2
        self.sigma_prev = sigma
        return sigma

    def set_current_orientation(self):
        """ Вычисление кватерниона текущей ориентации путем интегрирования уравнений движения """
        omega_prev = [self.omega_prev[0], self.omega_prev[1], self.omega_prev[2]]
        omega = [self.omega[0], self.omega[1], self.omega[2]]
        self.l_cur = integrate_angular_velocity(self.l_cur, omega, omega_prev, self.dt)

    def set_velocity(self, omega):
        """ Запоминание угловой скорости для интегрирования при вычислении ориентации """
        self.omega_prev = self.omega
        self.omega = omega
        if self.vel_key:
            self.omega = add_angular_velocity_noise(self.omega, self.t, self.dt)
        if self.corr_key:
            self.omega += self.omega_cor

    def update(self, vel, astrosensor):
        """
         Один шаг интегрирования для управляющего модуля
         Расчитывает управляющий момент и обновляет свое состояние

         Параметры
         ----------
         vel: np.array[1, 3] - текущая угловая скорость

         Возвращает
         ----------
         sigma: np.array[1, 3] - требуемые управляющие моменты по каждой оси
        """

        self.t += self.dt
        if self.corr_key:
            self.get_correction(astrosensor)
        self.set_velocity(vel)
        self.set_current_orientation()
        self.set_delta()
        self.set_gamma()
        sigma = self.get_control_moment()
        return sigma

    def get_parameters(self):
        """
            Получение текущих параметров модуля:
                l_cur
                l_delta
                gamma
                omega
        """
        return self.l_cur, self.l_delta, self.gamma, self.omega

