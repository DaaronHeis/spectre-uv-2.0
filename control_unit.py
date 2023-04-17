"""
    Модуль, вырабатывающий закон управления
"""

import numpy as np
import quaternion as qt
import random

pi = np.pi


def from_euler_to_quat(angles, order):
    """
        Перевод из углов Эйлера (или их аналогов) в кватернион

        Принимает
        ---------
            angles: list[float] или NumPy array углов
                gamma

                theta

                psi

            order: строка, определяющая порядок углов в формате XYZr

                r - оси вращаются вместе с телом

                s - оси неподвижны относительно тела

        Возвращает
        ---------
            q - кватернион
    """
    if angles is np.array:
        angles = angles.tolist()
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

            dL/dt =  0.5 * L * w

        Принимает:
        ---------
            L: quaternion - кватернион текущей ориентации
            w: list[float, float, float] -  текущее измерение угловой скорости
            w_prev: list[float, float, float] - предыдущее измерение угловой скорости
            dt: float - шаг интегрирования

        Возвращает
        ---------
            L_rotated: quaternion - кватернион ориентации после поворота
    """
    if w is np.array:
        w = w.tolist()
    if w_prev is np.array:
        w_prev = w_prev.tolist()
    avg_w = [0.0, 0.0, 0.0]
    avg_w[0] = 0.5 * (w[0] + w_prev[0])
    avg_w[1] = 0.5 * (w[1] + w_prev[1])
    avg_w[2] = 0.5 * (w[2] + w_prev[2])
    length = np.sqrt(avg_w[0]**2 + avg_w[1]**2 + avg_w[2]**2)
    theta = 0.5 * length * dt
    if length > 1.0e-12:
        w = np.cos(theta)
        s = np.sin(theta) / length
        q = qt.quaternion(w, avg_w[0]*s, avg_w[1]*s, avg_w[2]*s)
        if q.norm() > 1.01:
            q = q.normalized()
    else:
        q = qt.quaternion(1,0,0,0)

    L_rotated = L * q
    return L_rotated


class GIVUS:
    """
        Модель волоконно-оптического гироскопа Astrix-200
        Измеряет угловую скорость космического аппарата
        В отличие от астродатчика, не хранит измеренные значения скорости
    """
    def __init__(self, GIVUS_NOISE_KEY):

        self.name = 'Astrix-200'
        self.noise = 0.0001                             # погрешность измерения (град/sqrt(час))
        self.bias_stab = 0.0005                         # уход погрешности за час (град/час)
        self.scale_factor = 30 / 10000                  # нестабильность масштабного коэффициента
        self.key_noise = GIVUS_NOISE_KEY                # ключ учета погрешности измерений

    def add_angular_velocity_noise(self, w, t, dt):
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

        # добавление случайного дрейфа (angular random walk)
        step_size = self.bias_stab/3600/180*pi * dt
        sigma = (step_size**2) * t/dt
        w += random.choice([-1, 1]) * np.sqrt(sigma) + random.choice([-1, 1])*self.noise

        # добавление нестабильности масштабного коэффициента
        dw = self.scale_factor * w
        w += random.choice([-1, 1]) * dw

        return w

    def measure_velocity(self, w, t, dt):
        """ Измерение угловой скорости """
        w_meas = w.copy()
        if self.key_noise > 0:
            w_meas = self.add_angular_velocity_noise(w_meas, t, dt)
        
        return w_meas


class AstroSensor:
    """
        Матмодель астродатчика
        Измеряет (или каким-то другим образом получает) независимую от полученной блоком управления ориентацию КА,
        по которой, используя фильтр Чебышева 1-го рода (2 пор.), получает угловую скорость астрокоррекции
    """
    def __init__(self, L_meas, w_cor, KF, ERROR_KEY: bool):

        self.name = '348K'                              # имя астродатчика (из документации)
        self.L_meas = L_meas                            # текущий кватернион ориентации, измеренный астродатчиком
        self.L_cor = qt.quaternion(1,0,0,0)             # текущий кватернион коррекции
        self.w_cor = w_cor                              # текущая рассчитанная угловая скорость астрокоррекции
        self.delta = 12.0 / 60/60/2/3.1415926           # погрешность измерения астродатчика
        self.key_noise = ERROR_KEY                      # ключ учёта погрешностей измерения

        """
            Используется линейный фильтр Калмана
        """
        self.KF = KF
        self.k = 1                                     # коэффициент усиления

    def add_noise(self):
        """
            Добавление погрешностей измерения астродатчика

            Параметры
            ------------
            L: quaternion - измеренный кватернион ориентации
        """
        deltaL = qt.quaternion(1, random.choice([-1,1])*self.delta, random.choice([-1,1])*self.delta, random.choice([-1,1])*self.delta)
        self.L_meas = self.L_meas * deltaL

    def set_current_orientation(self, L, key_stab):
        """
            "Делает снимок неба", но на самом деле просто получает кватернион ориентации от handler-a 
        """
        self.L_meas = L.copy()
        if (self.key_noise is True) and (key_stab is False):
            self.add_noise()

    def filter(self, w):
        """
            Моделирование работы фильтра
            Разностное уравнение, составленное из передаточной функции
        """
        return self.KF.filter(self.L_meas, w)

    def set_correction(self, L_from_CU, w):
        """
            Определяет текущую ориентацию аппарата и
            соответствующую угловую скорость астрокоррекции

            Параметры
            ------------
            w: np.array 3х1 - текущая угловая скорость
            w_prev: np.array 3x1 - прошлая угловая скорость
            dt: float - шаг интегрироваия
        """

        """ Вычисление угловой скорости астрокоррекции """
        # определение кватерниона коррекции и углов расхождения
        self.L_cor = self.L_meas.inverse() * L_from_CU
        L_delta = qt.as_float_array(self.L_cor)
        fi = [0, 0, 0]
        fi[0] = 2 * L_delta[0] * L_delta[1]
        fi[1] = 2 * L_delta[0] * L_delta[2]
        fi[2] = 2 * L_delta[0] * L_delta[3]
        fi = np.array(fi)

        # получение угловой скорости астрокоррекции через фильтр
        self.w_cor = self.filter(w)

    def get_correction(self, L_from_CU, w):
        """
            Вычисляет угловую скорость астрокоррекции и возвращает её.
            Вызывается из модуля СУ
        """
        self.set_correction(L_from_CU, w)
        return self.w_cor

    def get_parameters(self):
        params = [self.L_meas, self.L_cor, self.w_cor, self.fi_prev]
        return params


class ControlUnit:
    """
        Матмодель расчёта управления
        Нужный управляемый момент рассчитывается с использованием кватернионов ориентации
    """
    def __init__(self, L_pr: qt.quaternion,
                       L_cur: qt.quaternion,
                       L_delta: qt.quaternion,
                       omega, gamma, omega_pr, w_bw, sigma_max,
                     ASTRO_CORR_KEY: bool, ARTIF_ERR_KEY: bool):
        """
            Кватернионы ориентации отсчитываются от ЭСК-2. То есть, высчитывается положение
            начальной точки поворота, конечной точки в ЭСК-2, и затем рассчитывается по этим известным
            углам кватернион поворота. Базовое положение (то, от которого отсчитываются повороты) - [1 0 0 0]
        """
        self.L_pr = L_pr                                       # кватернион программной ориентации
        self.L_cur = L_cur                                     # кватернион текущей ориентации
        self.L_delta = L_delta                                 # кватернион рассогласования

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

        self.K1_st = np.diag([0.025, 0.04375, 0.05875]) * 1
        self.K2_st = np.diag([1.755, 3.375, 4.365]) * 1

        self.threshold = 3 / 180 * np.pi                        # порог переключения

        self.omega_measured = omega                             # текущее измерение угловой скорости ГИВУСом
        self.omega = omega                                      # текущий вектор угловой скорости (со всеми шумами и коррекциями)
        self.gamma = gamma                                      # текущие вычисляемые углы поворота КА

        self.omega_pr = omega_pr                                # программный вектор угловой скорости
        self.omega_delta = np.array([0.0, 0.0, 0.0])            # разница между текущей и программной скоростью
        self.omega_cor = np.array([0.0, 0.0, 0.0])              # угловая скорость коррекции (от астродатчика)
        self.omega_prev = omega                                 # угловая скорость после коррекции на прошлом шаге
        self.sigma_prev = np.array([0.0, 0.0, 0.0])

        self.sigma_max = sigma_max
        self.w_bw = w_bw
        self.flag = False                                       # флаг перехода в режим стабилизации
        self.key_corr = ASTRO_CORR_KEY                          # ключ, подключающий астрокоррекцию
        self.key_orient_error = ARTIF_ERR_KEY                   # ключ, подключающий искусственную ощибку ориентации

    def set_gamma(self):
        """ Определение вектора углового положения КА """
        l_delta = qt.as_float_array(self.L_delta)
        self.gamma[0] = 2 * l_delta[0] * l_delta[1]
        self.gamma[1] = 2 * l_delta[0] * l_delta[2]
        self.gamma[2] = 2 * l_delta[0] * l_delta[3]
        self.gamma = np.array(self.gamma)

    def set_delta(self):
        """ Вычисление кватерниона рассогласования """
        self.L_delta = self.L_pr.inverse() * self.L_cur

    def get_correction(self, astrosensor):
        """ Получение угловой скорости коррекции от астродатчика """
        self.omega_cor = astrosensor.get_correction(self.L_cur, self.omega_measured)

    def get_control_moment(self, t):
        """ Расчёт и получение управляющего момента """
        self.set_delta()
        self.set_gamma()

        if (abs(self.gamma[0]) < self.threshold
            and abs(self.gamma[1]) < self.threshold
            and abs(self.gamma[2]) < self.threshold):

            # стабилизация
            if self.flag == False:
                self.flag = True
                print(t, self.flag)
            self.omega_delta = self.omega - self.omega_pr
            sigma = -(self.K1_st @ self.gamma + self.K2_st @ self.omega_delta)
            #sigma = sigma + np.sign(sigma)*0.0012
        else:

            # перенацеливание
            if self.flag == True:
                self.flag = False
                print(t, self.flag)
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

    def get_current_orientation(self):
        return self.L_cur

    def is_stabilisation(self):
        return self.flag

    def set_current_orientation(self, dt, astrosensor):
        """ Вычисление кватерниона текущей ориентации путем интегрирования уравнений движения """
        omega = self.omega
        if self.key_orient_error is True:
            self.add_error_in_orientation()
        if self.key_corr is True:
            self.get_correction(astrosensor)
            omega = self.omega + self.omega_cor
        
        omega_prev = [self.omega_prev[0], self.omega_prev[1], self.omega_prev[2]]
        omega = [omega[0], omega[1], omega[2]]
        self.L_cur = integrate_angular_velocity(self.L_cur, omega, omega_prev, dt)

    def set_current_velocity(self, w):
        """ Запоминание текущей угловой скорости для интегрирования при вычислении ориентации """
        self.omega_prev = self.omega
        self.omega = w

    def add_error_in_orientation(self):
        """
        Добавление искусственной ошибки при определении ориентации КА. Использовать для
        проверки работы астрокоррекции
        """
        k = 0.0001
        error_quat = qt.quaternion(1, k*random.choice([-1, 1]), k*random.choice([-1, 1]), k*random.choice([-1, 1]))
        error_quat = error_quat.normalized()
        self.L_cur = self.L_cur*error_quat

    def get_parameters(self):
        """
            Получение текущих параметров модуля:
                l_cur
                l_delta
                gamma
                omega
        """
        params = [self.L_cur, self.L_delta, self.gamma, self.omega]
        return params
