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
            q: quaternion - кватернион поворота
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
    step_size = 0.0005 / 3600 * dt
    sigma = (step_size**2) * t/dt
    w_noise += random.choice([-1, 1]) * np.sqrt(sigma)

    # добавление нестабильности масштабного коэффициента
    eps = 0.02 / 100
    dw = eps * w
    w_noise += random.choice([-1, 1]) * dw

    return w_noise


class AstroSensor:
    """
        Матмодель астродатчика
        Пока необходимо лишь рассчитывать угловую скорость коррекции
    """
    def __init__(self, omega_cor):

        self.name = '348K'
        self.omega_cor = omega_cor

    def set_correction(self, omega, omega_noise):
        """
            Пока принимает угловую скорость коррекции снаружи, а не рассчитывает сам, потому
            что пока она равняется угловой сорости, добавленной из-за шумов

            Параметры
            ------------
            omega: np.array 3х1 - чистая угловая скорость
            omega_noise: np.array 3x1 - угловая скорость с шумами
        """
        self.omega_cor = -(omega_noise - omega)

    def get_correction(self):

        return self.omega_cor


class ControlUnit:
    """
        Матмодель расчёта управления
        Нужный управляемый момент рассчитывается с использованием кватернионов ориентации
    """
    def __init__(self, l_pr: qt.quaternion,
                       l_cur: qt.quaternion,
                       l_delta: qt.quaternion,
                       omega, gamma, omega_pr, I, w_bw, sigma_max, t, dt, corr_key: bool):
        """
            Кватернионы ориентации отсчитываются от ЭСК-2. То есть, высчитывается положение
            начальной точки поворота, конечной точки в ЭСК-2, и затем рассчитывается по этим известным
            углам кватернион поворота. Базовое положение (то, от которого отсчитываются повороты) - [1 0 0 0]
        """
        self.l_pr = l_pr                                       # кватернион программной ориентации
        self.l_cur = l_cur                                     # кватернион текущей ориентации
        self.l_delta: l_delta                                  # кватернион рассогласования

        """
            Диагональные матрицы для расчёта управляющего момента
            Являются управляющими коэффициентами, подбираются вручную
            
            --------------------------------------------------------
                            Основной вариант
                self.K1 = np.diag([0.0475, 0.0875, 0.1075]) * 2
                self.K2 = np.diag([4.85, 4.85, 4.85]) * 3
            --------------------------------------------------------
            
            self.K1 = np.diag([0.0475, 0.0875, 0.1075]) * 3
            self.K2 = np.diag([4.85, 4.85, 4.85]) * 4
            
            Для угла в 30 градусов возникают колебания, но потом затухают
            self.K1 = np.diag([0.0475, 0.0675, 0.0675]) * 3
            self.K2 = np.diag([4.25, 2.15, 3.85]) * 4
        """

        self.K1 = np.diag([0.0675, 0.0875, 0.1075]) * 2
        self.K2 = np.diag([4.85, 4.85, 4.85]) * 3

        self.omega = omega                                      # текущий вектор угловой скорости
        self.gamma = gamma                                      # текущие вычисляемые углы поворота КА
        # угловая скорость берётся из показаний ГИВУСа (т.е. из основной матмодели)

        self.omega_pr = omega_pr                                # программный вектор угловой скорости
        self.omega_delta = np.array([0.0, 0.0, 0.0])            # разница между текущей и программной скоростью
        self.omega_cor = np.array([0.0, 0.0, 0.0])              # угловая скорость коррекции (от астродатчика)
        self.omega_prev = omega

        self.sigma_max = sigma_max
        self.dt = dt
        self.t = t
        self.corr_key = corr_key                                # ключ, подключающий коррекцию по астродатчику

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
        self.omega_delta = np.array(self.omega) - np.array(self.omega_pr)

    def get_correction(self, astrosensor):
        """ Получение угловой скорости коррекции от астродатчика """
        self.omega_cor = astrosensor.get_correction()

    def get_control_moment(self):
        """ Получение управляющего момента """
        sigma = -(self.K1 @ self.gamma + self.K2 @ self.omega_delta)
        for i in range(len(sigma)):
            if abs(sigma[i]) > self.sigma_max:
                sigma[i] = self.sigma_max * np.sign(sigma[i])
        return sigma

    def set_current_orientation(self, t0, t1):
        """ Вычисление кватерниона текущей ориентации путем интегрирования уравнений движения """
        omega_prev = [self.omega_prev[0], self.omega_prev[1], self.omega_prev[2]]
        omega = [self.omega[0], self.omega[1], self.omega[2]]
        self.l_cur = integrate_angular_velocity(self.l_cur, omega, omega_prev, self.dt)

    def set_velocity(self, omega):
        """ Запоминание угловой скорости для интегрирования при вычислении ориентации """
        self.omega_prev = self.omega
        self.omega = omega + self.omega_cor

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
        self.set_current_orientation(self.t-self.dt, self.t)
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

