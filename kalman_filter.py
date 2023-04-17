"""
    Моделирование линейного фильтра Калмана
    Используется при астрокоррекции

    Алгоритм:

    0) инициализация

    1) прогноз (predict):
    
        - состояние с прошлого шага умножается на матрицу процесса F

        - уточняется достоверность расчта для поправки на неточность модели

    2) обновление (update):

        - производится измерение и вычисляется его достоверность

        - вычисляется разница между предсказанным состоянием и измерением

        - вычисляется масштабирующий фактор на основе достоверностей

        - вычисляется итоговое состояние

        - обновляется достоверность состояния на основе точности измерения
"""

import numpy as np
import quaternion as qt
from scipy.linalg import inv


class LinearKalmanFilter:

    def __init__(self, x, P, Q, R, dt, L_pr, F, H, B):

        self.x = x              # вектор состояния
        self.x_predicted = x    # предсказанный вектор состояния
        self.P = P              # матрица ковариации 
        self.P_predicted = P    # предсказанная матрица ковариации
        self.Q = Q              # шум процесса
        self.R = R              # шум измерений
        self.n = len(x)         # размерность вектора состояния
        self.dt = dt            # такт работы системы
        self.L_pr = L_pr        # кватернион программмной ориентации
        self.F = F              # матрица процесса
        self.H = H              # матрица измерения
        self.B = B              # матрица управления

    def predict(self, u):
        """
            Прогноз
            u - подаваемое управление
        """
        self.x_predicted = self.F @ self.x + self.B @ u
        self.P_predicted = self.F @ self.P @ self.F.T + self.Q

    def update(self, z):
        """
            Коррекция
            z - результат измерений
        """
        S = self.H @ self.P @ self.H.T + self.R
        K = self.P_predicted @ self.H.T @ inv(S)
        y = z - self.H @ self.x_predicted
        self.x = self.x_predicted + K @ y
        self.P = self.P_predicted - K @ self.H @ self.P_predicted

    def get_measurements(self, L):
        """
            Метод, преобразовывающий полученный от астродатчика кватернион рассогласования в пространство состояния
        """
        L = self.L_pr.inverse() * L
        L = qt.as_float_array(L)
        z = [0., 0., 0.]
        z[0] = 2 * L[0] * L[1]
        z[1] = 2 * L[0] * L[2]
        z[2] = 2 * L[0] * L[3]
        return np.array(z)

    def filter(self, L_AS, w):
        """
            Метод, реализующий одну итерацию алгоритма фильтра
        """
        # прогноз
        self.predict(w)
        # получение измерения
        z = self.get_measurements(L_AS)
        # коррекция
        self.update(z)

        # расчёт скорости рассогласования
        w_cor = self.x / self.dt
        return w_cor
