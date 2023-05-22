"""
    Моделирование фильтра Калмана (UKF)
    Используется при астрокоррекции

    Алгоритм:

    0) инициализация

    1) предсказание (predict)

        - по модели процесса (f) рассчитывается состояние на следующем шаге

        - уточняется достоверность рассчёта для поправки на неточность модели

    2) обновление (update)

        - производится измерение и высчитывается его достоверность

        - высчитывается разница между рассчитыным (предсказанным) состоянием и измерением

        - вычисляется масштабирующий фактор на основе достоверностей предсказанного
        состояния и полученного измерения

        - рассчитывается итоговое состояние

        - обновляется достоверность состояния на основе точности измерения
"""


import numpy as np
import quaternion as qt


class LinearKalmanFilter:
    """
        Раньше был реализован нелинейный вариант, после совещания с научником
        было решено реализовать линейную модель
    """
    def __init__(self, x, P, Q, R, dt, L_pr):

        self.x = x              # вектор состояния
        self.x_prev = x
        self.x_predicted = x    # предсказанный вектор состояния
        self.P = P              # матрица ковариации 
        self.P_predicted = P    # предсказанная матрица ковариации
        self.Q = Q              # шум процесса
        self.R = R              # шум измерений
        self.n = x.shape[0]     # размерность вектора состояния
        self.dt = dt            # такт работы системы
        self.L_pr = L_pr        # кватернион программмной ориентации
        self.k = 1

        self.F = np.eye(3)          # матрица F процесса
        self.H = np.eye(3)          # матрица H измерений
        self.B = np.eye(3) * dt     # матрица B управления

    def predict(self, u):
        """
            Прогноз

            Обновляет предсказанные состояние и ковариацию
        """
        self.x_predicted = self.F @ self.x + self.B @ u
        self.P_predicted = self.F @ self.P @ self.F.T + self.Q

    def update(self, z):
        """
            Коррекция
        """
        y = z - self.H @ self.x_predicted
        S = self.H @ self.P_predicted @ self.H.T + self.R
        K = self.P_predicted @ self.H.T @ np.linalg.inv(S)

        self.x = self.x_predicted + K @ y
        self.P = (np.eye(3) - K @ self.H) @ self.P_predicted

    def get_measurements(self, L_AS, L_CU):
        """
            Метод, преобразовывающий полученный от астродатчика кватернион рассогласования в пространство состояния
        """
        #l = L.inverse() * L_CU
        l = self.L_pr.inverse() * L_AS
        l = qt.as_float_array(l)
        z = [0., 0., 0.]
        z[0] = 2 * l[0] * l[1]
        z[1] = 2 * l[0] * l[2]
        z[2] = 2 * l[0] * l[3]
        return np.array(z)

    def filter(self, L_AS, L_CU, w):
        """
            Метод, реализующий одну итерацию алгоритма фильтра

            Принимает
            ---------
                parameters: {'I': I, 'dH': dH, 'H': H, 'M': M, 'g': gamma} - параметры интегрирования
                L_AS: quaternion - кватернион коррекции, полученный от астродатчика
                w: list[float] - угловая скорость по каждой оси

            Возвращает
            ----------
                w_cor: list[float] - угловая скорость коррекции
        """
        # прогноз
        self.x_prev = self.x
        self.predict(w)
        # получение измерения
        z = self.get_measurements(L_AS, L_CU)
        # коррекция
        self.update(z)

        # расчёт скорости рассогласования
        l = qt.as_float_array(L_CU)
        a = 2 * l[0] * l[1]
        b = 2 * l[0] * l[2]
        c = 2 * l[0] * l[3]
        x_CU = np.array([a, b, c])
        w_cor = (self.x - x_CU) / self.dt - w
        return w_cor * self.k

    def change_R(self, R):
        self.R = R
