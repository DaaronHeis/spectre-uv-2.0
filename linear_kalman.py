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

# TODO: переписать как линейный (функционал уже как у линейного, осталось только оформить красиво)

import numpy as np
import quaternion as qt
from scipy.linalg import cholesky 
from control_unit import from_euler_to_quat
from control_unit import integrate_angular_velocity


class UnscentedKalmanFilter:
    """
        Раньше был реализован нелинейный вариант, после совещания с научником
        было решено реализовать линейную модель
    """
    def __init__(self, x, P, Q, R, dt, L_pr, alpha, beta, kappa):

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

        self.alpha = alpha      # параметр генерации сигма-точек
        self.beta = beta        # параметр генерации сигма-точек
        self.kappa = kappa      # параметр генерации сигма-точек

        lambda_ = self.alpha ** 2 * (self.n + self.kappa) - self.n
        self.Wm = np.full(2*self.n+1, 1. / (2*(self.n + lambda_)))
        self.Wc = np.full(2*self.n+1, 1. / (2*(self.n + lambda_)))
        self.sigmas_f = np.zeros((2*self.n+1, self.n))
        self.sigmas_h = np.zeros((2*self.n+1, self.n))
        self.num_sigmas = 2*self.n+1

    def compute_sigma_points(self, mu):
        """
            Метод вычисления сигма-точек для UKF (Van der Merwe's scaled Sigma Point Algorithm)

            0 <= a <= 1,
            b ~= 2,
            k = 3 - n
        """
        # weights = np.zeros(2*n+1)
        # sigmas = np.zeros((2*n+1, n))
        lambda_ = self.alpha**2 * (self.n + self.kappa) - self.n

        # вычисление весов 
        Wc = np.full(2*self.n+1, 1. / (2*(self.n + lambda_)))
        Wm = np.full(2*self.n+1, 1. / (2*(self.n + lambda_)))
        Wc[0] = lambda_ / (self.n+lambda_) + (1. - self.alpha**2 + self.beta)
        Wm[0] = lambda_ / (self.n + lambda_)

        # вычисление сигма-точек
        sigmas = np.zeros((2*self.n+1, self.n))
        U = cholesky((self.n+lambda_)*self.P)
        sigmas[0] = mu  
        for k in range(self.n):
            sigmas[k+1] = mu + U[k]
            sigmas[self.n+k+1] = mu - U[k]

        self.num_sigmas = sigmas.shape[0]
        return sigmas, Wm, Wc
    
    def unscented_transform(self, sigmas, Wm, Wc):
        """
        
        """
        # вычисление x
        x = np.dot(Wm, sigmas)

        # вычисление P
        kmax, n = np.shape(sigmas)
        P = np.zeros((n,n))
        for k in range(kmax):
            y = sigmas[k] - x
            P += Wc[k] * np.outer(y, y)
        P += self.Q

        return x, P

    def fx(self, x, u):
        """
            Модель системы

            Принимает
            ---------
            x: NumPy array - вектор состояния в момент времени tk-1
            u: list[float] - угловая скорость в момент времени tk

            Возвращает
            ---------
            х: NumPy array - вектор состояния в момент времени tk
        """
        L = self.L_pr * self.to_quat(x)
        # (2*w + 0)/2 = w
        L = integrate_angular_velocity(L, 2*u, [0., 0., 0.], self.dt)
        dL = self.L_pr.inverse() * L
        dL = qt.as_float_array(dL)
        x[0] = 2 * dL[0] * dL[1]
        x[1] = 2 * dL[0] * dL[2]
        x[2] = 2 * dL[0] * dL[3]
        return np.array(x)

    def hx(self, x):
        """
            Функция измерения (measurement function)

            По известному состоянию прогнозирует ожидаемое измерение
        """
        return x.copy()

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

    def get_measurements(self, L, L_CU):
        """
            Метод, преобразовывающий полученный от астродатчика кватернион рассогласования в пространство состояния
        """
        l = L.inverse() * L_CU
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
        w_cor = (self.x - self.x_prev) / self.dt
        return -w_cor * self.k

    def change_R(self, R):
        self.R = R
