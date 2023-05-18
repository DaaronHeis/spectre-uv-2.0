from control_unit import *


class Updater:
    """
        Класс, осуществляющий взаимодействие между модулями модели. По сути, показывает "истинные"
        значения скоростей и ориентаций
        Через этот класс должно осуществляться взаимодействие между БКУ и аппаратурой (ГИВУС, астродатчик и т.д.)
    """

    def __init__(self, t, dt, L_KA, w_KA):
        """
            t: текущее время интегрирования
            dt: такт управления
            L_KA: текущий кватернион ориентации КА
            w_KA: текущая угловая скорость КА
            L_KA_prev: кватернион ориентации КА на прошлом такте
            w_KA_prev: угловая скорость КА на прошлом такте
        """
        self.t = t
        self.dt = dt
        self.L_KA = L_KA
        self.L_KA_prev = L_KA.copy()
        self.w_KA = w_KA
        self.w_KA = w_KA.copy()
        self.w_KA_prev = w_KA

    def update(self, vel, CU:ControlUnit, Astro: AstroSensor, Gyro: GIVUS):
        """
            Один такт работы модели
            В начале такта имеются w_KA и L_KA, в конце интегрируются уравнения модели (в другом файле)
        """

        # получение текущей скорости и ориентации КА
        self.t += self.dt
        self.w_KA_prev = self.w_KA
        self.w_KA = vel
        self.L_KA_prev = self.L_KA
        self.L_KA = integrate_angular_velocity(self.L_KA, self.w_KA, self.w_KA_prev, self.dt)

        # измерение угловой скорости ГИВУСом и получение базовой (измеренной) ориентации ЗД
        w_meas = Gyro.measure_velocity(self.w_KA, self.t, self.dt)
        CU.set_current_velocity(w_meas)

        Astro.set_current_orientation(self.L_KA, CU.is_stabilisation(), self.t)

        # интегрирование уравнения кинематики и астрокоррекция
        CU.set_current_orientation(self.dt, Astro)
        sigma = CU.get_control_moment(self.t)

        # получение параметров модулей и их вывод
        CU_params = CU.get_parameters()
        Astro_params = Astro.get_parameters()
        params = [sigma, CU_params, Astro_params, self.w_KA]

        return params  
