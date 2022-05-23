"""
    Математическая модель двигателя-маховика "Агат-40С"

"""

import numpy as np

"""
    Константы, относящиеся ко всему комплексу
"""

# угол установки ДМ
alpha = 20.0 / 180 * np.pi

# матрица перехода от ДМ к ССК
Aw = np.array([[-np.sin(alpha), -np.sin(alpha), -np.sin(alpha), -np.sin(alpha)],
               [0.0, np.cos(alpha), 0.0, -np.cos(alpha)],
               [-np.cos(alpha), 0.0, np.cos(alpha), 0.0]])

# количество ДМ
N = 4

# минимальный электромагнитный момент (Нм)
Mem_min = -0.20
# максимальный электромагнитный момент (Нм)
Mem_max = 0.20


class Flywheel:
    """
        Класс, отражающий один ДМ
    """

    def __init__(self, w_self, H, HH, Md, w_bw):
        """
            Входные параметры:
                w_self: текущая собственная скорость вращения
                H: текущий собственный кинетический момент
                HH: текущая скорость изменения кинетического момента
                Md: текущий генерируемый динамический момент
                w_bw: максимальная скорость изменения угловой скорости
        """

        # динамические параметры ДМ
        self.w_self = w_self            # собственная скорость вращения ДМ
        self.H = H                      # собственный кин. момент ДМ
        self.HH = HH                    # собственная скорость изменения кин. момента ДМ
        self.Md = Md                    # получающийся динамический момент ДМ

        # статические параметры ДМ
        self.c = 0.020                  # демпфирующий коэффициент
        self.J = 1.5                    # момент инерции, кгм^2
        self.k_saturation = 2           # коэффициент ограничения
        self.HH_max = w_bw * self.J     # коэффициент насыщения
        self.HH_nonresponse = 0.02      # коэффициент нечувствительности
        self.Mc = 0.012                 # момент сопротивления (Нм; берется максимальный)
        self.Hm = -40                   # минимальный кинетический момент (Нмс)
        self.Hp = 40                    # максимальный кинетический момент (Нмс)

    # на вход нужно подавать уже приведенные к оси ДМ величины
    def change(self, Mem):

        # скорость изменения кинетического момента
        self.HH = Mem - self.Mc*np.sign(self.w_self)*abs(self.H/self.Hm) - self.c*self.w_self

        # насыщение
        if self.HH > self.HH_max:
            self.HH = self.HH_max
        if self.HH < -self.HH_max:
            self.HH = -self.HH_max

        """
        # нечувствительность
        if abs(self.HH) <= self.HH_nonresponse:
            self.HH = 0
        """

        self.w_self = self.w_self + self.HH/self.J

        # кинетический момент после поворота
        # self.H = self.H + self.HH
        self.H = self.J * self.w_self

        # ограничение по максимальному моменту
        if self.H > self.k_saturation * self.Hp:
            self.H = self.k_saturation * self.Hp
        if self.H < self.k_saturation * self.Hm:
            self.H = self.k_saturation * self.Hm

        # динамический момент
        self.Md = Mem - self.Mc*np.sign(self.H)*abs(self.H/self.Hm) - self.c*self.w_self

    def get_param(self):
        """
            Получение текущих динамических параметров ДМ
            Возвращает:
                Md
                H
                HH
                w
        """
        return self.Md, self.H, self.HH, self.w_self


"""
    Функции, описывающие работу комплекса ДМ
"""


# перевод параметров из ССК на оси ДМ
# par_in = [a, b, c]
def from_xyz_to_dm(par_in):

    par_out = [0.0, 0.0, 0.0, 0.0]
    par_x1 = par_in[0] / 2
    par_x2 = par_x1
    par_y = par_in[1]
    par_z = par_in[2]
    par_out[0] = -(par_z/2/np.cos(alpha) + par_x1/2/np.sin(alpha))
    par_out[1] = par_y/2/np.cos(alpha) - par_x2/2/np.sin(alpha)
    par_out[2] = par_z/2/np.cos(alpha) - par_x1/2/np.sin(alpha)
    par_out[3] = -(par_y/2/np.cos(alpha) + par_x2/2/np.sin(alpha))

    return par_out


# перевод параметров из осей ДМ в ССК
# par_in = [a, b, c, d]
def from_dm_to_xyz(par_in):

    par_out = Aw @ np.array(par_in).T
    """
    par_out = [0.0, 0.0, 0.0]

    par_out[0] = -sum(par_in)*np.sin(alpha)
    par_out[1] = (par_in[1] - par_in[3])*np.cos(alpha)
    par_out[2] = (par_in[2] - par_in[0])*np.cos(alpha)
    """

    return par_out


# одновременно обновляя каждый ДМ
# dm_all - список, содержащий все объекты ДМ-ов
def update_block(sigma, dm_all):

    Md = [0.0, 0.0, 0.0, 0.0]

    # генерирование управляющего момента
    Mem = sigma
    """
    for i in range(len(Mem)):
        if Mem[i] > Mem_max:
            Mem[i] = Mem_max
        if Mem[i] < Mem_min:
            Mem[i] = Mem_min
    """
    Mem_dm = from_xyz_to_dm(Mem)

    # расчет изменения параметров ДМ
    for i in range(N):
        dm_all[i].change(Mem_dm[i])
        [Md[i], _, _, _] = dm_all[i].get_param()
    Md = from_dm_to_xyz(Md)
    return Md


def get_all(dm_all):
    """
        Возвращает текущие динамические параметры каждого ДМ в блоке
        НЕ ПЕРЕВОДИТ в проекции на оси координат
        param - массив 4x4, строки - это ДМы,
                            столбцы - это list [Md, H, HH, w]
    """
    param = []
    for dm in dm_all:
        [Md, H, HH, w_self] = dm.get_param()
        param.append([Md, H, HH, w_self])

    return param
