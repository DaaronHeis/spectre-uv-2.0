"""
    Математическая модель двигателя-маховика "Агат-40С"

"""

import numpy as np

"""
    Константы
"""

# угол установки ДМ
alpha = 20.0 / 180 * np.pi

# количество ДМ
N = 4

# минимальный кинетический момент (Нмс)
Hm = -40.0
# максимальный кинетический момент (Нмс)
Hp = 40.0

# минимальный электромагнитный момент (Нм)
Mem_min = -0.20
# максимальный электромагнитный момент (Нм)
Mem_max = 0.20

# момент сопротивления (Нм; берется максимальный)
Mc = 0.012

# коэффициенты управления
k_angle = 0.02
k_velocity = 1

# коэффициенты насыщения и инерционности
k_saturation = 1
HH_max = 1

# демпфирующий коэффициент
c = 0


"""
    Класс, отражающий один ДМ
"""


class Flywheel:

    def __init__(self, w_self, H, HH, Md):

        self.J = 1000           # момент инерции
        self.w_self = w_self    # собственная скорость вращения ДМ
        self.H = H              # собственный кин. момент ДМ
        self.HH = HH            # собственная скорость изменения кин. момента ДМ
        self.Md = Md            # получающийся динамический момент ДМ

    # на вход нужно подавать уже приведенные к оси ДМ величины
    def change(self, Mem):

        # скорость изменения кинетического момента
        self.HH = Mem - Mc*np.sign(self.w_self) - c*self.w_self
        self.w_self = self.w_self + self.HH/self.J

        # инерционность
        if self.HH > HH_max:
            self.HH = HH_max
        if self.HH < -HH_max:
            self.HH = -HH_max

        # кинетический момент после поворота
        self.H = self.H + self.HH

        # насыщение
        if self.H > k_saturation * Hp:
            self.H = k_saturation * Hp
        if self.H < k_saturation * Hm:
            self.H = k_saturation * Hm

        # динамический момент
        self.Md = Mem - Mc*np.sign(self.H)


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

    par_out = [0.0, 0.0, 0.0]

    par_out[0] = -sum(par_in)*np.sin(alpha)
    par_out[1] = (par_in[1] - par_in[3])*np.cos(alpha)
    par_out[2] = (par_in[2] - par_in[0])*np.cos(alpha)

    return par_out


# генерирование управляющего момента Mem
# возвращает управляющий момент по каждой оси
def get_Mem(x, w):

    Mem = [0.0, 0.0, 0.0]
    Mem[0] = -float(k_angle*x[0] + k_velocity*w[0])
    Mem[1] = -float(k_angle*x[1] + k_velocity*w[1])
    Mem[2] = -float(k_angle*x[2] + k_velocity*w[2])
    for i in range(3):
        if Mem[i] > Mem_max:
            Mem[i] = Mem_max
        if Mem[i] < Mem_min:
            Mem[i] = Mem_min
    return Mem


# возвращает все требуемые параметры,
# одновременно обновляя каждый ДМ
# dm_all - список, содержащий все объекты ДМ-ов
def get_all(x, w, dm_all):

    Md = [0.0, 0.0, 0.0, 0.0]
    H = [0.0, 0.0, 0.0, 0.0]
    HH = [0.0, 0.0, 0.0, 0.0]

    # генерирование управляющего момента
    Mem = get_Mem(x, w)
    Mem_dm = from_xyz_to_dm(Mem)

    # расчет изменения параметров ДМ
    for i in range(N):

        dm_all[i].change(Mem_dm[i])
        Md[i] = dm_all[i].Md
        H[i] = dm_all[i].H
        HH[i] = dm_all[i].HH

    Md = from_dm_to_xyz(Md)
    H = from_dm_to_xyz(H)
    HH = from_dm_to_xyz(HH)

    return [Md, H, HH]

