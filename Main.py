"""
    Точка запуска моделирования
    Задание начальныз условий, отрисовка графиков
"""

import numpy as np
import quaternion as qt
from init_functions import init_start_orientation, init_target_orientation
from updater import run, OutputData
import matplotlib.pyplot as plt
from cycler import cycler

"""
    ---------------------------------------------
    Выставление начальных условий
    ---------------------------------------------
"""

# коэффициент перехода от градусов к радианам
k = np.pi / 180

t_span_variant = 9

if t_span_variant == 0:
    t_span = [0, 300]
elif t_span_variant == 1:
    t_span = [0, 900]
elif t_span_variant == 2:
    t_span = [0, 1800]
elif t_span_variant == 3:
    t_span = [0, 6000]
else:
    t_span = [0, 500]

# тензор инерции (данные из документации, не менять)
I = np.array([[3379.4, 25.6, 3.2], [25.6, 9283.9, 19.6], [3.2, 19.6, 10578.5]])

# внешний возмущающий постоянный момент
M = [0.0, 0.0, 0.0]

# такт вычислений
dt = 0.25

# начальные и конечные углы и скорости
angles_0 = np.array([00.0, 90.0, 90.0]) * k
angles_end = np.array([0.0, 0.0, 0.0]) * k
vel_0 = np.array([0.0, 0.0, 0.0]) * k
vel_end = np.array([0.0, 0.0, 0.0]) * k

# угол наклона солнечных панелей
gamma = 0

# начальная и конечная ориентация относительно ИСК
L_0 = init_start_orientation(angles_0)
# L_pr = qt.quaternion(1, 0, 0, 0)
L_pr = qt.quaternion(0.86472620667614, 0.256908589358277, 0.334920502035886, 0.272166900113631)
#L_pr = init_target_orientation(angles_end, vel_end)

initial_conditions = [t_span, dt, L_0, L_pr, vel_0, vel_end, M, I, gamma]
required_parameters = ["Углы отклонения от заданной ориентации",
                        #"Угол отклонения оси х аппарата от оси х ИСК"
                        "Проекции вектора угловой скорости",
                        "Проекции измеренной угловой скорости",
                        #"Кинетические моменты каждого ДМ",
                        #"Кинетические моменты в проекциях на оси ССК",
                        #"Скорости вращения ДМ",
                        "Проекции вектора управляющего момента",
                        #"Кватернион рассогласования",
                        "Кватернион текущей ориентации"
                        ]

t, data_ideal, data_full = run(initial_conditions)

n = len(t)

a = data_full.results["Углы отклонения от заданной ориентации"]
a = a[len(a)-1]
print("------------------------------------------------")
print("        Ошибка ориентации по углам              ")
print("По углу крена: ", a[0] * 3600, "угл. мин.")
print("По углу тангажа: ", a[1] * 3600, "угл. мин.")
print("По углу рысканья: ", a[2] * 3600, "угл. мин.")
print("------------------------------------------------")

i = 1
colors3 = cycler(color=['red', 'green', 'blue'])
colors4 = cycler(color=['red', 'green', 'blue', 'yellow'])
for key in required_parameters:
    fig, ax = plt.subplots()
    a_full = data_full.results[key]
    a_ideal = data_ideal.results[key]
    if key == "Углы отклонения от заданной ориентации" or key == "Проекции вектора угловой скорости" or key == "Проекции измеренной угловой скорости":
        a_full = [vector/k for vector in a_full]
        a_ideal = [vector/k for vector in a_ideal]
    # с этого момента все углы в градусах
    b_full = data_full.handles[key]
    b_ideal = data_ideal.handles[key]
    n = len(b_full)
    labels_full = []
    labels_ideal = []
    for j in range(n-1):
        labels_full.append(b_full[j+1])
        labels_ideal.append(b_ideal[j+1]+'_ideal')
    if len(a_full[0]) == 4:
        ax.set_prop_cycle(colors4)
    else:
        ax.set_prop_cycle(colors3)
    ax.plot(t, a_full, label=labels_full)
    ax.plot(t, a_ideal, linestyle='--', label=labels_ideal)
    ax.set_title(key)
    ax.legend()
    ax.set_xlabel("Время, с")
    ax.set_ylabel(b_full[0])
    ax.grid(True)
    i += 1

plt.show()
