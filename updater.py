"""
    Основная модель
    Ideal - идеальная модель без коррекции и ошибок измерений, по сути то, с чем сравнивается полная модель
    Full - полная модель, с фильтром Калмана и ошибками измерений
"""

from init_functions import *
from dynamic_model_equations import runge_kutta
import numpy as np
import quaternion as qt
import flywheel_engine as dm


def to_nparray(data):
    """
        Конвертирует любые данные в подходящий по размеру NumPy-array
    """
    if data is np.array:
        return data
    else:
        return np.array(data)


class OutputData:
    """
        Класс, описывающий параметры, которые необходимо выводить для обработки,
        и все методы обработки этих параметров
        Здесь описаны все необходимые для вывода параметры: сами значения, 
        переменные для отрисовки графиков (оси, легенды, названия и т.д.)

        Идея данного класса в том, чтобы все манипуляции по обработке и представлению
        выходных данных происходили здесь, и во всех остальных местах необходимо было бы 
        предоставлять только ключи и соответствующие им данные
    """

    def __init__(self):
        
        self.results = {}    # значения выходных параметров в виде списка NumPy массивов
        self.handles = {}    # значения для обработки графиков (ед. измерения оси y, легенда)
        self.results["Углы отклонения от заданной ориентации"] = []
        self.results["Угол отклонения оси х аппарата от оси х ИСК"] = []
        self.results["Проекции вектора угловой скорости"] = []
        self.results["Проекции измеренной угловой скорости"] = []
        self.results["Кинетические моменты каждого ДМ"] = []
        self.results["Кинетические моменты в проекциях на оси ССК"] = []
        self.results["Скорости вращения ДМ"] = []
        self.results["Проекции вектора управляющего момента"] = []
        self.results["Кватернион рассогласования"] = []
        self.results["Кватернион текущей ориентации"] = []

        self.handles["Углы отклонения от заданной ориентации"] = ["град", "gamma", "theta", "psi"]
        self.handles["Угол отклонения оси х аппарата от оси х ИСК"] = ["град", "угол"]
        self.handles["Проекции вектора угловой скорости"] = ["град/с","w_x", "w_y", "w_z"]
        self.handles["Проекции измеренной угловой скорости"] = ["град/с","w_x", "w_y", "w_z"]
        self.handles["Кинетические моменты каждого ДМ"] = ["Нм", "H_1", "H_2", "H_3", "H_4"]
        self.handles["Кинетические моменты в проекциях на оси ССК"] = ["Нм", "H_x", "H_y", "H_z"]
        self.handles["Скорости вращения ДМ"] = ["w_1", "w_2", "w_3", "w_4"]
        self.handles["Проекции вектора управляющего момента"] = ["Нм", "sigma_x", "sigma_y", "sigma_z"]
        self.handles["Кватернион рассогласования"] = ["","lambda_0", "lambda_1", "lambda_2", "lambda_3"]
        self.handles["Кватернион текущей ориентации"] = ["","$/lambda$_0", "lambda_1", "lambda_2", "lambda_3"]

    def add_data(self, key, data):
        """
            Добавляет data в нужные массивы данных results в соответствие с key
            Здесь же происходит преобразование каждого типа данных в нужный для добавления в results вид
            Все данные в results лежат в виде списка NumPy-массивов
        """
        if key == "Кватернион рассогласования" or key == "Кватернион текущей ориентации":
            data = to_nparray(qt.as_float_array(data))
            data = data.reshape(4)
        elif data is list:
            data = to_nparray(data)
            data = np.array(data).reshape(len(data))
        else:
            data = np.array(data).reshape(len(data))
        self.results[key].append(data)


def run(initial_conditions):
    """
        initial_conditions:
            - t_span = [t0 t_end] - временной промежуток
            - dt - шаг моделирования
            - L_0, L_pr - начальная и целевая (конечная) ориентация отн. ИСК
            - w_0, w_end - начальные и конечные угловые скорости КА
            - M, I - постоянный внешний момент и тензор инерции
            - gamma - угол наклона солнечных панелей

        required_parameters - список с именами параметров, которые нужно вывести
    """

    # TODO: все выводы проверить
    t_span = initial_conditions[0]
    dt = initial_conditions[1]
    L_0 = initial_conditions[2]
    L_pr = initial_conditions[3]
    w_0 = initial_conditions[4]
    w_pr = initial_conditions[5]
    M0 = initial_conditions[6]
    M_ideal = M0.copy()
    M_full = M0.copy()
    I = initial_conditions[7]
    gamma = initial_conditions[8]

    # инициализация времени
    [t, t_curr, h, t_end] = init_time(t_span, dt)

    # инициализация моделей ДМ
    sigma_max = 0.2
    w_bw = 5
    n = 4
    dmAll_ideal = init_flywheels(n, w_bw, sigma_max)
    dmAll_full = init_flywheels(n, w_bw, sigma_max)

    # инициализация БКУ
    angles_ideal, w_ideal, ctrlUnit_ideal = init_control_unit(L_0, L_pr, w_0, w_pr, w_bw, sigma_max, 
                                                              CORR_KEY=False, ARTIF_ERR_KEY=False)
    angles_full, w_full, ctrlUnit_full = init_control_unit(L_0, L_pr, w_0, w_pr, w_bw, sigma_max, 
                                                           CORR_KEY=True, ARTIF_ERR_KEY=False)

    # инициализация астродатчиков 
    astrosensor_ideal = init_astrosensor(L_0, L_pr, dt, A_S_ERR_KEY=False)
    astrosensor_full = init_astrosensor(L_0, L_pr, dt,  A_S_ERR_KEY=True)

    # инициализация ГИВУСов
    givus_ideal = init_GIVUS(GIVUS_ERR_KEY=False)
    givus_full = init_GIVUS(GIVUS_ERR_KEY=True)

    # инициализация структуры выходных данных
    data_ideal = OutputData()
    data_full = OutputData()

    # -------------------------------------
    # TODO: переписать в функцию
    # получение текущей скорости КА
    wMeas_ideal = givus_ideal.measure_velocity(w_ideal, t_curr, h)
    wMeas_full = givus_full.measure_velocity(w_full, t_curr, h)
    ctrlUnit_ideal.set_current_velocity(wMeas_ideal)
    ctrlUnit_full.set_current_velocity(wMeas_full)

    data_ideal.add_data("Проекции измеренной угловой скорости", wMeas_ideal)
    data_full.add_data("Проекции измеренной угловой скорости", wMeas_full)

    data_ideal.add_data("Проекции вектора угловой скорости", w_ideal)
    data_full.add_data("Проекции вектора угловой скорости", w_full)

    """
        Измерение ориентации астродатчиком производится только для полной модели,
        и представляет из себя кватернион ориентации. полученный идеальной моделью
    """
    # измерение текущей ориентации для идеальной модели
    ctrlUnit_ideal.set_current_orientation(h, astrosensor_ideal)
    L_KA = ctrlUnit_ideal.get_current_orientation()
    data_ideal.add_data("Кватернион текущей ориентации", L_KA)
    # передача полученного кватерниона как "снимка неба" астродатчику полной модели
    astrosensor_full.set_current_orientation(L_KA, ctrlUnit_full.is_stabilisation())
    ctrlUnit_full.set_current_orientation(h, astrosensor_full)
    data_full.add_data("Кватернион текущей ориентации", ctrlUnit_full.get_current_orientation())

    data_ideal.add_data("Углы отклонения от заданной ориентации", ctrlUnit_ideal.get_parameters()[2])
    data_full.add_data("Углы отклонения от заданной ориентации", ctrlUnit_full.get_parameters()[2])

    # формирование управляющих импульсов
    sigma_ideal = ctrlUnit_ideal.get_control_moment(t_curr)
    sigma_full = ctrlUnit_full.get_control_moment(t_curr)
    data_ideal.add_data("Проекции вектора управляющего момента", sigma_ideal)
    data_full.add_data("Проекции вектора управляющего момента", sigma_full)

    # расчёт динамического момента КУДМ
    Md_ideal = dm.update_block(sigma_ideal, dmAll_ideal)
    Md_full = dm.update_block(sigma_full, dmAll_full)

    # получение параметров КУДМ
    dmParam_ideal = dm.get_all(dmAll_ideal)
    dmParam_full = dm.get_all(dmAll_full)
    Md_ideal = dmParam_ideal[0]
    H_ideal = dmParam_ideal[1]
    HH_ideal = dmParam_ideal[2]
    w_dm_ideal = dmParam_ideal[3]
    Md_full = dmParam_full[0]
    H_full = dmParam_full[1]
    HH_full = dmParam_full[2]
    w_dm_full = dmParam_full[3]

    # расчёт действующего момента
    for j in range(3):
        M_ideal[j] = M0[j] + Md_ideal[j]
        M_full[j] = M0[j] + Md_full[j]

    data_ideal.add_data("Кинетические моменты каждого ДМ", H_ideal)
    data_full.add_data("Кинетические моменты каждого ДМ", H_full)
    data_ideal.add_data("Кинетические моменты в проекциях на оси ССК", dm.from_dm_to_xyz(H_ideal))
    data_full.add_data("Кинетические моменты в проекциях на оси ССК", dm.from_dm_to_xyz(H_full))
    data_ideal.add_data("Скорости вращения ДМ", w_dm_ideal)
    data_full.add_data("Скорости вращения ДМ", w_dm_full)
    # ----------------------------------------------

    # инициализация начальных колебательных координат
    q_ideal = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    dq_ideal = q_ideal.copy()
    q_full = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    dq_full = q_full.copy()

    # основной цикл моделирования
    while t_curr < t_end:

        t_curr += h
        t.append(t_curr)

        # получение текущей скорости КА
        wMeas_ideal = givus_ideal.measure_velocity(w_ideal, t_curr, h)
        wMeas_full = givus_full.measure_velocity(w_full, t_curr, h)
        ctrlUnit_ideal.set_current_velocity(wMeas_ideal)
        ctrlUnit_full.set_current_velocity(wMeas_full)

        data_ideal.add_data("Проекции измеренной угловой скорости", wMeas_ideal)
        data_full.add_data("Проекции измеренной угловой скорости", wMeas_full)

        data_ideal.add_data("Проекции вектора угловой скорости", w_ideal)
        data_full.add_data("Проекции вектора угловой скорости", w_full)

        """
            Измерение ориентации астродатчиком производится только для полной модели,
            и представляет из себя кватернион ориентации. полученный идеальной моделью
        """
        # измерение текущей ориентации для идеальной модели
        ctrlUnit_ideal.set_current_orientation(h, astrosensor_ideal)
        L_KA = ctrlUnit_ideal.get_current_orientation()
        data_ideal.add_data("Кватернион текущей ориентации", L_KA)
        # передача полученного кватерниона как "снимка неба" астродатчику полной модели
        astrosensor_full.set_current_orientation(L_KA, ctrlUnit_full.is_stabilisation())
        ctrlUnit_full.set_current_orientation(h, astrosensor_full)
        data_full.add_data("Кватернион текущей ориентации", ctrlUnit_full.get_current_orientation())

        data_ideal.add_data("Углы отклонения от заданной ориентации", ctrlUnit_ideal.get_parameters()[2])
        data_full.add_data("Углы отклонения от заданной ориентации", ctrlUnit_full.get_parameters()[2])

        # формирование управляющих импульсов
        sigma_ideal = ctrlUnit_ideal.get_control_moment(t_curr)
        sigma_full = ctrlUnit_full.get_control_moment(t_curr)
        data_ideal.add_data("Проекции вектора управляющего момента", sigma_ideal)
        data_full.add_data("Проекции вектора управляющего момента", sigma_full)

        # расчёт динамического момента КУДМ
        Md_ideal = dm.update_block(sigma_ideal, dmAll_ideal)
        Md_full = dm.update_block(sigma_full, dmAll_full)

        # получение параметров КУДМ
        dmParam_ideal = dm.get_all(dmAll_ideal) 
        dmParam_full = dm.get_all(dmAll_full)
        Md_ideal = dmParam_ideal[0]
        H_ideal = dmParam_ideal[1]
        HH_ideal = dmParam_ideal[2]
        w_dm_ideal = dmParam_ideal[3]
        Md_full = dmParam_full[0]
        H_full = dmParam_full[1]
        HH_full = dmParam_full[2]
        w_dm_full = dmParam_full[3]

        # расчёт действующего момента
        for j in range(3):
            M_ideal[j] = M0[j] + Md_ideal[j]
            M_full[j] = M0[j] + Md_full[j]

        data_ideal.add_data("Кинетические моменты каждого ДМ", H_ideal)
        data_full.add_data("Кинетические моменты каждого ДМ", H_full)
        data_ideal.add_data("Кинетические моменты в проекциях на оси ССК", dm.from_dm_to_xyz(H_ideal))
        data_full.add_data("Кинетические моменты в проекциях на оси ССК", dm.from_dm_to_xyz(H_full))
        data_ideal.add_data("Скорости вращения ДМ", w_dm_ideal)
        data_full.add_data("Скорости вращения ДМ", w_dm_full)

        # интегрирование уравнений
        [w_ideal, q_ideal, dq_ideal] = runge_kutta(h, wMeas_ideal, q_ideal, dq_ideal, I, HH_ideal, H_ideal, M_ideal, gamma)
        w_ideal = np.array(w_ideal).reshape(3)
        q_ideal = np.array(q_ideal).reshape(6)
        dq_ideal = np.array(dq_ideal).reshape(6)
        [w_full, q_full, dq_full] = runge_kutta(h, wMeas_full, q_full, dq_full, I, HH_full, H_full, M_full, gamma)
        w_full = np.array(w_full).reshape(3)
        q_full = np.array(q_full).reshape(6)
        dq_full = np.array(dq_full).reshape(6)

        # вывод текущего времени каждые 1000 секунд, чтобы знать, что цикл не завис
        if t_curr % 1000 == 0:
            print('t_curr = ', t_curr)

    return t, data_ideal, data_full
      




