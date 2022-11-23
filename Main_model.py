"""
    Основные уравнения динамики КА и их решение
"""

import numpy as np
import quaternion as qt
import matplotlib.pyplot as plt
import flywheel_engine as dm
import control_unit as ctrl
import solar_panels_constants as spc
import updater as upd


def f(w, M, H, HH, I):
    # TODO: добавить колебательность
    """
        Диффуры динамики КА
        Представляют из себя следующую систему:
            w' = A * [a] @ [sigma]
            phi' = w
            [a] и [sigma] - матрица 3х3 и 3х1 соответственно
            [a] - матрица моментов инерции
            [sigma] - вектор-столбец действующих моментов
    """
    # ww = np.array([[0.0], [0.0], [0.0], [0.0], [0.0], [0.0]])
    A = 1/(I[0,0]*I[1,1]*I[2,2] - I[0,0]*I[1,2]**2 - I[1,1]*I[0,2]**2 - I[2,2]*I[0,1]**2 - 2*I[0,1]*I[0,2]*I[1,2])
    sigma = np.array([0.0, 0.0, 0.0])
    sigma[0] = M[0] - (HH[0] + H[2]*w[1] - H[1]*w[2])
    sigma[1] = M[1] - (HH[1] + H[0]*w[2] - H[2]*w[0])
    sigma[2] = M[2] - (HH[2] + H[1]*w[0] - H[0]*w[1])
    a = np.array([[I[1,1]*I[2,2]-I[1,2]**2, I[2,2]*I[0,1]+I[0,2]*I[1,2], I[1,1]*I[0,2]+I[0,1]*I[1,2]],
                  [I[2,2]*I[0,1]+I[0,2]*I[1,2], I[0,0]*I[2,2]-I[0,2]**2, I[0,0]*I[1,2]+I[0,1]*I[0,2]],
                  [I[1,1]*I[0,2]+I[0,1]*I[1,2], I[0,0]*I[1,2]+I[0,1]*I[0,2], I[0,0]*I[1,1]-I[0,1]**2]])
    ww = A * (a @ sigma.T)
    return ww


def init_time(t_span, dt):
    """ Инициализация временного отрезка моделирования """
    t_begin = t_span[0]
    t_end = t_span[1]
    h = dt
    t = [t_begin]
    t_curr = t_begin + h
    return t, t_curr, h, t_end


def init_flywheels(n, w_bw, sigma_max):
    """ Инициализация двигателей-маховиков """
    dm_all = []
    for j in range(n):
        dm_all.append(dm.Flywheel(0,0,0,0,w_bw))
    return dm_all


def init_target_orientation(angles_end, vel_end):
    """ Выставление требуемой ориентации через заданные углы """
    omega_pr = vel_end.copy()
    gamma_pr = angles_end.copy()
    l_k = ctrl.from_euler_to_quat(gamma_pr, 'YZXr')
    # l_pr = l_k.inverse() * l_0
    #l_pr = qt.quaternion(1, 0, 0, 0)
    l_pr = qt.quaternion(0.86472620667614, 0.256908589358277, 0.334920502035886, 0.272166900113631)
    print('l_pr = ', l_pr)
    return l_pr, omega_pr


def init_start_orientation(angles_0):
    """ Выставление начальной ориентации (сразу в виде кватерниона) """
    angles = angles_0
    l_0 = ctrl.from_euler_to_quat(angles, 'YZXr')
    print('l_0 = ', l_0)
    return l_0


def init_control_unit(l_0, l_pr, vel_0, omega_pr, w_bw, sigma_max,
                      CORR_KEY: bool, ARTIF_ERR_KEY: bool):
    """ Инициализация модуля блока управления """
    vel = vel_0.copy()
    l_cur = l_0.copy()
    l_delta = l_pr.inverse() * l_cur
    print('l_delta = ', l_delta)
    angles = np.array([2 * l_delta.w * l_delta.x, 2 * l_delta.w * l_delta.y, 2 * l_delta.w * l_delta.z])
    ctrl_unit = ctrl.ControlUnit(l_pr, l_cur, l_delta, vel, angles, omega_pr,
                                 w_bw, sigma_max, CORR_KEY, ARTIF_ERR_KEY)
    return angles, vel, ctrl_unit


def init_GIVUS(GIVUS_ERR_KEY: bool):
    """ Инициализация модуля ГИВУСа """
    givus = ctrl.GIVUS(GIVUS_ERR_KEY)
    return givus


def init_astrosensor(l_0, A_S_ERR_KEY: bool):
    """ Инициализация модуля астродатчика """
    astrosensor = ctrl.AstroSensor(l_0, np.array([0,0,0]), A_S_ERR_KEY)
    return astrosensor

def init_handler(t, dt, L_KA, w_KA):
    """ Инициализация интерфейса - updater-a """
    handler = upd.Updater(0,dt,L_KA,w_KA)
    return handler


def run(t_span, dt, angles_0, angles_end, vel_0, vel_end, M0, I, CORR_KEY, A_S_ERR_KEY, GIVUS_ERR_KEY, ARTIF_ERR_KEY):
    """
        t_span: float[1x2] - начальный и конечный моменты времени
        dt: float
        angles_0: np.array(float[1x3]) - начальные углы в радианах
        angles_end: np.array(float[1x3]) - конечные углы в радианах (пока не используются)
        vel_0: np.array(float[1x3]) - начальные угловые скорости
        vel_end: np.array(float[1x3]) - конечные угловые скорости
        M0: float[1x3] - внешние возмущающие моменты
        I: np.array(float[3x3]) - тензор инерции
    """

    # инициализация времени
    [t, t_curr, h, t_end] = init_time(t_span, dt)

    # инициализация моделей ДМ
    # TODO: учесть ограничения на скорость из-за астродатчика (не более 0.2 град/с и 0.1 град/с^2 отн. ЭСК)
    sigma_max = 0.2
    w_bw = 5
    n = 4
    dm_all = init_flywheels(n, w_bw, sigma_max)

    # инициалиация ориентации
    [l_pr, omega_pr] = init_target_orientation(angles_end, vel_end)
    l_0 = init_start_orientation(angles_0)

    """
        -------------------------------------
        Инициализация модуля управления здесь
        -------------------------------------
    """
    [angles, vel, ctrl_unit] = init_control_unit(l_0, l_pr, vel_0, omega_pr, w_bw, sigma_max,
                                                 CORR_KEY, ARTIF_ERR_KEY)

    astrosensor = init_astrosensor(l_0, A_S_ERR_KEY)
    givus = init_GIVUS(GIVUS_ERR_KEY)

    handler = init_handler(t, dt, l_0, vel_0)

    # начальные условия
    k = 0

    """
        ---------------------------------------------
        Инициализирование словаря с выходными данными
        ---------------------------------------------
    """

    H_dm_out = np.array([0.0, 0.0, 0.0, 0.0])
    H_xyz_out = np.array([0.0, 0.0, 0.0])
    HH_dm_out = np.array([0.0, 0.0, 0.0, 0.0])
    HH_xyz_out = np.array([0.0, 0.0, 0.0])
    sigma_out = np.array([0.0, 0.0, 0.0])

    params = ctrl_unit.get_parameters()
    a = qt.as_float_array(params[0])
    l_cur_out = np.array([a[0], a[1], a[2], a[3]])
    b = qt.as_float_array(params[1])
    l_delta_out = np.array([b[0], b[1], b[2], b[3]])

    # vel_noise_out = np.array(vel)
    angles_out = np.array(angles)
    vel_out = np.array(vel)
    w_dm_out = np.array([0.0, 0.0, 0.0, 0.0])

    results = {}
    results["Углы отклонения от заданной ориентации"] = [angles_out]
    results["Проекции вектора угловой скорости на оси ССК"] = [vel_out]
    # results["Зашумленные проекции вектора угловой скорости на оси ССК"] = [vel_noise_out]
    # results["Кинетические моменты каждого ДМ"] = [H_dm_out]
    # results["Кинетические моменты в проекциях на оси ССК"] = [H_xyz_out]
    # results["Скорость изменения кин. момента каждого ДМ"] = [HH_dm_out]
    # results["Скорость изменения кин. момента в проекциях на оси ССК"] = [HH_xyz_out]
    # results["Скорости вращения ДМ"] = [w_dm_out]
    # results["Проекции вектора управляющего момента на оси ССК"] = [sigma_out]
    results["Кватернион рассогласования"] = [l_delta_out]
    # results["Кватернион текущей ориентации"] = [l_cur_out]

    handles = {}
    handles["Углы отклонения от заданной ориентации"] = ["град", "gamma", "theta", "psi"]
    handles["Проекции вектора угловой скорости на оси ССК"] = ["град/с","w_x", "w_y", "w_z"]
    # handles["Зашумленные проекции вектора угловой скорости на оси ССК"] = ["град/с","w_x", "w,y", "w_z"]
    # handles["Кинетические моменты каждого ДМ"] = ["H_1", "H_2", "H_3", "H_4"]
    # handles["Кинетические моменты в проекциях на оси ССК"] = ["H_x", "H_y", "H_z"]
    # handles["Скорость изменения кин. момента каждого ДМ"] = [HH_dm_out]
    # handles["Скорость изменения кин. момента в проекциях на оси ССК"] = [HH_xyz_out]
    # handles["Скорости вращения ДМ"] = ["w_1", "w_2", "w_3", "w_4"]
    # handles["Проекции вектора управляющего момента на оси ССК"] = ["Нм", "sigma_x", "sigma_y", "sigma_z"]
    handles["Кватернион рассогласования"] = ["","lambda_0", "lambda_1", "lambda_2", "lambda_3"]
    # handles["Кватернион текущей ориентации"] = ["","lambda_0", "lambda_1", "lambda_2", "lambda_3"]

    # интегрирование
    # TODO: переписать модульно, с учетом добавления новых модулей в будущем

    while t_curr <= t_end:

        t.append(t_curr)

        params = handler.update(vel, ctrl_unit, astrosensor, givus)
        sigma = params[0]
        sigma = [sigma[0], sigma[1], sigma[2]]

        # получение текущих параметров модуля управления
        ctrl_unit_params = params[1].copy()
        a = qt.as_float_array(ctrl_unit_params[0])
        l_cur = np.array([a[0], a[1], a[2], a[3]])
        b = qt.as_float_array(ctrl_unit_params[1])
        l_delta = np.array([b[0], b[1], b[2], b[3]])
        angles = ctrl_unit_params[2].copy()

        # расчет динамического момента комплекса двигателей-маховиков
        Md = dm.update_block(sigma, dm_all)

        # получение параметров комплекса двигателей-маховиков
        dm_param = dm.get_all(dm_all)
        H_dm = []
        HH_dm = []
        w_self_dm = []
        for param in dm_param:
            H_dm.append(param[1])
            HH_dm.append(param[2])
            w_self_dm.append(param[3])
        H_xyz = dm.from_dm_to_xyz(H_dm)
        HH_xyz = dm.from_dm_to_xyz(HH_dm)
        H = H_xyz.copy()
        HH = HH_xyz.copy()

        # расчет действующего момента
        for j in range(3):
            M[j] = M0[j] + Md[j]

        # коэффициенты для Рунге-Кутты
        k1 = f(vel, M, H, HH, I)
        k2 = f(vel + h*k1/2, M, H, HH, I)
        k3 = f(vel + h*k2/2, M, H, HH, I)
        k4 = f(vel + h*k3, M, H, HH, I)

        # формирование выходных массивов
        vel = vel + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        results["Углы отклонения от заданной ориентации"].append(angles)
        results["Проекции вектора угловой скорости на оси ССК"].append(vel)
        # results["Зашумленные проекции вектора угловой скорости на оси ССК"].append(vel_noise)
        # results["Кинетические моменты каждого ДМ"].append(np.array(H_dm))
        # results["Кинетические моменты в проекциях на оси ССК"].append(H_xyz)
        # results["Скорости вращения ДМ"].append(np.array(w_self_dm))
        # results["Проекции вектора управляющего момента на оси ССК"].append(np.array(sigma))
        results["Кватернион рассогласования"].append(l_delta)
        # results["Кватернион текущей ориентации"].append(l_cur)

        k += 1
        t_curr += h
        # вывод текущего времени каждые 1000 секунд, чтобы знать, что цикл не завис
        if t_curr % 1000 == 0:
            print('t_curr = ', t_curr)

    return t, results, handles


if __name__ == '__main__':

    # время интегрирования
    # TODO: написать нормальный интерфейс вместо постоянной правки кода

    t_span_variant = 1

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

    # такт вычислений
    dt = 0.25

    # коэффициент перехода от градусов к радианам
    k = np.pi / 180

    """
        Небольшое пояснение по поводу физического смысла вектора состояния
        Вектор состояния состоит из трех угловых скоростей и трех углов (по каждой оси)
        Углы получаются путем интегрирования угловой скорости при решении 
        кинематических уравнений, а значит, они являются углами поворота вокруг каждой оси
        Скорости измеряются БИУС, которые, по идее, измеряют их относительно ИСК, а значит, и углы
        будут являться углами вращения вокруг каждой из осей ИСК
        Изначальные и конечные углы задаются поворотами в порядке 'YZXr', то есть сначала
        поворот по тангажу (ось Y), затем по рысканью (вокруг оси Z'), и затем по крену (вокруг оси X'')  
    
        Значения внутри векторов - в град/с и град, но сразу же после этого они
        пересчитываются в радианы  
        
        "Красивые" начальные углы: [35.0, 40.0, 20.0]   
    """
    angles_0 = np.array([00.0, 90.0, 90.0]) * k
    angles_end = np.array([0.0, 0.0, 0.0]) * k
    vel_0 = np.array([0.0, 0.0, 0.0]) * k
    vel_end = np.array([0.0, 0.0, 0.0]) * k

    # внешний возмущающий постоянный момент
    M = [0.0, 0.0, 0.0]

    # с этого момента все углы в радианах

    # тензор инерции (данные из документации, не менять)
    I = np.array([[3337.94, 25.6, 3.2], [25.6, 9283.9, 19.6], [3.2, 19.6, 10578.5]])

    """
            !!!ВАЖНО!!!
        Выходные углы - это углы, которые получаются из текущего кватерниона ориентации
        Они (по идее) отображают углы, на которые осталось повернуться, чтобы прийти к итоговой ориентации
    """
    [t, results, handles] = run(t_span, dt, angles_0, angles_end, vel_0, vel_end, M, I,
                                CORR_KEY=True,
                                A_S_ERR_KEY=True,
                                GIVUS_ERR_KEY=True,
                                ARTIF_ERR_KEY=True)

    # отображение графиков
    n = len(t)
    # tmp = [vector/k for vector in results["angles"]]
    results["Углы отклонения от заданной ориентации"] = \
        [vector/k for vector in results["Углы отклонения от заданной ориентации"]]
    results["Проекции вектора угловой скорости на оси ССК"] = \
        [vector/k for vector in results["Проекции вектора угловой скорости на оси ССК"]]
    """ results["Зашумленные проекции вектора угловой скорости на оси ССК"] = \
        [vector/k for vector in results["Зашумленные проекции вектора угловой скорости на оси ССК"]]"""

    # с этого момента все углы в градусах

    # включение LaTeX
    # plt.rcParams.update({"text.usetex":True})

    a = results["Углы отклонения от заданной ориентации"]
    a = a[len(a)-1]
    print("------------------------------------------------")
    print("        Ошибка ориентации по углам              ")
    print("По углу крена: ", a[0]*60, "угл. мин.")
    print("По углу тангажа: ", a[1] * 60, "угл. мин.")
    print("По углу рысканья: ", a[2] * 60, "угл. мин.")
    print("------------------------------------------------")

    i = 1
    for u in results.keys():
        fig = plt.figure(i)
        a = results[u]
        b = handles[u]
        n = len(b)
        labels = []
        for j in range(n-1):
            labels.append(b[j+1])
        plt.plot(t, a, label=labels)
        plt.title(u)
        plt.legend()
        plt.xlabel("Время, с")
        plt.ylabel(b[0])
        plt.grid(True)
        i += 1
    plt.show()
