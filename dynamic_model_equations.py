"""
    Набор функций, описывающих интегрирование уравнений динамики КА
"""

import numpy as np
from numpy.linalg import inv, det
from solar_panel_constants import *


def f(x, parameters):
    """
        Решение уравнений динамики КА в матричном виде

                A@dx - B = 0

                A@dx = B

                dx = A^-1 @ B
        
        Вектор х = [wx wy wz q11 q12 q13 q21 q22 q23 p11 p12 p13 p21 p22 p23]'

        Уравнения динамики:

        Ix*dwx - Ixy*dwy - Ixz*dwz + s41*dp11 + s51*dp12 + s61*dp13 + s71*dp21 + s81*dp22 + s91*dp23 = -s1;

        -Ixy*dwx + Iy*dwy - Iyz*dwz + s42*dp11 + s52*dp12 + s62*dp13 + s72*dp21 + s82*dp22 + s92*dp23 = -s2;

        -Ixz*dwx - Iyz*dwy + Iz*dwz + s43*dp11 + s53*dp12 + s63*dp13 + s73*dp21 + s83*dp22 + s93*dp23 = -s3;

        s41*dwx + d42*dwy + s43*dwz + eps11*dq11 + dp11 = -o11*q11;

        s51*dwx + d52*dwy + s53*dwz + eps12*dq12 + dp12 = -o12*q12;

        s61*dwx + d62*dwy + s63*dwz + eps13*dq13 + dp13 = -o13*q13;

        s71*dwx + d72*dwy + s73*dwz + eps21*dq21 + dp21 = -o21*q21;

        s81*dwx + d82*dwy + s83*dwz + eps22*dq22 + dp22 = -o22*q22;

        s91*dwx + d92*dwy + s93*dwz + eps23*dq23 + dp23 = -o23*q23;

        dq11 = p11;

        dq12 = p12;

        dq13 = p13;

        dq21 = p21;

        dq22 = p22;

        dq23 = p23;
        
        Возвращает
        ----------
            dx - производные компонент вектора х
    """

    """
        Переименовывание для читаемости
    """
    """
    wx = w[0]
    wy = w[1]
    wz = w[2]
    q11 = q[0][0]
    q12 = q[0][1]
    q13 = q[0][2]
    q21 = q[1][0]
    q22 = q[1][1]
    q23 = q[1][2]
    p11 = p[0][0]
    p12 = p[0][1]
    p13 = p[0][2]
    p21 = p[1][0]
    p22 = p[1][1]
    p23 = p[1][2]
    """
    wx = x[0, 0]
    wy = x[1, 0]
    wz = x[2, 0]
    q11 = x[3, 0]
    q12 = x[4, 0]
    q13 = x[5, 0]
    q21 = x[6, 0]
    q22 = x[7, 0]
    q23 = x[8, 0]
    p11 = x[9, 0]
    p12 = x[10, 0]
    p13 = x[11, 0]
    p21 = x[12, 0]
    p22 = x[13, 0]
    p23 = x[14, 0]

    I = parameters['I']
    dH = parameters['dH']
    H = parameters['H']
    M = parameters['M']
    g = parameters['g']

    Ix = I[0, 0]
    Iy = I[1, 1]
    Iz = I[2, 2]
    Ixy = I[0, 1]
    Ixz = I[0, 2]
    Iyz = I[1, 2]
    dHx = dH[0]
    dHy = dH[1]
    dHz = dH[2]
    Hx = H[0]
    Hy = H[1]
    Hz = H[2]
    Mx = M[0]
    My = M[1]
    Mz = M[2]

    """ Промежуточные коэффициенты """
    s1 = (Iz-Iy)*wz*wy + Iyz*(wz**2-wy**2) + Ixy*wx*wz-Ixz*wx*wy + dHx+Hz*wy-Hy*wz - Mx
    s2 = (Ix-Iz)*wx*wz + Ixz*(wx**2-wy**2) - Ixy*wy*wz+Iyz*wx*wy + dHy+Hx*wz-Hz*wx - My
    s3 = (Iy-Ix)*wx*wy + Ixy*(wy**2-wx**2) + Ixz*wy*wz-Iyz*wx*wz + dHz+Hy*wx-Hx*wy - Mz

    s41 = bx11*np.cos(g) + bz11*np.sin(g)
    s42 = by11*np.sin(g)
    s43 = -bx11*np.sin(g) + bz11*np.cos(g)

    s51 = bx12*np.cos(g) + bz12*np.sin(g)
    s52 = by21*np.sin(g)
    s53 = -bx12*np.sin(g) + bz12*np.cos(g)

    s61 = bx13*np.cos(g) + bz13*np.sin(g)
    s62 = by12*np.cos(g)
    s63 = -bx13*np.sin(g) + bz13*np.cos(g)

    s71 = bx21*np.cos(g) + bz21*np.sin(g)
    s72 = by22*np.cos(g)
    s73 = -bx21*np.sin(g) + bz21*np.cos(g)

    s81 = bx22*np.cos(g) + bz22*np.sin(g)
    s82 = by13
    s83 = -bx22*np.sin(g) + bz22*np.cos(g)

    s91 = bx23*np.cos(g) + bz23*np.sin(g)
    s92 = by23
    s93 = -bx23*np.sin(g) + bz23*np.cos(g)

    """ Запись матриц А и В """
    # х = [wx wy wz q11 q12 q13 q21 q22 q23 p11 p12 p13 p21 p22 p23]'
    A = np.array([
        [Ix, -Ixy, -Ixz, 0, 0, 0, 0, 0, 0, s41, s51, s61, s71, s81, s91],
        [-Ixy, Iy, -Iyz, 0, 0, 0, 0, 0, 0, s42, s52, s62, s72, s82, s92],
        [-Ixz, -Iyz, Iz, 0, 0, 0, 0, 0, 0, s43, s53, s63, s73, s83, s93],
        [s41, s42, s43, eps11, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
        [s51, s52, s53, 0, eps12, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
        [s61, s62, s63, 0, 0, eps13, 0, 0, 0, 0, 0, 1, 0, 0, 0],
        [s71, s72, s73, 0, 0, 0, eps21, 0, 0, 0, 0, 0, 1, 0, 0],
        [s81, s82, s83, 0, 0, 0, 0, eps22, 0, 0, 0, 0, 0, 1, 0],
        [s91, s92, s93, 0, 0, 0, 0, 0, eps23, 0, 0, 0, 0, 0, 1],
        [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]
        ])

    B = np.array([
       [-s1],
       [-s2],
       [-s3],
       [-o11*q11],
       [-o12*q12],
       [-o13*q13],
       [-o21*q21],
       [-o22*q22],
       [-o23*q23],
       [p11],
       [p12],
       [p13],
       [p21],
       [p22],
       [p23]
       ])

    dx = inv(A) @ B

    return dx


def runge_kutta(h, w, q, dq, I, dH, H, M, gamma):
    """ 
        Интегрирование уравнений динамики методом Рунге-Кутты 

        Возвращает
        ----------
        проинтегрированные
            w - вектор угловой скорости КА

            q - колебательные координаты КА

            dq - скорость изменения колебательных координат КА    
    """

    # параметры КА, остающиеся постоянными при интегрировании, т.е. считаются константами
    const_parameters = {'I': I, 'dH': dH, 'H': H, 'M': M, 'g': gamma}
    x = np.array([[w[0]],[w[1]],[w[2]],
                [q[0]],[q[1]],[q[2],],[q[3]],[q[4]],[q[5]],
                [dq[0]],[dq[1]],[dq[2],],[dq[3]],[dq[4]],[dq[5]]])
    k1 = f(x, const_parameters)
    k2 = f(x + h*k1/2, const_parameters)
    k3 = f(x + h*k2/2, const_parameters)
    k4 = f(x + h*k3, const_parameters)

    x = x + h/6 * (k1 + 2*k2 + 2*k3 + k4)
    x = x.tolist()
    w = x[0:3]
    q = x[3:9]
    dq = x[9:15]

    return [w, q, dq]
