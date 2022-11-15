"""
    Параметры солнечных батарей, использующиеся при расчете их колебательного
    воздействия на движение КА
    Являются константами
    Структура таблиц (для каждой панели)
    j | w^2 | eps | bx | by | bz
    1 |
    2 |
    3 |
"""

solar_panel_1 = [[1.72,  0.013,     0, -0.12, 23.65],
                 [30.06, 0.052, 25.33, -7.42, -0.08],
                 [35.51, 0.057, -0.04,  3.68, -0.19]]

solar_panel_2 = [[1.83,  0.013,     0, 0.09, 23.09],
                 [31.52, 0.054,-24.73,-7.29, 0.08],
                 [36.83, 0.058,     0, 3.61, 0.35]]