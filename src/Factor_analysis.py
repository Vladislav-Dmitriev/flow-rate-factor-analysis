from src.use_cases.FA_calculator import FactorAnalysis
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors

import numpy as np
import pandas as pd
import json

# input_data_plan = {
#     "base_prop": {
#         "calc_type": 'optimal',  # тип расчета (optimal упрощенная модель, desuperposition с эффектом суперпозиции
#         # давления, segmentation сегментация скважины на интервалы)
#         "segment_count": 20,  # кол-во сегментов разбиения скважины
#         "border_type_x": "n",  # тип граничного условия (c - непроницаемая граница, n - свободный поток через границы)
#         "border_type_y": "n",  # тип граничного условия
#         "border_type_z_up": "n",  # тип граничного условия
#         "border_type_z_down": "n",  # тип граничного условия
#     },
#     "target": {
#         "cumulative_time": 2400,  # накопленное время работы, в сутках или часах?
#         "target_values": {
#             "q_liq": 50,  # целевое значение дебита
#             "p_bhp": 150  # целевое значение забойного давления
#         },
#     },
#     "unit": {
#         "skin": 0,  # скин-фактор
#         "h_eff": 14,  # Эффективная толщина пласта
#         "vertical_offset": 0.5,  #
#         "perfres_ratio": 1,  # Вертикальное положение (в долях от мощности) (для горизонтальных)
#         "wellbore_prop": {
#             "wellbore_type": "horizontal",  # Тип скважины: вертикальная, горизонтальная или многозабойная
#             "wellbore_r": 0.1,  # Радиус скважины, м
#             "horizontal_wellbore_length": 624,  # Длина горизонтального ствола
#             "horizontal_wellbore_perf_ratio": 1,  # Доля перфораций горизонтальной скважины
#             "horizontal_perf_count": 5,  # Количество перфорированных интервалов
#             "permeability": 100000000  # Проницаемость скважины
#         },
#         "layer_prop": {
#             "permeability": 12,  # Проницаемость пласта, мД
#             "water_cut": 15,  # Обводненность, %
#             "p_bubble": 100,  # Давление насыщения, атм
#             "kv_kh_ratio": 0.1,  # Отношение вертикальной к горизонтальной проницаемости
#             "compressibility": 0.00015,  # Общая сжимаемость
#             "porosity": 0.2,  # Пористость
#             "p_res_init": 260,  # Начальное пластовое давление
#             "viscosity_eff": 1,  # Вязкость нефти
#             "b_eff": 1.2,  # Объемный коэффициент
#             "res_model_type": "Homogeneous",  # Модель пласта
#             "f_compressibility": 0.0001,  # Сжимаемость трещин
#             "f_porosity": 0.001,  # Пористость трещин
#             "lambda": 0.00001,  # Коэффициент пропускания
#             "internal_r": None,  # Внутренний радиус (для составной модели пласта)
#             "kmu_in_out_ratio": 2,  # (k/mu)_in / (k/mu)_out
#             "kmuphict_in_out_ratio": 2,  # (k/mu_phi_ct)_in / (k/mu_phi_ct)_out
#             "grp_flag": False,  # Параметр, обозначающий, был ли ГРП на скважине
#             "grp_prop": {
#                 "xf": 60,  # Полудлина трещины, м
#                 "kf": 1500,  # Проницаемость трещины, дарси
#                 "wellbore_wf": 4.5,  # Ширина трещины у ствола скважины, мм
#                 "res_wf": 4.31,  # Ширина трещины на конце, мм
#                 "grade": 1,  # Степень аппроксимационной функции
#                 "skin_border": 0,  # Скин на стенке трещины
#                 "skin_ch": 0,  # Скин за счет схождения (Choke Skin)
#                 "fracture_grow_t": 0,  # Время роста трещины, ч
#             },
#             "mgrp_flag": False,  # Параметры МГРП
#             "mgrp_prop": {
#                     "grp_count": 2,  # Количество ГРП
#                     "f_direction": "transverse"  # Направление трещин
#                 },
#             "xe": 1000,  # Длина прямоугольника, м (в напр. Запад - Восток)
#             "ye": 1000,  # Ширина прямоугольника, м (в напр. Север - Юг)
#             "lc_ratio": 0.5,  # Доля от длины до центра
#             "wc_rectangle_ratio": 0.5  # Доля от ширины до центра
#         }
#     }
# }
#
# input_data_fact = {
#     "base_prop": {
#         "calc_type": 'optimal',  # тип расчета (optimal упрощенная модель, desuperposition с эффектом суперпозиции
#         # давления, segmentation сегментация скважины на интервалы)
#         "segment_count": 20,  # кол-во сегментов разбиения скважины
#         "border_type_x": "n",  # тип граничного условия (c - непроницаемая граница, n - свободный поток через границы)
#         "border_type_y": "n",  # тип граничного условия
#         "border_type_z_up": "n",  # тип граничного условия
#         "border_type_z_down": "n",  # тип граничного условия
#     },
#     "target": {
#         "cumulative_time": 2400,  # накопленное время работы, в часах
#         "target_values": {
#             "q_liq": 50,  # целевое значение дебита
#             "p_bhp": 160  # целевое значение забойного давления
#         },
#     },
#     "unit": {
#         "skin": 0,  # скин-фактор
#         "h_eff": 10,  # Эффективная толщина пласта
#         "vertical_offset": 0.5,  #
#         "perfres_ratio": 1,  # Вертикальное положение (в долях от мощности) (для горизонтальных)
#         "wellbore_prop": {
#             "wellbore_type": "horizontal",  # Тип скважины: вертикальная, горизонтальная или многозабойная
#             "wellbore_r": 0.1,  # Радиус скважины, м
#             "horizontal_wellbore_length": 624,  # Длина горизонтального ствола
#             "horizontal_wellbore_perf_ratio": 1,  # Доля перфораций горизонтальной скважины
#             "horizontal_perf_count": 5,  # Количество перфорированных интервалов
#             "permeability": 100000000  # Проницаемость скважины
#         },
#         "layer_prop": {
#             "permeability": 11,  # Проницаемость пласта, мД
#             "water_cut": 20,  # Обводненность, %
#             "p_bubble": 100,  # Давление насыщения, атм
#             "kv_kh_ratio": 0.1,  # Отношение вертикальной к горизонтальной проницаемости
#             "compressibility": 0.00014,  # Общая сжимаемость
#             "porosity": 0.2,  # Пористость
#             "p_res_init": 255,  # Начальное пластовое давление
#             "viscosity_eff": 1,  # Вязкость нефти
#             "b_eff": 1.1,  # Объемный коэффициент
#             "res_model_type": "Homogeneous",  # Модель пласта
#             "f_compressibility": 0.0001,  # Сжимаемость трещин
#             "f_porosity": 0.001,  # Пористость трещин
#             "lambda": 0.00001,  # Коэффициент пропускания
#             "internal_r": None,  # Внутренний радиус (для составной модели пласта)
#             "kmu_in_out_ratio": 2,  # (k/mu)_in / (k/mu)_out
#             "kmuphict_in_out_ratio": 2,  # (k/mu_phi_ct)_in / (k/mu_phi_ct)_out
#             "grp_flag": False,  # Параметр, обозначающий, был ли ГРП на скважине
#             "grp_prop": {
#                 "xf": 50,  # Полудлина трещины, м
#                 "kf": 1400,  # Проницаемость трещины, дарси
#                 "wellbore_wf": 4.4,  # Ширина трещины у ствола скважины, мм
#                 "res_wf": 4.31,  # Ширина трещины на конце, мм
#                 "grade": 1,  # Степень аппроксимационной функции
#                 "skin_border": 0,  # Скин на стенке трещины
#                 "skin_ch": 0,  # Скин за счет схождения (Choke Skin)
#                 "fracture_grow_t": 0,  # Время роста трещины, ч
#             },
#             "mgrp_flag": False,  # Параметры МГРП
#             "mgrp_prop": {
#                     "grp_count": 2,  # Количество ГРП
#                     "f_direction": "transverse"  # Направление трещин
#                 },
#             "xe": 1000,  # Длина прямоугольника, м (в напр. Запад - Восток)
#             "ye": 1000,  # Ширина прямоугольника, м (в напр. Север - Юг)
#             "lc_ratio": 0.5,  # Доля от длины до центра
#             "wc_rectangle_ratio": 0.5  # Доля от ширины до центра
#         }
#     }
# }

# Запись данных в файлы json
# with open("data_plan.json", "w", encoding="utf-8") as f:
#     json.dump(input_data_plan, f, indent=4)
#
# with open("data_fact.json", "w", encoding="utf-8") as f:
#     json.dump(input_data_fact, f, indent=4)

# Считывание входных данных из файлов json
with open("250_GW_fact.json", "r") as f:
    input_data_fact = json.load(f)

with open("250_GW_plan.json", "r") as f:
    input_data_plan = json.load(f)

# for i in range(2, 9):
#     input_data_fact["unit"]["layer_prop"]["mgrp_prop"]["grp_count"] = i
#     # print(input_data_fact)
#     factor_analysis = FactorAnalysis(input_data_fact, input_data_plan, 100)
#     q_fact = factor_analysis.calc_flow_rate(input_data_fact)
#     print(q_fact)

factor_analysis = FactorAnalysis(input_data_fact, input_data_plan, 100)
q_fact = factor_analysis.calc_flow_rate(input_data_fact)
q_plan = factor_analysis.calc_flow_rate(input_data_plan)
print(f" Фактический дебит: {q_fact}", "\n",
      f"Плановый дебит: {q_plan}")
dict_factors = factor_analysis.calc_factors()
print(dict_factors)
print(f"Значения факторов: {dict_factors}", "\n",
      f"Сумма факторов: {sum(value for value in dict_factors.values() if isinstance(value, (int, float)))}", "\n",
      f"Разница дебитов: {q_plan - q_fact}", "\n",
      f"Погрешность ФА: {q_plan - q_fact -
                         sum(value for value in dict_factors.values() if isinstance(value, (int, float)))}")

# import numpy as np
#
#
# def plot_waterfall_chart(data, title="Водопадный график факторного анализа"):
#     """
#     Построение водопадного графика на основе словаря факторного анализа.
#
#     :param data: Словарь {фактор: изменение}, где ключи - названия факторов, значения - их изменения.
#     :param title: Заголовок графика.
#     """
#     labels = list(data.keys())  # Названия факторов
#     values = np.array(list(data.values()))  # Значения факторов
#     cumulative = np.cumsum(values)  # Каскадная сумма факторов
#
#     # Определение точек начала и конца колонок
#     starts = np.insert(cumulative[:-1], 0, 0)
#     ends = cumulative
#
#     # Определение цветов (рост - зеленый, падение - красный, итог - синий)
#     colors = ['green' if v > 0 else 'red' for v in values]
#     colors.append('blue')  # Итоговое значение
#
#     # Добавляем финальное значение (итог)
#     labels.append("Итог")
#     starts = np.append(starts, starts[-1] + values[-1])
#     ends = np.append(ends, starts[-1])
#
#     # Создание графика
#     fig, ax = plt.subplots(figsize=(10, 6))
#     ax.bar(labels, ends - starts, bottom=starts, color=colors, edgecolor="black")
#
#     # Добавляем текстовые подписи
#     for i, (start, end) in enumerate(zip(starts, ends)):
#         ax.text(i, end if values[i - 1] >= 0 else start, f"{values[i - 1]:.2f}", ha='center',
#                 va='bottom' if values[i - 1] >= 0 else 'top')
#
#     ax.set_title(title)
#     ax.set_ylabel("Изменение дебита")
#     ax.set_xticklabels(labels, rotation=45, ha='right')
#
#     plt.show()
#
#
# # Пример данных факторного анализа
# factor_changes = {
#     "Пермеабильность": 10,
#     "Толщина пласта": -5,
#     "Давление пласта": 7,
#     "Давление забоя": -4,
#     "Вязкость": -3,
#     "Обводненность": 6
# }
#
# plot_waterfall_chart(factor_changes)

# # Сортируем по абсолютным значениям
# sorted_data = dict(sorted(dict_factors.items(), key=lambda item: abs(item[1]), reverse=True))
#
# # Преобразуем данные в списки
# labels = list(sorted_data.keys())
# values = np.array(list(sorted_data.values()))
#
# # Нормализация значений для цветовой карты
# norm = mcolors.Normalize(vmin=min(values), vmax=max(values))
# cmap = cm.RdBu  # Выбор цветовой карты (можно поменять)
# colors = [cmap(norm(val)) for val in values]
#
# # Построение торнадо-чарта
# fig, ax = plt.subplots(figsize=(8, 5))
# bars = ax.barh(labels, values, color=colors, edgecolor="black")
#
# # Добавляем цветовую шкалу
# sm = cm.ScalarMappable(cmap=cmap, norm=norm)
# sm.set_array([])
# cbar = plt.colorbar(sm, ax=ax)
# # cbar.set_label("Impact Magnitude")
#
# # Оформление графика
# ax.set_xlabel("Отклонение")
# ax.set_title("Tornado Chart факторный анализ дебита")
# ax.axvline(0, color="black", linewidth=1)  # Центральная ось
#
# plt.show()