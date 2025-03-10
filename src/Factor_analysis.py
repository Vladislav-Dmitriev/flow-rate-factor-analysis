from src.use_cases.FA_calculator import FactorAnalysis
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors

import numpy as np
import json

# Считывание входных данных из файлов jsonnoc
with open("250_GW_fact.json", "r") as f:
    input_data_fact = json.load(f)

with open("250_GW_plan.json", "r") as f:
    input_data_plan = json.load(f)

# #  Варьирование кол-ва трещин ГРП
# for i in range(2, 9):
#     input_data_fact["unit"]["layer_prop"]["mgrp_prop"]["grp_count"] = i
#     # print(input_data_fact)
#     factor_analysis = FactorAnalysis(input_data_fact, input_data_plan, 100)
#     q_fact = factor_analysis.calc_flow_rate(input_data_fact)
#     print(q_fact)

# #  Варьирование полудлины трещин ГРП (уменьшение полудлины на определенный % при фиксированном кол-ве стадий)
# xf0 = input_data_fact["unit"]["layer_prop"]["grp_prop"]["xf"]
# for i in range(0, 64, 4):
#     input_data_fact["unit"]["layer_prop"]["grp_prop"]["xf"] = xf0 * (1.1 - i / 100)
#     factor_analysis = FactorAnalysis(input_data_fact, input_data_plan, 100)
#     q_fact = factor_analysis.calc_flow_rate(input_data_fact)
#     print(f"{q_fact}")

# Факторный анализ
factor_analysis = FactorAnalysis(input_data_fact, input_data_plan, 100)
q_fact = factor_analysis.calc_flow_rate(input_data_fact)
q_plan = factor_analysis.calc_flow_rate(input_data_plan)
print(f" Фактический дебит нефти: {q_fact}", "\n",
      f"Плановый дебит нефти: {q_plan}")
dict_factors = factor_analysis.calc_factors()
print(f" Значения факторов: {dict_factors}", "\n",
      f"Сумма факторов: {sum(value for value in dict_factors.values() if isinstance(value, (int, float)))}", "\n",
      f"Разница дебитов: {q_plan - q_fact}", "\n",
      f"Абсолютная погрешность ФА: {abs(q_plan - q_fact -
                                        sum(value for value in dict_factors.values()
                                            if isinstance(value, (int, float))))}")
#
# # Сортируем по абсолютным значениям
# sorted_data = dict(sorted(dict_factors.items(), key=lambda item: abs(item[1])))
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
# ax.set_xlabel("Значение фактора")
# ax.set_title("Tornado Chart факторный анализ дебита")
# ax.axvline(0, color="black", linewidth=0.5)  # Центральная ось
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
#
# plt.show()
