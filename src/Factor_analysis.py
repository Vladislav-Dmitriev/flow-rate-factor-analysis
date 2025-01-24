input_data_plan = {
    "base_prop": {
        "calc_type": 'optimal',  # тип расчета (optimal упрощенная модель, desuperposition с эффектом суперпозиции
        # давления, segmentation сегментация скважины на интервалы)
        "segment_count": 20,  # кол-во сегментов разбиения скважины
        "border_type_x": "n",  # тип граничного условия (c - непроницаемая граница, n - свободный поток через границы)
        "border_type_y": "n",  # тип граничного условия
        "border_type_z_up": "n",  # тип граничного условия
        "border_type_z_down": "n",  # тип граничного условия
    },
    "target": {
        "cumulative_time": 2400,  # накопленное время работы, в сутках
        "target_values": {
            "q_liq": 50,  # целевое значение дебита
            "p_bhp": 150 # целевое значение забойного давления
        },
    },
    "unit": {
        "skin": 0,  # скин-фактор
        "h_eff": 12.5,  # Эффективная толщина пласта
        "vertical_offset": 0.5,  #
        "perfres_ratio": 1,  # Вертикальное положение (в долях от мощности) (для горизонтальных)
        "wellbore_prop": {
            "wellbore_type": "horizontal",  # Тип скважины: вертикальная, горизонтальная или многозабойная
            "wellbore_r": 0.1,  # Радиус скважины, м
            "horizontal_wellbore_length": 624,  # Длина горизонтального ствола
            "horizontal_wellbore_perf_ratio": 1,  # Доля перфораций горизонтальной скважины
            "horizontal_perf_count": 5,  # Количество перфорированных интервалов
            "permeability": 100000000  # Проницаемость скважины
        },
        "layer_prop": {
            "permeability": 12,  # Проницаемость пласта, мД
            "water_cut": 12,  # Обводненность, %
            "p_bubble": 100,  # Давление насыщения, атм
            "kv_kh_ratio": 0.1,  # Отношение вертикальной к горизонтальной проницаемости
            "compressibility": 0.00015,  # Общая сжимаемость
            "porosity": 0.2,  # Пористость
            "p_res_init": 260,  # Начальное пластовое давление
            "viscosity_oil": 1,  # Вязкость нефти
            "b_oil": 1.2,  # Объемный коэффициент
            "res_model_type": "Homogeneous",  # Модель пласта
            "f_compressibility": 0.0001,  # Сжимаемость трещин
            "f_porosity": 0.001,  # Пористость трещин
            "lambda": 0.00001,  # Коэффициент пропускания
            "internal_r": None,  # Внутренний радиус (для составной модели пласта)
            "kmu_in_out_ratio": 2,  # (k/mu)_in / (k/mu)_out
            "kmuphict_in_out_ratio": 2,  # (k/mu_phi_ct)_in / (k/mu_phi_ct)_out
            "grp_flag": False,  # Параметр, обозначающий, был ли ГРП на скважине
            "grp_prop": {
                "hf": 60,  # Полудлина трещины, м
                "kf": 1500,  # Проницаемость трещины, дарси
                "wellbore_wf": 4.5,  # Ширина трещины у ствола скважины, мм
                "res_wf": 4.31,  # Ширина трещины на конце, мм
                "grade": 1,  # Степень аппроксимационной функции
                "skin_border": 0,  # Скин на стенке трещины
                "skin_ch": 0,  # Скин за счет схождения (Choke Skin)
                "fracture_grow_t": 0,  # Время роста трещины, ч
            },
            "mgrp_flag": False,  # Параметры МГРП
            "mgrp_prop": {
                    "grp_count": 2,  # Количество ГРП
                    "f_direction": "transverse"  # Направление трещин
                },
            "xe": 1000,  # Длина прямоугольника, м (в напр. Запад - Восток)
            "ye": 1000,  # Ширина прямоугольника, м (в напр. Север - Юг)
            "lc_ratio": 0.5,  # Доля от длины до центра
            "wc_rectangle_ratio": 0.5  # Доля от ширины до центра
        }
    }
}

input_data_fact = {
    "base_prop": {
        "calc_type": 'optimal',  # тип расчета (optimal упрощенная модель, desuperposition с эффектом суперпозиции
        # давления, segmentation сегментация скважины на интервалы)
        "segment_count": 20,  # кол-во сегментов разбиения скважины
        "border_type_x": "n",  # тип граничного условия (c - непроницаемая граница, n - свободный поток через границы)
        "border_type_y": "n",  # тип граничного условия
        "border_type_z_up": "n",  # тип граничного условия
        "border_type_z_down": "n",  # тип граничного условия
    },
    "target": {
        "cumulative_time": 2400,  # накопленное время работы, в сутках
        "target_values": {
            "q_liq": 50,  # целевое значение дебита
            "p_bhp": 150 # целевое значение забойного давления
        },
    },
    "unit": {
        "skin": 0,  # скин-фактор
        "h_eff": 12,  # Эффективная толщина пласта
        "vertical_offset": 0.5,  #
        "perfres_ratio": 1,  # Вертикальное положение (в долях от мощности) (для горизонтальных)
        "wellbore_prop": {
            "wellbore_type": "horizontal",  # Тип скважины: вертикальная, горизонтальная или многозабойная
            "wellbore_r": 0.1,  # Радиус скважины, м
            "horizontal_wellbore_length": 624,  # Длина горизонтального ствола
            "horizontal_wellbore_perf_ratio": 1,  # Доля перфораций горизонтальной скважины
            "horizontal_perf_count": 5,  # Количество перфорированных интервалов
            "permeability": 100000000  # Проницаемость скважины
        },
        "layer_prop": {
            "permeability": 11,  # Проницаемость пласта, мД
            "water_cut": 12,  # Обводненность, %
            "p_bubble": 100,  # Давление насыщения, атм
            "kv_kh_ratio": 0.1,  # Отношение вертикальной к горизонтальной проницаемости
            "compressibility": 0.00015,  # Общая сжимаемость
            "porosity": 0.2,  # Пористость
            "p_res_init": 260,  # Начальное пластовое давление
            "viscosity_oil": 1,  # Вязкость нефти
            "b_oil": 1.2,  # Объемный коэффициент
            "res_model_type": "Homogeneous",  # Модель пласта
            "f_compressibility": 0.0001,  # Сжимаемость трещин
            "f_porosity": 0.001,  # Пористость трещин
            "lambda": 0.00001,  # Коэффициент пропускания
            "internal_r": None,  # Внутренний радиус (для составной модели пласта)
            "kmu_in_out_ratio": 2,  # (k/mu)_in / (k/mu)_out
            "kmuphict_in_out_ratio": 2,  # (k/mu_phi_ct)_in / (k/mu_phi_ct)_out
            "grp_flag": False,  # Параметр, обозначающий, был ли ГРП на скважине
            "grp_prop": {
                "hf": 60,  # Полудлина трещины, м
                "kf": 1500,  # Проницаемость трещины, дарси
                "wellbore_wf": 4.5,  # Ширина трещины у ствола скважины, мм
                "res_wf": 4.31,  # Ширина трещины на конце, мм
                "grade": 1,  # Степень аппроксимационной функции
                "skin_border": 0,  # Скин на стенке трещины
                "skin_ch": 0,  # Скин за счет схождения (Choke Skin)
                "fracture_grow_t": 0,  # Время роста трещины, ч
            },
            "mgrp_flag": False,  # Параметры МГРП
            "mgrp_prop": {
                    "grp_count": 2,  # Количество ГРП
                    "f_direction": "transverse"  # Направление трещин
                },
            "xe": 1000,  # Длина прямоугольника, м (в напр. Запад - Восток)
            "ye": 1000,  # Ширина прямоугольника, м (в напр. Север - Юг)
            "lc_ratio": 0.5,  # Доля от длины до центра
            "wc_rectangle_ratio": 0.5  # Доля от ширины до центра
        }
    }
}

import pandas as pd
import re
from src.use_cases.FA_calculator import FactorAnalysis

factor_analysis = FactorAnalysis(input_data_fact, input_data_plan, 100)
q_fact = factor_analysis.calc_flow_rate(input_data_fact)
q_plan = factor_analysis.calc_flow_rate(input_data_plan)
print(f" Фактический дебит: {q_fact}", "\n",
      f"Плановый дебит: {q_plan}")
dict_factors = factor_analysis.calc_qkliq()
print(f"Значения факторов: {dict_factors}", "\n",
      f"Сумма факторов: {sum(value for value in dict_factors.values() if isinstance(value, (int, float)))}", "\n",
      f"Разница дебитов: {q_fact - q_plan}", "\n",
      f"Погрешность ФА: {q_fact - q_plan - 
                         sum(value for value in dict_factors.values() if isinstance(value, (int, float)))}")

# df = pd.read_excel("C:\\Users\\USER\\GPN\\Факторный анализ дебита\\Документы\\РБ_Факторный_Анализ_Итог.xlsm",
#                    sheet_name="Рейтинг бурения по пластам", skiprows=3,
#                    header=[0, 1, 2, 3, 4], engine="openpyxl").fillna('')
# new_columns = []
# for col in df.columns:
#     col = list(col)
#     # Добавляем элементы, которые не содержат "Unnamed" и либо содержат буквы, либо не содержат цифр
#     new_columns.append("_".join(
#         [str(i) for i in col if
#          ("Unnamed" not in str(i) and (re.search(r'[a-zA-Zа-яА-Я]', str(i))))]
#     ))
#
# df.columns = new_columns
# print(df)