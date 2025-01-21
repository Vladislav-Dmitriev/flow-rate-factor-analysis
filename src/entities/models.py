from pydantic import BaseModel, Field
from typing import Optional


class WellboreProperty(BaseModel):
    """
    Параметры скважины
    """
    wellbore_type: str = Field('', title='Тип скважины: вертикальная, горизонтальная или многозабойная', examples=['vertical, horizontal', 'multilateral'])
    wellbore_r: float = Field(0.1, title='Радиус скважины, м', gt=0)
    horizontal_wellbore_length: Optional[int] = Field(title='Длина горизонтального ствола', ge=0)
    horizontal_wellbore_perf_ratio: Optional[float] = Field(None, title='Доля перфораций горизонтальной скважины', ge=0, le=1)
    horizontal_perf_count: Optional[int] = Field(None, title='Количество перфорированных интервалов', ge=0)
    permeability: Optional[float] = Field(None, title='Проницаемость скважины', ge=0)


class Grp_Prop(BaseModel):
    """
    Параметры ГРП
    """
    hf: float = Field(None, title='Полудлина трещины, м', ge=0)
    kf: float = Field(None, title='Проницаемость трещины, дарси', ge=0)
    wellbore_wf: float = Field(None, title='Ширина трещины у ствола скважины, мм', ge=0)
    res_wf: float = Field(None, title='Ширина трещины на конце, мм', ge=0)
    grade: int = Field(None, title='Степень аппроксимационной функции', gt=0)
    skin_border: float = Field(None, title='Скин на стенке трещины', ge=0)
    skin_ch: float = Field(None, title='Скин за счет схождения (Choke Skin)', ge=0)
    fracture_grow_t: float = Field(None, title='Время роста трещины, ч', ge=0)


class Mgrp_Prop(BaseModel):
    """
    Параметры МГРП
    """
    grp_count: int = Field(None, title='Количество ГРП', gt=0)
    f_direction: str = Field('', title='Направление трещин', examples=['transverse', 'longitudal'])


class Mlt_Prop(BaseModel):
    """
    Параметры многозабойной скважины
    """
    l_lateral: int = Field(None, title='Длина бокового ствола', gt=0)
    n_lateral: int = Field(None, title='Количество боковых стволов', gt=0)
    psi_lateral: float = Field(None, title='Угол наклона боковых стволов', gt=0)


class LayerProperty(BaseModel):
    """
    Параметры пласта
    """
    permeability: float = Field(None, title='Проницаемость пласта, мД', gt=0)
    water_cut: float = Field(None, title='Обводненность, %', ge=0, le=100)
    p_b: float = Field(None, title='Давление насыщения, атм', gt=0, alias='p_bubble')
    kv_kh_ratio: float = Field(None, title='Отношение вертикальной к горизонтальной проницаемости', ge=0, le=1)
    compressibility: float = Field(None, title='Общая сжимаемость', gt=0)
    porosity: float = Field(None, title='Пористость', ge=0, le=1)
    p_res_init: float = Field(None, title='Начальное пластовое давление', gt=0)
    viscosity_oil: float = Field(None, title='Вязкость нефти', gt=0)
    b_oil: float = Field(None, title='Объемный коэффициент', gt=0)
    res_model_type: str = Field('', title='Модель пласта', examples=['Homogeneous', 'DP PSS', 'DP Slab', 'composite'])
    f_compressibility: float = Field(None, title='Сжимаемость трещин', ge=0)
    f_porosity: float = Field(None, title='Пористость трещин', ge=0)
    lambda_ratio: float = Field(None, title='Коэффициент пропускания', alias='lambda', ge=0)
    internal_r: Optional[float] = Field(None, title='Внутренний радиус (для составной модели пласта)', ge=0)
    kmu_in_out_ratio: float = Field(None, title='(k/mu)_in / (k/mu)_out', ge=0)
    kmuphict_in_out_ratio: float = Field(None, title='(k/mu_phi_ct)_in / (k/mu_phi_ct)_out', ge=0)
    grp_flag: bool = Field(default=False, title='Параметр, обозначающий, был ли ГРП на скважине')
    grp_prop: Optional[Grp_Prop] = Field(None, title='Параметры ГРП')
    mgrp_flag: bool = Field(default=False, title='Параметр, обозначающий, был ли МГРП на скважине')
    mgrp_prop: Optional[Mgrp_Prop] = Field(None, title='Параметры МГРП')
    multilateral_prop: Optional[Mlt_Prop] = Field(None, title='Параметры многозабойной скважины')
    xe: int = Field(None, title='Длина прямоугольника, м (в напр. Запад - Восток)', gt=0)
    ye: int = Field(None, title='Ширина прямоугольника, м (в напр. Север - Юг)', gt=0)
    lc_ratio: float = Field(None, title='Доля от длины до центра', ge=0, le=1)
    wc_rectangle_ratio: float = Field(None, title='Доля от ширины до центра', ge=0, le=1)
    extend_reflections_mode: Optional[bool] = Field(default=False, title='Активировать расширенный режим отражений')


class WellboreDefinition(BaseModel):
    """
    Описание скважины
    """
    skin: float = Field(None, title='Скин-фактор')
    h_eff: float = Field(None, title='Эффективная толщина пласта, м', gt=0)
    vertical_offset: float = Field(None, title='Вертикальное положение (в долях от мощности) (для горизонтальных)', ge=0, le=1)
    perfres_ratio: float = Field(None, title='Величина вскрытия пласта (в долях от мощности) (для всех)', ge=0, le=1)
    wellbore_prop: WellboreProperty = Field(title="Параметры скважины")
    layer_prop: LayerProperty = Field(title="Параметры пласта")


class BaseTask(BaseModel):
    """
    Описание базовой задачи
    """
    calc_type: str = Field("", title="Способ расчета", examples=['optimal', 'segmentation', 'desuperposition'])
    segment_count: int = Field(None, title='Число сегментов', gt=0, multiple_of=10)
    border_type_x: str = Field("", title="Тип границы по x", examples=['n', 'c'])
    border_type_y: str = Field("", title="Тип границы по y", examples=['n', 'c'])
    border_type_z_up: str = Field("", title="Тип верхней границы по z", examples=['n', 'c'])
    border_type_z_down: str = Field("", title="Тип нижней границы по z", examples=['n', 'c'])


class TargetValues(BaseModel):
    """
    Целевые значения
    """
    p_bhp: float = Field([], title="Целевое забойное давление", gt=0)
    q_liq: float = Field("", title="Целевой дебит жидкости, m3/день", gt=0)


class Target(BaseModel):
    """
    Описание задачи расчета по дебиту/по забойному
    """
    cumulative_time: float = Field(None, title="Накопленное время работы, в сутках", ge=0)
    target_values: TargetValues = Field(title='Целевые значения')


class LaplaceAccuracyOptions(BaseModel):
    """
    Параметры точности вычислений Лапласа
    """
    tiny: Optional[float] = Field(default=0.00000001, title='Параметр управления сходимостью рядов')
    tiny_2: Optional[float] = Field(default=1E-10, title='Параметр управления сходимостью рядов')
    large_s: Optional[float] = Field(default=1E-300, title='Параметр улучшения сходимости')
    small: Optional[float] = Field(default=1E-300, title='Параметр аппроксимации времени')
    number_of_lapl_coeff: Optional[int] = Field(default=10, title=' Количество коэффициентов Лапласа')
    part_sum_num: Optional[int] = Field(default=5, title='Минимальное приращение')
    max_it: Optional[int] = Field(default=2000, title='Максимальное количество итераций при суммировании рядов')
    small_bessel_arg: Optional[float] = Field(default=0.0001, title='Промежуток времени для функции Бесселя')
    part_sum_num_b2_2: Optional[int] = Field(default=10, title='Минимальное приращение')
    part_sum_num_f_4: Optional[int] = Field(default=10000, title='Минимальное приращение')


class TaskTarget(BaseModel):
    """
    Класс для расчета дебита жидкости / давления
    """
    base_prop: BaseTask = Field(title='Описание базовой задачи')
    target: Target = Field(title='Описание задачи расчета по дебиту/по забойному')
    unit: WellboreDefinition = Field(title='Описание скважины')
    options: Optional[LaplaceAccuracyOptions] = Field(None, title='Параметры точности вычислений Лапласа')
