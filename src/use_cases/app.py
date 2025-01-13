import time
from fastapi import FastAPI
from fastapi.responses import JSONResponse
from logger import logger

from src.entities.models import TaskTarget
import src.use_cases.calc as calc

app = FastAPI(
    title="NGT UNIFLOC - kprod_services API",
    description="API для расчета коэффициента продуктивности",
)


@app.get("/")
async def root():
    """Проверка подключения.

    :return: Результат подключения
    :rtype: fastapi.responses.JSONResponse
    """

    start_time = time.time()

    res = {"message": "Hello! It's NGT Unifloc | kprod_services API."}

    logger.debug(f"Ping kprod_services | Time: {time.time() - start_time}")

    return JSONResponse(res)


@app.post("/calc_flow_rate/")
def send_to_calc_flow_rate(input_data: TaskTarget) -> JSONResponse:
    """Отправка параметров для расчета дебита жидкости

    :param input_data: Входные данные для расчета параметров
    :type input_data: src.entities.models.TaskTarget

    :return: Дебит жидкости (в сутки)
    :rtype: fastapi.responses.JSONResponse
    """

    start_time = time.time()

    response = calc.send_to_calc_flow_rate(input_data)

    logger.debug(f"Send to calc_flow_rate | Time: {time.time() - start_time}")

    return JSONResponse(response)


@app.post("/calc_p_bhp/")
def send_to_calc_p_bhp(input_data: TaskTarget) -> JSONResponse:
    """Отправка параметров для расчета давления

    :param input_data: Входные данные для расчета параметров
    :type input_data: src.entities.models.TaskTargetPbhp

    :return: Результат расчета параметров
    :rtype: fastapi.responses.JSONResponse
    """

    start_time = time.time()

    response = calc.send_to_calc_p_bhp(input_data)

    logger.debug(f"Send to calc_p_bhp | Time: {time.time() - start_time}")

    return JSONResponse(response)


if __name__ == "__main__":
    pass
