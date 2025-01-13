from tests.test_calc import test_calc_flow_rate, test_calc_p_bhp


def test_calc_flow_benchmark(benchmark):
    benchmark(test_calc_flow_rate)


def test_calc_p_bhp_benchmark(benchmark):
    benchmark(test_calc_p_bhp)
