#!/usr/bin/env python
"""
单元测试：验证 GCMCAcceptance 的详细对称公式实现。
重点覆盖 cavity bias 与 Λ³ 因子，确保 log-space 计算与解析公式一致。
"""

from __future__ import annotations

import math

import pygcmc

KB = 8.314e-3  # kJ/(mol·K)


def _clamp_probability(log_ratio: float) -> float:
    """与 C++ safeExp 保持一致的夹逼策略。"""
    if log_ratio >= 0.0:
        return 1.0
    return math.exp(max(log_ratio, -700.0))


def _expected_insertion(
    *,
    beta: float,
    mu: float,
    volume_nm3: float,
    cavity_fraction: float,
    n_before: int,
    lambda_nm: float,
    rosenbluth_weight: float,
    proposal_log_ratio: float,
    delta_e: float,
) -> float:
    log_ratio = (
        proposal_log_ratio
        - beta * delta_e
        + beta * mu
        + math.log(volume_nm3)
        + math.log(cavity_fraction)
        - math.log(n_before + 1.0)
        + math.log(rosenbluth_weight)
        - 3.0 * math.log(lambda_nm)
    )
    return _clamp_probability(log_ratio)


def _expected_deletion(
    *,
    beta: float,
    mu: float,
    volume_nm3: float,
    cavity_fraction: float,
    n_before: int,
    lambda_nm: float,
    rosenbluth_weight: float,
    proposal_log_ratio: float,
    delta_e: float,
) -> float:
    log_ratio = (
        proposal_log_ratio
        - beta * delta_e
        - beta * mu
        + math.log(n_before)
        - (math.log(volume_nm3) + math.log(cavity_fraction))
        - math.log(rosenbluth_weight)
        + 3.0 * math.log(lambda_nm)
    )
    return _clamp_probability(log_ratio)


def test_insertion_probability_detailed_matches_formula() -> None:
    acc = pygcmc.GCMCAcceptance()
    temperature = 298.15
    volume_nm3 = 42.0
    type_id = 0
    mu = -5.0  # kJ/mol
    beta = 1.0 / (KB * temperature)
    cavity_fraction = 0.35
    lambda_nm = 0.6
    n_before = 4
    delta_e = 1.2
    rosenbluth_weight = 2.5  # W_new / K
    cbmc_trials = 6
    proposal_ratio = 1.8  # p_delete / p_insert
    proposal_log_ratio = math.log(proposal_ratio)

    acc.setTemperature(temperature)
    acc.setVolume(volume_nm3)
    acc.setChemicalPotential(type_id, mu)
    acc.setThermalLambda(type_id, lambda_nm)

    expected = _expected_insertion(
        beta=beta,
        mu=mu,
        volume_nm3=volume_nm3,
        cavity_fraction=cavity_fraction,
        n_before=n_before,
        lambda_nm=lambda_nm,
        rosenbluth_weight=rosenbluth_weight,
        proposal_log_ratio=proposal_log_ratio,
        delta_e=delta_e,
    )

    expected_log_ratio = (
        proposal_log_ratio
        - beta * delta_e
        + beta * mu
        + math.log(volume_nm3)
        + math.log(cavity_fraction)
        - math.log(n_before + 1.0)
        + math.log(rosenbluth_weight)
        - 3.0 * math.log(lambda_nm)
    )

    result_cpp = acc.calculate_insertion_probability_detailed(
        type_id,
        n_before,
        delta_e,
        cavity_fraction,
        lambda_nm,
        rosenbluth_weight,
        cbmc_trials,
        proposal_log_ratio,
    )

    assert math.isclose(
        result_cpp["probability"], expected, rel_tol=1e-12, abs_tol=1e-12
    )

    assert math.isclose(
        result_cpp["logRatio"], expected_log_ratio, rel_tol=1e-12, abs_tol=1e-12
    )


def test_deletion_probability_detailed_matches_formula() -> None:
    acc = pygcmc.GCMCAcceptance()
    temperature = 310.0
    volume_nm3 = 30.0
    type_id = 1
    mu = -3.2
    beta = 1.0 / (KB * temperature)
    cavity_fraction = 0.42
    lambda_nm = 0.75
    n_before = 5
    delta_e = -0.8  # 删除后能量 - 删除前能量
    rosenbluth_weight = 1.7  # W_old / K
    cbmc_trials = 4
    proposal_ratio = 1.5  # p_delete / p_insert
    proposal_log_ratio = math.log(proposal_ratio)

    acc.setTemperature(temperature)
    acc.setVolume(volume_nm3)
    acc.setChemicalPotential(type_id, mu)
    acc.setThermalLambda(type_id, lambda_nm)

    expected = _expected_deletion(
        beta=beta,
        mu=mu,
        volume_nm3=volume_nm3,
        cavity_fraction=cavity_fraction,
        n_before=n_before,
        lambda_nm=lambda_nm,
        rosenbluth_weight=rosenbluth_weight,
        proposal_log_ratio=proposal_log_ratio,
        delta_e=delta_e,
    )

    expected_log_ratio = (
        proposal_log_ratio
        - beta * delta_e
        - beta * mu
        + math.log(n_before)
        - (math.log(volume_nm3) + math.log(cavity_fraction))
        - math.log(rosenbluth_weight)
        + 3.0 * math.log(lambda_nm)
    )

    result_cpp = acc.calculate_deletion_probability_detailed(
        type_id,
        n_before,
        delta_e,
        cavity_fraction,
        lambda_nm,
        rosenbluth_weight,
        cbmc_trials,
        proposal_log_ratio,
    )

    assert math.isclose(
        result_cpp["probability"], expected, rel_tol=1e-12, abs_tol=1e-12
    )

    assert math.isclose(
        result_cpp["logRatio"], expected_log_ratio, rel_tol=1e-12, abs_tol=1e-12
    )


def test_deletion_log_ratio_decreases_with_positive_delta_energy() -> None:
    acc = pygcmc.GCMCAcceptance()
    temperature = 300.0
    beta = 1.0 / (KB * temperature)

    type_id = 0
    mu = 0.0  # keep other terms fixed
    volume_nm3 = 100.0
    cavity_fraction = 1.0
    lambda_nm = 1.0
    rosenbluth_weight = 1.0
    cbmc_trials = 1
    proposal_log_ratio = 0.0
    n_before = 5

    acc.setTemperature(temperature)
    acc.setVolume(volume_nm3)
    acc.setChemicalPotential(type_id, mu)
    acc.setThermalLambda(type_id, lambda_nm)

    res_low = acc.calculate_deletion_probability_detailed(
        type_id,
        n_before,
        0.0,
        cavity_fraction,
        lambda_nm,
        rosenbluth_weight,
        cbmc_trials,
        proposal_log_ratio,
    )

    res_high = acc.calculate_deletion_probability_detailed(
        type_id,
        n_before,
        5.0,  # kJ/mol
        cavity_fraction,
        lambda_nm,
        rosenbluth_weight,
        cbmc_trials,
        proposal_log_ratio,
    )

    assert res_high["logRatio"] < res_low["logRatio"]
    assert math.isclose(
        res_high["logRatio"] - res_low["logRatio"], -beta * 5.0, rel_tol=1e-12, abs_tol=1e-12
    )


def test_activity_semantics_do_not_double_count_lambda_in_detailed_terms() -> None:
    """
    P1 回归：当使用 setActivity(z) 直接指定 grand-canonical activity（z=exp(βμ)/Λ³）时，
    detailed 接受率中的 Λ³ 因子不应再被额外扣一次（否则会退化成 z/Λ³）。

    这个测试刻意设置 thermalLambda!=1，并在调用 detailed 接口时传入同样的 lambdaNm，
    以模拟 engine 路径（engine 会把 getThermalLambda(typeId) 传入 terms）。
    """
    type_id = 0
    temperature = 298.15
    beta = 1.0 / (KB * temperature)

    volume_nm3 = 8.0
    cavity_fraction = 0.25
    lambda_nm = 0.6
    mu = -6.0  # kJ/mol

    z = math.exp(beta * mu) / (lambda_nm**3)

    n_before = 6
    delta_e = 1.3
    rosenbluth_weight = 1.4
    cbmc_trials = 5
    proposal_log_ratio = math.log(0.8)

    acc = pygcmc.GCMCAcceptance()
    acc.setTemperature(temperature)
    acc.setVolume(volume_nm3)
    acc.setActivity(type_id, z)
    acc.setThermalLambda(type_id, lambda_nm)

    expected_ins_log_ratio = (
        proposal_log_ratio
        - beta * delta_e
        + math.log(z)
        + math.log(volume_nm3)
        + math.log(cavity_fraction)
        - math.log(n_before + 1.0)
        + math.log(rosenbluth_weight)
    )

    ins = acc.calculate_insertion_probability_detailed(
        typeId=type_id,
        currentNumber=n_before,
        deltaE=delta_e,
        cavityFraction=cavity_fraction,
        lambdaNm=lambda_nm,
        rosenbluthWeight=rosenbluth_weight,
        cbmcTrials=cbmc_trials,
        proposalLogRatio=proposal_log_ratio,
    )
    assert math.isclose(
        ins["logRatio"], expected_ins_log_ratio, rel_tol=1e-12, abs_tol=1e-12
    )

    expected_del_log_ratio = (
        proposal_log_ratio
        - beta * delta_e
        - math.log(z)
        + math.log(float(n_before))
        - (math.log(volume_nm3) + math.log(cavity_fraction))
        - math.log(rosenbluth_weight)
    )

    dele = acc.calculate_deletion_probability_detailed(
        typeId=type_id,
        currentNumber=n_before,
        deltaE=delta_e,
        cavityFraction=cavity_fraction,
        lambdaNm=lambda_nm,
        rosenbluthWeight=rosenbluth_weight,
        cbmcTrials=cbmc_trials,
        proposalLogRatio=proposal_log_ratio,
    )
    assert math.isclose(
        dele["logRatio"], expected_del_log_ratio, rel_tol=1e-12, abs_tol=1e-12
    )
