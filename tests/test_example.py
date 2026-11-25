"""Numerical tests: run the bundled example and compare with reference values.

The reference values are the results of the pure Python implementation
(openCOSMO-RS_py, tests/test_cpp_port.py) for the same COSMO files and options.
"""
import os
import runpy

import numpy as np
import pytest

import openCOSMORS  # noqa: F401  (the example imports it as well)

EXAMPLE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bindings", "run_example.py")


@pytest.fixture(scope="module")
def calculations():
    cwd = os.getcwd()  # the example changes the working directory
    try:
        namespace = runpy.run_path(EXAMPLE, run_name="__main__")
    finally:
        os.chdir(cwd)
    return namespace["calculations"]


def test_residual_ln_gamma_cosmo_reference_state(calculations):
    np.testing.assert_allclose(
        calculations[0]["ln_gamma_x_SR_residual_calc"], [[4.781831, 2.211788, 0.955943]], atol=1e-4
    )


def test_residual_ln_gamma_pure_component_reference_state(calculations):
    np.testing.assert_allclose(
        calculations[2]["ln_gamma_x_SR_residual_calc"], [[2.325362, -0.037878, 1.313109]], atol=1e-4
    )


def test_average_surface_energies_are_in_joule_per_mol(calculations):
    energies = calculations[0]["average_surface_energies"][0, 0, :, :].sum(1)
    np.testing.assert_allclose(energies, [7759.787, -1692.607, -6864.229], rtol=1e-4, atol=0.5)


def test_partial_molar_energies(calculations):
    energies = calculations[0]["partial_molar_energies"][0, 0, :]
    np.testing.assert_allclose(energies, [7502.846, -1166.276, -7570.112], rtol=1e-4, atol=0.5)
