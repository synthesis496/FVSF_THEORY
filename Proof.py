#Copyright (c) 2026 Chutiphong Bunloed. All Rights Reserved.
#DOI 10.5281/zenodo.23087025
import math

phi = (1 + math.sqrt(5)) / 2

D_total = 248
D_SM = 36
D_base = 16
D_SO32 = 496

alpha_inv_fvs = D_total / phi - D_base - phi**-3
dm_fraction_fvs = (D_total - D_SM) / D_total
e_recoil_fvs = 510.99895000 / phi**(3/2)
sin2_theta_w_fvs = phi**-3
n_s_fvs = 1 - 9 / (D_total - D_SM)
exponent_lambda = 2 * D_total + int((D_total - D_SM) / D_total * 100)
lambda_fvs = phi**-exponent_lambda
r_fvs = phi**-12
mu_fvs = 512 * phi**2 + D_SO32
G_inv_norm = (D_total - D_SM) / D_total * 160 + 4 - phi**-3

assert exponent_lambda == 581
assert abs(alpha_inv_fvs - 137.0363612) < 1e-6
assert abs(dm_fraction_fvs - 0.8548387) < 1e-6
assert abs(e_recoil_fvs - 248.276) < 0.01
assert abs(sin2_theta_w_fvs - 0.2360679) < 1e-6
assert abs(n_s_fvs - 0.9575471) < 1e-6
assert abs(lambda_fvs - 3.78e-122) < 0.01e-122
assert abs(mu_fvs - 1836.4334) < 0.01
assert abs(G_inv_norm - 140.538) < 0.001

result = {
    "alpha_inv": alpha_inv_fvs,
    "dm_fraction": dm_fraction_fvs,
    "e_recoil_keV": e_recoil_fvs,
    "sin2_theta_w": sin2_theta_w_fvs,
    "n_s": n_s_fvs,
    "lambda_exponent": exponent_lambda,
    "lambda_fvs": lambda_fvs,
    "r_fvs": r_fvs,
    "mu_fvs": mu_fvs,
    "G_inv_norm": G_inv_norm
}

result
