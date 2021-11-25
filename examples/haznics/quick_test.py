"""
Computes the rational approximation of (alpha*x^s + beta*x^t)^(-1) on interval [0, bound]:
    * -1 <= s,t <= 1
    * upper bound of x, should be an upper bound of lambda_max(A), e.g. l_infty norm of A
"""
import numpy as np
from baryrat import aaa
import haznics


def rational_approx_haznics(x, f, AAA_tol):
    # Approximation of (alpha*x^s + beta*x^t)^(-1) on interval [0, 1]
    # import pdb; pdb.set_trace()
    res = haznics.ra_aaa(x, f, AAA_tol)

    res = res.to_ndarray()
    k = int((res.size + 2)/4)

    poles_r = res[:k-1]
    poles_i = res[(k-1):(2*k-2)]
    residues_r = res[(2*k-2):(3*k-2)]
    residues_i = res[(3*k-2):]

    if np.any(poles_i > 1E-12) or np.any(residues_i > 1E-12):
        poles = poles_r + 1j * poles_i
        residues = residues_r + 1j * residues_i
    else:
        poles = poles_r
        residues = residues_r

    return poles, residues


if __name__ == "__main__":
    # params
    s = -0.2
    t = 0.6
    alpha = 1.0
    beta = 1.0
    coeff = 1./alpha if alpha > beta else 1/beta

    bound = 1.0
    tol = 1E-12
    npoints = (8*1024)+1

    # interpolation points
    x = np.linspace(0.0001, bound, npoints)

    # interpolation values at the points
    F = np.power((alpha * np.power(x, s) + beta * np.power(x, t)), -1)

    print("\nRational approx of a function f(x) = (a * x^s + b * x^t)^-1 on [0, 1]")
    print("s = %.1f, t = %.1f, alpha = %.1f, beta = %.1f\n" % (s, t, alpha, beta))

    # compute poles, res with haznics
    pol_haz, res_haz = rational_approx_haznics(x, F, tol)

    # compute poles, res with baryrat
    r = aaa(x, F, tol=tol)
    pol_br, res_br = r.polres()

    # compute the constant term c0
    c0 = np.power((alpha * np.power(1.0, s) + beta * np.power(1.0, t)), -1) - np.sum(res_br / (1.0 - pol_br))
    res_br = np.append(c0, res_br)

    # different points
    y = np.linspace(0.0001, bound, 6001)
    FF = np.power((alpha * np.power(y, s) + beta * np.power(y, t)), -1)

    rff = coeff * res_haz[0] * np.ones_like(FF)
    for pk, rk in zip(pol_haz, res_haz[1:]):
        rff[:] += coeff * rk / (y - pk)

    # error in interpolation points
    er = F - r(x)
    print(" BR  ---- Rational approx error in interp points: ", np.linalg.norm(er, np.inf) / np.linalg.norm(F, np.inf))

    # er_haz = F - rf

    # error in different points
    err = FF - r(y)
    err_haz = FF - rff

    print("\n HAZ ---- Number of poles: ", pol_haz.size)
    print(" HAZ ---- Poles: \n", pol_haz)
    print(" HAZ ---- Residues: \n", res_haz)

    print()
    print(" BR ---- Number of poles: ", pol_br.size)
    print(" BR ---- Poles: \n", pol_br)
    print(" BR ---- Residues: \n", res_br)

    # print(" HAZ ---- Rational approx error in interp points: ", np.linalg.norm(er_haz, np.inf) / np.linalg.norm(F, np.inf))
    print(" HAZ ---- Rational approx error in different points: ", np.linalg.norm(err_haz, np.inf) / np.linalg.norm(FF, np.inf))
    print(" BR  ---- Rational approx error in different points: ", np.linalg.norm(err, np.inf) / np.linalg.norm(FF, np.inf))


