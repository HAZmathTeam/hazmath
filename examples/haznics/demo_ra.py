"""
\file examples/haznics/demo_ra.py
Created by Ana Budisa on 2021-07-28.

We look for rational approximation of the function

     (alpha * x^s + beta * x^t)^-1

on interval [0, bound], i.e. we want to find scalars (c_k, d_k), k = 0, 1, ... such that

     (alpha * x^s + beta * x^t)^-1 \approx c_0 + sum( c_k (x - d_k)^-1 )

We assume -1 <= s,t <= 1, alpha, beta >= 0, and bound >= 1.
We use the algorithms implemented in hazmath library to compute (c_k, d_k).
In this example, we use AAA algorithm from

Y. Nakatsukasa, O. Sète, and L. N. Trefethen, The AAA Algorithm for Rational Approximation.
SIAM Journal on Scientific Computing, 40 (2018), pp. A1494-A1522.

The tolerance for the approximation is set to 1E-12.
"""
import haznics


def rational_approx(s=-0.5, t=0.5, alpha=1., beta=1., scaling_alpha=1., scaling_beta=1.):
    # Approximation of (alpha*x^s + beta*x^t)^(-1) on interval [0, 1]
    res = haznics.compute_ra_aaa(s, t, alpha, beta, scaling_alpha, scaling_beta)

    res = res.to_ndarray()
    k = int((res.size - 1)/2)
    poles = res[:k]
    residues = res[k+1:]

    return poles, residues


if __name__ == "__main__":
    s = -0.5
    t = 0.0
    alpha = 1E-2
    beta = 1E-6

    print("Rational approximation without scalings")
    poles, residues = rational_approx(s, t, alpha, beta)

    print("---- Number of poles: ", poles.size)
    print("---- Poles:")
    print(poles)
    print("---- Residues:")
    print(residues)

    scaling_alpha = 10
    scaling_beta = 0.1

    print("Rational approximation with scalings")
    poles, residues = rational_approx(s, t, alpha, beta, scaling_alpha, scaling_beta)

    print("---- Number of poles: ", poles.size)
    print("---- Poles:")
    print(poles)
    print("---- Residues:")
    print(residues)

    # poles and residues can be saved for possible input in hazmath
    # np.savetxt('aaa_data/poles_s_05_t_05.out', poles)
    # np.savetxt('aaa_data/residues_s_05_t_05.out', residues)

