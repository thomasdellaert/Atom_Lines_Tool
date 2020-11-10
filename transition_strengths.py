import sympy.physics.wigner as w
import py3nj as nj
import math

"""This code only currently works for E1 transitions. The math is cribbed from
 https://demonstrations.wolfram.com/TransitionStrengthsOfAlkaliMetalAtoms/"""


def rel_transition_strength_sym(I, q, J1, F0, M0, F1, M1):
    return float(3 * (2 * J1 + 1) * (2 * F0 + 1) * (2 * F1 + 1) *
                 w.wigner_6j(1, 0, 1, 0.5, J1, 0.5, prec=8) ** 2 *
                 w.wigner_6j(J1, 0.5, 1, F0, F1, I, prec=8) ** 2 *
                 w.wigner_3j(F1, 1, F0, -M1, q, M0) ** 2)


def tkq_LS_transition_strength_sym(I, k, q, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1):
    if S0 == S1:
        return float((2 * J0 + 1) * (2 * J1 + 1) * (2 * F0 + 1) * (2 * F1 + 2) *
                     w.wigner_6j(L0, L1, k, J1, J0, S0, prec=8) ** 2 *
                     w.wigner_6j(J0, J1, k, F1, F0, I, prec=8) ** 2 *
                     w.wigner_3j(F1, k, F0, -M1, q, M0) ** 2)
    else:
        return 0.0


def rel_transition_strength(I, q, J1, F0, M0, F1, M1):
    return float(3 * (2 * J1 + 1) * (2 * F0 + 1) * (2 * F1 + 1) *
                 nj.wigner6j(2, 0,           2,
                             1, int(J1 * 2), 1) ** 2 *
                 nj.wigner6j(int(J1 * 2), 1,           2,
                             int(F0 * 2), int(F1 * 2), int(I * 2)) ** 2 *
                 nj.wigner3j(int(F1 * 2),  2,          int(F0 * 2),
                             int(-M1 * 2), int(q * 2), int(M0 * 2)) ** 2)


def tkq_LS_transition_strength(I, k, q, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1):
    if S0 == S1:
        return float((2 * J0 + 1) * (2 * J1 + 1) * (2 * F0 + 1) * (2 * F1 + 2) *
                     nj.wigner6j(int(L0 * 2), int(L1 * 2), int(k * 2),
                                 int(J1 * 2), int(J0 * 2), int(S0 * 2)) ** 2 *
                     nj.wigner6j(int(J0 * 2), int(J1 * 2), int(k * 2),
                                 int(F1 * 2), int(F0 * 2), int(I * 2)) ** 2 *
                     nj.wigner3j(int(F1 * 2),  int(k * 2), int(F0 * 2),
                                 int(-M1 * 2), int(q * 2), int(M0 * 2)) ** 2)
    else:
        return 0.0


def E1_transition_strength(eps, I, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1):
    tot = 0
    tot += tkq_LS_transition_strength(I, 1, -1, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1) * 0.5 * math.sin(eps[2]) ** 2
    tot += tkq_LS_transition_strength(I, 1,  0, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1) * math.cos(eps[2]) ** 2
    tot += tkq_LS_transition_strength(I, 1,  1, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1) * 0.5 * math.sin(eps[2]) ** 2
    return tot


def M1_transition_strength(eps, I, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1):
    tot = 0
    if S1 == S0 and L1 == L0:
        tot += tkq_LS_transition_strength(I, 1, -1, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1) * 0.5 * (eps[0] + eps[1]) ** 2
        tot += tkq_LS_transition_strength(I, 1,  0, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1) * eps[2] ** 2
        tot += tkq_LS_transition_strength(I, 1,  1, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1) * 0.5 * (eps[0] + eps[1]) ** 2
        return tot
    else:
        return 0


def E2_transition_strength(eps, k, I, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1):
    tot = 0
    tot += tkq_LS_transition_strength(I, 2, -2, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1) * \
        (eps[0] ** 2 + eps[1] ** 2) * (k[0] ** 2 + k[1] ** 2)
    tot += tkq_LS_transition_strength(I, 2, -1, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1) * \
        (eps[2] * k[0] + eps[0] * k[2]) ** 2 + (eps[1] * k[0] + eps[0] * k[1]) ** 2
    tot += tkq_LS_transition_strength(I, 2, 0, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1) * \
        (2./3.) * (3 * k[0] * eps[0] + 3 * k[1] * eps[1] + 2 * k[2] * eps[2]) ** 2
    tot += tkq_LS_transition_strength(I, 2, 1, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1) * \
        (eps[2] * k[0] + eps[0] * k[2]) ** 2 + (eps[1] * k[0] + eps[0] * k[1]) ** 2
    tot += tkq_LS_transition_strength(I, 2, 2, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1) * \
        (eps[0] ** 2 + eps[1] ** 2) * (k[0] ** 2 + k[1] ** 2)
    return tot

# TODO: Implement spatial/non-polarized averages for these transition strengths, as they're currently dependent on polarization
#  and k-vector

if __name__ == '__main__':
    for (Mg, Me) in [(-3, -4), (-2, -3), (-3, -3), (-1, -2), (-2, -2), (0, -1), (-1, -1), (1, 0), (-2, -1), (2, 1), (-1, 0),
                     (3, 2), (0, 1), (1, 2), (2, 3), (3, 4)]:
        s = float(rel_transition_strength(2.5, Me - Mg, 1.5, 3, Mg, 4, Me))
        print "{0:.6f} {1:} {2:}".format(s, Mg, Me)
