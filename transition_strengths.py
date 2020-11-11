import py3nj as nj
import numpy as np


# noinspection PyTypeChecker
def tkq_LS_transition_strength(I, k, q, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1):
    if S0 == S1:
        return float((2 * J0 + 1) * (2 * J1 + 1) * (2 * F0 + 1) * (2 * F1 + 2) *
                     # nj.wigner6j(int(L0 * 2), int(L1 * 2), int(k * 2),
                     #             int(J1 * 2), int(J0 * 2), int(S0 * 2)) ** 2 *
                     nj.wigner6j(int(J0 * 2), int(J1 * 2), int(k * 2),
                                 int(F1 * 2), int(F0 * 2), int(I * 2)) ** 2 *
                     nj.wigner3j(int(F1 * 2), int(k * 2), int(F0 * 2),
                                 int(-M1 * 2), int(q * 2), int(M0 * 2)) ** 2)
    else:
        return 0.0


# The tensor math for E2 (and somewhat E1) transitions is adapted from Tony's thesis (Ransford 2020)

def E1_transition_strength_geom(eps, I, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1):
    eps = eps / np.linalg.norm(eps)
    tot = 0
    tot += tkq_LS_transition_strength(I, 1, -1, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1) * 0.5 * (eps[0] + eps[1]) ** 2
    tot += tkq_LS_transition_strength(I, 1, 0, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1) * eps[2] ** 2
    tot += tkq_LS_transition_strength(I, 1, 1, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1) * 0.5 * (eps[0] + eps[1]) ** 2
    return tot


def E1_transition_strength_avg(I, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1):
    if (L1 - L0) % 2 == 1:
        tot = 0
        tot += tkq_LS_transition_strength(I, 1, -1, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1) * (1.0 / 3.0)
        tot += tkq_LS_transition_strength(I, 1, 0, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1) * (1.0 / 3.0)
        tot += tkq_LS_transition_strength(I, 1, 1, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1) * (1.0 / 3.0)
        return tot
    else:
        return 0


def M1_transition_strength_geom(eps, I, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1):
    eps = eps / np.linalg.norm(eps)
    tot = 0
    if S1 == S0 and L1 == L0:
        tot += tkq_LS_transition_strength(I, 1, -1, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1) * 0.5 * (eps[0] + eps[1]) ** 2
        tot += tkq_LS_transition_strength(I, 1,  0, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1) * eps[2] ** 2
        tot += tkq_LS_transition_strength(I, 1,  1, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1) * 0.5 * (eps[0] + eps[1]) ** 2
        return tot
    else:
        return 0


def M1_transition_strength_avg(I, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1):
    tot = 0
    if S1 == S0 and L1 == L0:
        tot += tkq_LS_transition_strength(I, 1, -1, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1) * (1.0 / 3.0)
        tot += tkq_LS_transition_strength(I, 1,  0, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1) * (1.0 / 3.0)
        tot += tkq_LS_transition_strength(I, 1,  1, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1) * (1.0 / 3.0)
        return tot
    else:
        return 0


def E2_transition_strength_geom(eps, k, I, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1):
    eps = eps / np.linalg.norm(eps)
    k = k / np.linalg.norm(k)
    if np.dot(eps, k) != 0:
        print "k-vector and polarization are not orthogonal"
    tot = 0
    tot += tkq_LS_transition_strength(I, 2, -2, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1) * \
        (eps[0] ** 2 + eps[1] ** 2) * (k[0] ** 2 + k[1] ** 2)
    tot += tkq_LS_transition_strength(I, 2, -1, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1) * \
        (eps[2] * k[0] + eps[0] * k[2]) ** 2 + (eps[1] * k[0] + eps[0] * k[1]) ** 2
    tot += tkq_LS_transition_strength(I, 2, 0, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1) * \
        (2. / 3.) * (3 * k[0] * eps[0] + 3 * k[1] * eps[1] + 2 * k[2] * eps[2]) ** 2
    tot += tkq_LS_transition_strength(I, 2, 1, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1) * \
        (eps[2] * k[0] + eps[0] * k[2]) ** 2 + (eps[1] * k[0] + eps[0] * k[1]) ** 2
    tot += tkq_LS_transition_strength(I, 2, 2, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1) * \
        (eps[0] ** 2 + eps[1] ** 2) * (k[0] ** 2 + k[1] ** 2)
    return tot


def E2_transition_strength_avg(I, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1):
    tot = 0
    tot += tkq_LS_transition_strength(I, 2, -2, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1) * (3.0 / 29.0)
    tot += tkq_LS_transition_strength(I, 2, -1, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1) * (9.0 / 58.0)
    tot += tkq_LS_transition_strength(I, 2, 0, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1) * (14.0 / 29.0)
    tot += tkq_LS_transition_strength(I, 2, 1, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1) * (9.0 / 58.0)
    tot += tkq_LS_transition_strength(I, 2, 2, L0, S0, J0, F0, M0, L1, S1, J1, F1, M1) * (4.0 / 15.0)
    return tot

# TODO: Check the geometric factors on the averages for the E1 and M1 transitions

# TODO: Make this work outside of LS coupling. Probably involves converting
#  between LS, JJ, and JK couplings in a different .py file


if __name__ == '__main__':
    I_0 = 2.5
    L_0 = 3
    L_1 = 3
    S_0 = 0.5
    S_1 = 0.5
    J_0 = 3.5
    J_1 = 3.5
    F_0 = 4
    F_1 = 3
    G = [-4, -3, -2, -1, 0, 1, 2, 3, 4]
    E = [-3, -2, -1, 0, 1, 2, 3]
    print "M1"
    flg = False
    for mg in G:
        for me in E:
            s = M1_transition_strength_avg(I_0, L_0, S_0, J_0, F_0, mg, L_1, S_1, J_1, F_1, me)
            if s != 0:
                flg = True
                print "{0:.6f} {1:} {2:}".format(s, mg, me)
    if not flg:
        print "no allowed transitions"
    print "E1"
    flg = False
    for mg in G:
        for me in E:
            s = E1_transition_strength_avg(I_0, L_0, S_0, J_0, F_0, mg, L_1, S_1, J_1, F_1, me)
            if s != 0:
                flg = True
                print "{0:.6f} {1:} {2:}".format(s, mg, me)
    if not flg:
        print "no allowed transitions"
