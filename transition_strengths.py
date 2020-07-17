from sympy.physics.wigner import wigner_3j, wigner_6j

def rel_transiton_strength(I, q, J1, F0, M0, F1, M1):
    return float(3*(2*J1+1)*(2*F0+1)*(2*F1+1)*\
           wigner_6j(1, 0, 1, 0.5, J1, 0.5, prec=8)**2*\
           wigner_6j(J1, 0.5, 1, F0, F1, I, prec=8)**2*\
           wigner_3j(F1, 1, F0, -M1, q, M0)**2)

if __name__ == '__main__':
    for (Mg, Me) in [(-3, -4), (-2, -3), (-3, -3), (-1, -2), (-2, -2),  (0, -1), (-1, -1), (1, 0), (-2, -1), (2, 1), (-1, 0), (3, 2), (0, 1), (1, 2), (2, 3), (3, 4)]:
        s = float(rel_transiton_strength(2.5, Me-Mg, 1.5, 3, Mg, 4, Me))
        print "{0:.6f} {1:} {2:}".format(s, Mg, Me)