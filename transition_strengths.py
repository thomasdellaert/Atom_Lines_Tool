from sympy.physics.wigner import wigner_3j, wigner_6j

def rel_transiton_strength(I, q, J1, F0, M0, F1, M1):
    return float(3*(2*J1+1)*(2*F0+1)*(2*F1+1)*\
           wigner_6j(1, 0, 1, 0.5, J1, 0.5, prec=8)**2*\
           wigner_6j(J1, 0.5, 1, F0, F1, I, prec=8)**2*\
           wigner_3j(F1, 1, F0, -M1, q, M0)**2)

if __name__ == '__main__':
    for Me in {-4, -3, -2, -1, 0, 1, 2, 3, 4}:
        for Mg in {-3, -2, -1, 0, 1, 2, 3}:
            for q  in {-1, 0, 1}:
                s = float(rel_transiton_strength(2.5, q, 1.5, 3, Mg, 4, Me))
                if s != 0:
                    print "{0:.6f} {1:} {2:}".format(s, Mg, Me)