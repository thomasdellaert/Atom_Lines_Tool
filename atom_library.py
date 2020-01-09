from atoms import *
from parsers import parse_NIST_levels

df_Yb = df = parse_NIST_levels("YbII_NIST_levels.csv")

# TODO: allow this to take its A and B coeffs from a csv.
def populate_levels(df, atom, I=0, default_A=0, default_B = 0, nlevels="all"):
    ls = []
    if nlevels == "all":
        nlevels = len(df["index"])
    for i in range(nlevels):
        try:
            level = EnergyLevel(df, i, I=I, A_coeff=default_A, B_coeff=default_B)
            ls.append(level)
        except:
            pass
    atom.add_level(ls)


Yb_173 = Atom(name="173Yb")
S12_173 = EnergyLevel(df, 0, I=2.5, A_coeff=-3.4975e-3, B_coeff=0.0)  # 2S1/2
P12_173 = EnergyLevel(df, 8, I=2.5, A_coeff=-0.5812e-3, B_coeff=0.0)  # 2P1/2
P32_173 = EnergyLevel(df, 9, I=2.5, A_coeff=-0.245e-3, B_coeff=1.46e-3)  # 2P3/2
F72_173 = EnergyLevel(df, 1, I=2.5, A_coeff=-0.24e-3, B_coeff=-4.762e-3)  # 2F7/2
D32_173 = EnergyLevel(df, 3, I=2.5, A_coeff=-0.11031e-3, B_coeff=0.09514e-3)  # 2D3/2
D52_173 = EnergyLevel(df, 4, I=2.5, A_coeff=0.00347e-3, B_coeff=1.1904e-3)  # 2D5/2
B32_173 = EnergyLevel(df, 6, I=2.5, A_coeff=-0.01e-3, B_coeff=0.0)  # 1D[3/2]5/2 random guess
B52_173 = EnergyLevel(df, 5, I=2.5, A_coeff=-0.055e-3, B_coeff=-1.72e-3)  # 1D[5/2]5/2

levels_173 = (S12_173, P12_173, P32_173, F72_173, D32_173, D52_173, B32_173, B52_173)
Yb_173.add_level(levels_173)

repump_173 =           Transition(level_0=D52_173, F_0=4, m_F_0=0, level_1=P32_173, F_1=3, m_F_1=0)
cycling_173 =          Transition(level_0=S12_173, F_0=2, m_F_0=0, level_1=P12_173, F_1=3, m_F_1=0)
four_sixty_seven_173 = Transition(level_0=S12_173, F_0=2, m_F_0=0, level_1=F72_173, F_1=5, m_F_1=0)
four_eleven_173 =      Transition(level_0=S12_173, F_0=2, m_F_0=0, level_1=D52_173, F_1=3, m_F_1=0)

transitions_173 = (repump_173, cycling_173, four_sixty_seven_173, four_eleven_173)
Yb_173.add_transition(transitions_173)

Yb_171 = Atom(name="171Yb")
S12_171 = EnergyLevel(df, 0, I=0.5, A_coeff=12.645e-3)  # 2S1/2
P12_171 = EnergyLevel(df, 8, I=0.5, A_coeff=2.1079e-3)  # 2P1/2
P32_171 = EnergyLevel(df, 9, I=0.5, A_coeff=0.8754e-3)  # 2P3/2
F72_171 = EnergyLevel(df, 1, I=0.5, A_coeff=1.105e-3)  # 2F7/2
D32_171 = EnergyLevel(df, 3, I=0.5, A_coeff=0.40048e-3)  # 2D3/2
D52_171 = EnergyLevel(df, 4, I=0.5, A_coeff=-0.01258e-3)  # 2D5/2
B32_171 = EnergyLevel(df, 6, I=0.5, A_coeff=0.01e-3)  # 1D[3/2]5/2 random guess
B52_171 = EnergyLevel(df, 5, I=0.5, A_coeff=-0.199e-3)  # 1D[5/2]5/2

levels_171 = (S12_171, P12_171, P32_171, F72_171, D32_171, D52_171, B32_171, B52_171)
Yb_171.add_level(levels_171)

repump_171 =      Transition(level_0=D52_171, F_0=2, m_F_0=0, level_1=P32_171, F_1=1, m_F_1=0)
cycling_171 =     Transition(level_0=S12_171, F_0=0, m_F_0=0, level_1=P12_171, F_1=1, m_F_1=0)
four_eleven_171 = Transition(level_0=S12_171, F_0=0, m_F_0=0, level_1=D52_171, F_1=3, m_F_1=0)
test =            Transition(level_0=S12_171, F_0=0, m_F_0=0, level_1=D32_171, F_1=2, m_F_1=0)

transitions_171 = (repump_171, cycling_171, four_eleven_171, test)
Yb_171.add_transition(transitions_171)


Yb_174 = Atom(name="174Yb")
S12_174 = EnergyLevel(df, 0)
P12_174 = EnergyLevel(df, 1)
P32_174 = EnergyLevel(df, 8)
F72_174 = EnergyLevel(df, 9)
D32_174 = EnergyLevel(df, 3)
D52_174 = EnergyLevel(df, 4)
B32_174 = EnergyLevel(df, 6)
B52_174 = EnergyLevel(df, 5)

levels_174 = (S12_174, P12_174, P32_174, F72_174, D32_174, D52_174, B32_174, B52_174)
Yb_174.add_level(levels_174)

# Yb_174.add_transition([Transition(level_0=S12_174, F_0=0.5, m_F_0=0.5, level_1=P12_174, F_1=0.5, m_F_1=0.5)])

Yb_174_big = Atom(name="174Yb")
populate_levels(df, Yb_174_big, I=0, nlevels=30)

# TODO: these As currently aren't being overwritten.
Yb_171_big = Atom(name="171Yb")
populate_levels(df, Yb_171_big, I=0.5, nlevels=30, default_A=1e-3)
Yb_171_big.levels["2S1/2"].set_coeffs(A_coeff=12.642812e-3)
Yb_171_big.levels["2P*1/2"].set_coeffs(A_coeff=2.1079e-3)
Yb_171_big.levels["2P*3/2"].set_coeffs(A_coeff=0.8754e-3)
Yb_171_big.levels["2F*7/2"].set_coeffs(A_coeff=1.105e-3)
Yb_171_big.levels["2D3/2"].set_coeffs(A_coeff=0.40048e-3)
Yb_171_big.levels["2D5/2"].set_coeffs(A_coeff=-0.01258e-3)
Yb_171_big.levels["3[3/2]*5/2"].set_coeffs(A_coeff=0.199e-3)

t = []
t.append(Transition(level_0=Yb_171_big.levels["2D5/2"], F_0=2, m_F_0=0, level_1=Yb_171_big.levels["2P*3/2"], F_1=1, m_F_1=0))
t.append(Transition(level_0=Yb_171_big.levels["2S1/2"], F_0=0, m_F_0=0, level_1=Yb_171_big.levels["2P*1/2"], F_1=1, m_F_1=0))
t.append(Transition(level_0=Yb_171_big.levels["2S1/2"], F_0=0, m_F_0=0, level_1=Yb_171_big.levels["2D5/2"], F_1=3, m_F_1=0))
t.append(Transition(level_0=Yb_171_big.levels["2S1/2"], F_0=0, m_F_0=0, level_1=Yb_171_big.levels["2D3/2"], F_1=2, m_F_1=0))
Yb_171_big.add_transition(t)

Yb_173_big = Atom(name="173Yb")
populate_levels(df, Yb_173_big, I=2.5, nlevels=30)


if __name__ == "__main__":
    print Yb_171_big.levels.keys()
