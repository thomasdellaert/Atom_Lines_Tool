from atoms import *
from parsers import parse_NIST_levels

df_Yb = df = parse_NIST_levels("YbII_NIST_levels.csv")


# TODO: allow this to take its A and B coeffs from a csv.
def populate_levels(df, atom, I=0, default_A=0, default_B=0, nlevels="all"):
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


Yb_174 = Atom(name="174Yb")
populate_levels(df, Yb_174, I=0, nlevels=30)

Yb_173 = Atom(name="173Yb")
populate_levels(df, Yb_173, I=2.5, nlevels=30, default_A=0.01e-3, default_B=0.0)
Yb_173.levels["2S1/2"].set_coeffs(A_coeff=-3.4975e-3, B_coeff=0.0)
Yb_173.levels["2P*1/2"].set_coeffs(A_coeff=-0.5812e-3, B_coeff=0.0)
Yb_173.levels["2P*3/2"].set_coeffs(A_coeff=-0.245e-3, B_coeff=1.46e-3)
Yb_173.levels["2F*7/2"].set_coeffs(A_coeff=-0.24e-3, B_coeff=4.762e-3)
Yb_173.levels["2D3/2"].set_coeffs(A_coeff=-0.11031e-3, B_coeff=0.09514e-3)
Yb_173.levels["2D5/2"].set_coeffs(A_coeff=0.00347e-3, B_coeff=1.1904e-3)
Yb_173.levels["3[3/2]*5/2"].set_coeffs(A_coeff=-0.055e-3, B_coeff=-1.72e-3)

t = [
    Transition(level_0=Yb_173.levels["2D5/2"], F_0=4, m_F_0=0, level_1=Yb_173.levels["2P*3/2"], F_1=3, m_F_1=0),
    Transition(level_0=Yb_173.levels["2S1/2"], F_0=3, m_F_0=0, level_1=Yb_173.levels["2P*1/2"], F_1=3, m_F_1=0),
    Transition(level_0=Yb_173.levels["2S1/2"], F_0=3, m_F_0=0, level_1=Yb_173.levels["2D5/2"], F_1=3, m_F_1=0),
    Transition(level_0=Yb_173.levels["2S1/2"], F_0=3, m_F_0=0, level_1=Yb_173.levels["2F*7/2"], F_1=5, m_F_1=0)]
Yb_173.add_transition(t)

Yb_171 = Atom(name="171Yb")
populate_levels(df, Yb_171, I=0.5, nlevels=30, default_A=0.01e-3)
Yb_171.levels["2S1/2"].set_coeffs(A_coeff=12.642812e-3)
Yb_171.levels["2P*1/2"].set_coeffs(A_coeff=2.1079e-3)
Yb_171.levels["2P*3/2"].set_coeffs(A_coeff=0.8754e-3)
Yb_171.levels["2F*7/2"].set_coeffs(A_coeff=1.105e-3)
Yb_171.levels["2D3/2"].set_coeffs(A_coeff=0.40048e-3)
Yb_171.levels["2D5/2"].set_coeffs(A_coeff=-0.01258e-3)
Yb_171.levels["3[3/2]*5/2"].set_coeffs(A_coeff=0.199e-3)

t = [
    Transition(level_0=Yb_171.levels["2D5/2"], F_0=2, m_F_0=0, level_1=Yb_171.levels["2P*3/2"], F_1=1, m_F_1=0),
    Transition(level_0=Yb_171.levels["2S1/2"], F_0=0, m_F_0=0, level_1=Yb_171.levels["2P*1/2"], F_1=1, m_F_1=0),
    Transition(level_0=Yb_171.levels["2S1/2"], F_0=0, m_F_0=0, level_1=Yb_171.levels["2D5/2"], F_1=3, m_F_1=0),
    Transition(level_0=Yb_171.levels["2S1/2"], F_0=0, m_F_0=0, level_1=Yb_171.levels["2D3/2"], F_1=2, m_F_1=0)]
Yb_171.add_transition(t)

Yb_173_big = Atom(name="173Yb")
populate_levels(df, Yb_173_big, I=2.5, nlevels=30)

if __name__ == "__main__":
    print Yb_173.levels.keys()
