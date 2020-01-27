from atoms import *
from parsers import parse_NIST_levels

df_Yb = parse_NIST_levels("YbII_NIST_levels.csv")


# TODO: allow this to take its A and B coeffs from a csv.
def _populate_levels(df, atom, I=0.0, default_A=0.0, default_B=0.0, n_levels=-1):
    ls = []
    if n_levels == -1:
        n_levels = len(df["index"])
    for i in range(n_levels):
        try:
            level = EnergyLevel(df, i, I=I, A_coeff=default_A, B_coeff=default_B)
            ls.append(level)
        except Exception as e:
            print e.__doc__
            print e.message
    atom.add_level(ls)


# region define 174Yb
Yb_174 = Atom(name="174Yb")
_populate_levels(df_Yb, Yb_174, I=0, n_levels=30)
Yb_174.rezero()
# endregion

# region define 173Yb
Yb_173 = Atom(name="173Yb")
_populate_levels(df_Yb, Yb_173, I=2.5, n_levels=30, default_A=0.01e-3, default_B=0.0)
Yb_173.levels["2S1/2"].set_coeffs(A_coeff=-3.4975e-3, B_coeff=0.0)
Yb_173.levels["2P*1/2"].set_coeffs(A_coeff=-0.5812e-3, B_coeff=0.0)
Yb_173.levels["2P*3/2"].set_coeffs(A_coeff=-0.245e-3, B_coeff=1.46e-3)
Yb_173.levels["2F*7/2"].set_coeffs(A_coeff=-0.24e-3, B_coeff=4.762e-3)
Yb_173.levels["2D3/2"].set_coeffs(A_coeff=-0.11031e-3, B_coeff=0.09514e-3)
Yb_173.levels["2D5/2"].set_coeffs(A_coeff=0.00347e-3, B_coeff=1.1904e-3)
Yb_173.levels["3[3/2]*5/2"].set_coeffs(A_coeff=-0.055e-3, B_coeff=-1.72e-3)
Yb_173.rezero()

_t = [
    Transition(level_0=Yb_173.levels["2D5/2"], F_0=4, m_F_0=0, level_1=Yb_173.levels["2P*3/2"], F_1=3, m_F_1=0),
    Transition(level_0=Yb_173.levels["2S1/2"], F_0=3, m_F_0=0, level_1=Yb_173.levels["2P*1/2"], F_1=3, m_F_1=0),
    Transition(level_0=Yb_173.levels["2S1/2"], F_0=3, m_F_0=0, level_1=Yb_173.levels["2D5/2"], F_1=3, m_F_1=0),
    Transition(level_0=Yb_173.levels["2S1/2"], F_0=3, m_F_0=0, level_1=Yb_173.levels["2F*7/2"], F_1=5, m_F_1=0)]
Yb_173.add_transition(_t)
# endregion

# region define 171Yb
Yb_171 = Atom(name="171Yb")
_populate_levels(df_Yb, Yb_171, I=0.5, n_levels=30, default_A=0.01e-3)
Yb_171.levels["2S1/2"].set_coeffs(A_coeff=12.642812e-3)
Yb_171.levels["2P*1/2"].set_coeffs(A_coeff=2.1079e-3)
Yb_171.levels["2P*3/2"].set_coeffs(A_coeff=0.8754e-3)
Yb_171.levels["2F*7/2"].set_coeffs(A_coeff=1.105e-3)
Yb_171.levels["2D3/2"].set_coeffs(A_coeff=0.40048e-3)
Yb_171.levels["2D5/2"].set_coeffs(A_coeff=-0.01258e-3)
Yb_171.levels["3[3/2]*5/2"].set_coeffs(A_coeff=0.199e-3)
Yb_171.levels["1[3/2]*3/2"].set_coeffs(A_coeff=4.45e-3)
Yb_171.rezero()

_t = [
    Transition(level_0=Yb_171.levels["2S1/2"], F_0=1, m_F_0=0, level_1=Yb_171.levels["2P*1/2"], F_1=1, m_F_1=0),
    Transition(level_0=Yb_171.levels["2S1/2"], F_0=1, m_F_0=0, level_1=Yb_171.levels["2D5/2"], F_1=3, m_F_1=0),
    Transition(level_0=Yb_171.levels["3[3/2]*1/2"], F_0=0, m_F_0=0, level_1=Yb_171.levels["2D3/2"], F_1=2, m_F_1=0),
    Transition(level_0=Yb_171.levels["2F*7/2"], F_0=4, m_F_0=0, level_1=Yb_171.levels["1[3/2]*3/2"], F_1=2, m_F_1=0),
    Transition(level_0=Yb_171.levels["2F*7/2"], F_0=3, m_F_0=0, level_1=Yb_171.levels["1[3/2]*3/2"], F_1=1, m_F_1=0),
    Transition(level_0=Yb_171.levels["2D5/2"], F_0=3, m_F_0=0, level_1=Yb_171.levels["2F*7/2"], F_1=3, m_F_1=0),
    Transition(level_0=Yb_171.levels["1[3/2]*3/2"], F_0=2, m_F_0=0, level_1=Yb_171.levels["2S1/2"], F_1=1, m_F_1=0),
    Transition(level_0=Yb_171.levels["1[3/2]*3/2"], F_0=2, m_F_0=0, level_1=Yb_171.levels["2D5/2"], F_1=3, m_F_1=0)]
Yb_171.add_transition(_t)
# endregion


if __name__ == "__main__":
    print Yb_171.levels.keys()
