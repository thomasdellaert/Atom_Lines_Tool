from atoms import *
from parsers import parse_NIST_levels
from pandas import read_csv


def populate_levels(df, atom, I=0.0, default_A=0.0, default_B=0.0, n_levels=0, hf_source=None):
    ls = []
    hfs = None
    print(n_levels)
    if n_levels < 1:
        print(len(df))
        n_levels = len(df)
    print(n_levels)
    if I != 0 and hf_source is not (None or ""):
        hfs = read_csv(hf_source, index_col=0)
    for i in range(n_levels):
        try:
            level = EnergyLevel(df, i, I=I)
            if hfs is not None:
                try:
                    level.A_coeff = hfs.at[level.name, 'A_coeff']
                    level.B_coeff = hfs.at[level.name, 'B_coeff']
                except KeyError:
                    level.A_coeff = default_A
                    level.B_coeff = default_B
            ls.append(level)
        except Exception as e:
            print("Error parsing line {}".format(i))
            print(e.__doc__)
            print(str(e))
    atom.add_level(ls)
    atom.rezero()

# TODO: populate_transitions()


if __name__ == "__main__":

    df_Yb = parse_NIST_levels('YbII_NIST_levels.csv')

    # region define 174Yb
    print('Initializing 174Yb')
    Yb_174 = Atom(name='174Yb')
    populate_levels(df_Yb, Yb_174, I=0, n_levels=30)
    Yb_174.rezero()
    print('174Yb Initialized')
    # endregion

    # region define 173Yb
    print('Initializing 173Yb')
    Yb_173 = Atom(name='173Yb')
    populate_levels(df_Yb, Yb_173, I=2.5, n_levels=10, default_A=0.01e-3, default_B=0.0, hf_source='173Yb_Hyperfine.csv')
    Yb_173.add_transition([
        Transition(level_0=Yb_173.levels['2D5/2'], F_0=4, m_F_0=0, level_1=Yb_173.levels['2P*3/2'], F_1=3, m_F_1=0),
        Transition(level_0=Yb_173.levels['2S1/2'], F_0=3, m_F_0=0, level_1=Yb_173.levels['2P*1/2'], F_1=3, m_F_1=0),
        Transition(level_0=Yb_173.levels['2S1/2'], F_0=3, m_F_0=0, level_1=Yb_173.levels['2D5/2'], F_1=3, m_F_1=0),
        Transition(level_0=Yb_173.levels['2S1/2'], F_0=3, m_F_0=0, level_1=Yb_173.levels['2F*7/2'], F_1=5, m_F_1=0)]
    )
    print('173Yb Initialized')

    # endregion

    # region define 171Yb
    print('Initializing 171Yb')
    Yb_171 = Atom(name='171Yb')
    populate_levels(df_Yb, Yb_171, I=0.5, n_levels=30, default_A=0.01e-3, hf_source='171Yb_Hyperfine.csv')

    Yb_171.add_transition([
        Transition(level_0=Yb_171.levels['2S1/2'], F_0=1, m_F_0=0, level_1=Yb_171.levels['2P*1/2'], F_1=1, m_F_1=0),
        Transition(level_0=Yb_171.levels['2S1/2'], F_0=1, m_F_0=0, level_1=Yb_171.levels['2D5/2'], F_1=3, m_F_1=0),
        Transition(level_0=Yb_171.levels['3[3/2]*1/2'], F_0=0, m_F_0=0, level_1=Yb_171.levels['2D3/2'], F_1=2, m_F_1=0),
        Transition(level_0=Yb_171.levels['2F*7/2'], F_0=4, m_F_0=0, level_1=Yb_171.levels['1[3/2]*3/2'], F_1=2, m_F_1=0),
        Transition(level_0=Yb_171.levels['2F*7/2'], F_0=3, m_F_0=0, level_1=Yb_171.levels['1[3/2]*3/2'], F_1=1, m_F_1=0),
        Transition(level_0=Yb_171.levels['2D5/2'], F_0=3, m_F_0=0, level_1=Yb_171.levels['2F*7/2'], F_1=3, m_F_1=0),
        Transition(level_0=Yb_171.levels['1[3/2]*3/2'], F_0=2, m_F_0=0, level_1=Yb_171.levels['2S1/2'], F_1=1, m_F_1=0),
        Transition(level_0=Yb_171.levels['1[3/2]*3/2'], F_0=2, m_F_0=0, level_1=Yb_171.levels['2D5/2'], F_1=3, m_F_1=0)]
    )
    print('171Yb Initialized')
    # endregion


