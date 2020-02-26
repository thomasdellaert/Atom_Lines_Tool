from pandas import DataFrame
import pandas as pd
from parsers import term_frac


class EnergyLevel:
    # TODO: isotope shifts?
    def __init__(self, df, df_index, name="term", I=0.0, A_coeff=0.0, B_coeff=0.0):
        self._initialized = False
        # get parameters from datafile
        my_row = df.iloc[df_index]

        (self.configuration, self.term, self.S, self.L, self.K,
         self.j1, self.j2, self.J, self.parity, self.level, self.lande) = my_row.values
        if name == "term":
            self.name = self.term + term_frac(self.J)
        elif name == "full":
            self.name = self.configuration + " " + self.term + term_frac(self.J)
        else:
            self.name = name
        self.A_coeff = A_coeff
        self.B_coeff = B_coeff
        self.coupling = self.get_coupling()
        self.I = I
        self.Fs = self.get_Fs()
        self.hf_levels, self.hf_shifts = self.get_hyperfine_data()
        self.z_levels, self.z_shifts = self.get_zeeman_data()

        self._initialized = True

    def __setattr__(self, key, value):
        self.__dict__[key] = value
        if self._initialized:
            if key in ['level', 'lande', 'A_coeff', 'B_coeff']:
                self.hf_levels, self.hf_shifts = self.get_hyperfine_data()
            if key in ['level', 'lande']:
                self.z_levels, self.z_shifts = self.get_zeeman_data()

    def get_Fs(self):
        Fs = []
        F = abs(self.I - self.J)
        while F <= abs(self.I + self.J):
            Fs.append(F)
            F += 1
        return Fs

    def get_coupling(self):
        if self.L is not None:
            coupling = "LS"
        elif self.K is not None:
            coupling = "JK"
        else:
            coupling = "jj"
        return coupling

    def get_hyperfine_data(self):
        hf_levels = {}
        hf_shifts = {}
        if self.I == 0:
            hf_levels[self.J] = self.level
            hf_shifts[self.J] = 0
        elif self.I <= 1:
            for F in self.Fs:
                K1 = (F * (F + 1) - self.I * (self.I + 1) - self.J * (self.J + 1))
                shift = (self.A_coeff * K1 / 2)
                hf_shifts[F] = shift
                hf_levels[F] = shift + self.level
        else:
            for F in self.Fs:
                K1 = (F * (F + 1) - self.I * (self.I + 1) - self.J * (self.J + 1))
                K2 = (3 * K1 * (K1 + 1) / 2 - 2 * self.J * (self.J + 1) * self.I * (self.I + 1)) / (
                        self.J * (2 * self.J - 1) * self.I * (2 * self.I - 1))
                shift = (self.A_coeff * K1 / 2) + (self.B_coeff * K2 / 4)
                hf_shifts[F] = shift
                hf_levels[F] = shift + self.level
        return hf_levels, hf_shifts

    def get_zeeman_data(self):
        z_levels = {}
        z_shifts = {}
        for F in self.Fs:
            try:
                g_F = self.lande * (F * (F + 1) + self.J * (self.J + 1) - self.I * (self.I + 1)) / (2 * F * (F + 1))
            except ZeroDivisionError:
                g_F = 0
            z_levels[F] = {}
            z_shifts[F] = {}
            m_F = -F
            while m_F <= F:
                z_shifts[F][m_F] = g_F * m_F * 1.4e-6
                z_levels[F][m_F] = self.hf_levels[F] + z_shifts[F][m_F]
                m_F += 1
        return z_levels, z_shifts

    def set_name(self, name):
        """Sets the name of the level. Can be 'term', 'full', or a custom name"""
        if name == "term":
            self.name = self.term + term_frac(self.J)
        elif name == "full":
            self.name = self.configuration + " " + self.term + term_frac(self.J)
        else:
            self.name = name

    def data_table(self, hf=True, zeeman=True):
        table = DataFrame(columns=["configuration", "term", "level",
                                   "J", "F", "m_F", "J_frac", "F_frac", "m_F_frac", "hf", "z"])
        if not hf:
            line = DataFrame(data={"configuration": [self.configuration], "term": [self.term], "level": self.level,
                                   "J": [self.J], "F": [None], "m_F": [None],
                                   "J_frac": [term_frac(self.J)], "F_frac": [None], "m_F_frac": [None],
                                   "hf": [0.0], "z": [0.0]})
            table = table.append(line, ignore_index=True)
        else:
            for F in self.Fs:
                if not zeeman:
                    hyperfine = self.hf_shifts[F]
                    level = self.hf_levels[F]
                    line = DataFrame(data={"configuration": [self.configuration], "term": [self.term], "level": level,
                                           "J": [self.J], "F": [F], "m_F": [None],
                                           "J_frac": [term_frac(self.J)], "F_frac": [term_frac(F)], "m_F_frac": [None],
                                           "hf": [hyperfine], "z": [0.0]})
                    table = table.append(line, ignore_index=True)
                else:
                    for m_F in self.z_shifts[F].keys():
                        z = self.z_shifts[F][m_F]
                        hyperfine = self.hf_shifts[F]
                        level = self.level + self.hf_shifts[F] + self.z_shifts[F][m_F]
                        line = DataFrame(
                            data={"configuration": [self.configuration], "term": [self.term], "level": level,
                                  "J": [self.J], "F": [F], "m_F": [m_F],
                                  "J_frac": [term_frac(self.J)], "F_frac": [term_frac(F)], "m_F_frac": [term_frac(m_F)],
                                  "hf": [hyperfine], "z": [z]})
                        table = table.append(line, ignore_index=True)
        return table


class Transition:
    def __init__(self, level_0, F_0, m_F_0, level_1, F_1, m_F_1):
        self.level_0, self.level_1 = level_0, level_1
        self.F_0, self.m_F_0, self.F_1, self.m_F_1 = F_0, m_F_0, F_1, m_F_1
        self.J_0, self.J_1 = level_0.J, level_1.J
        self.L_0, self.L_1 = level_0.L, level_1.L
        self.parity_0, self.parity_1 = level_0.parity, level_1.parity
        self.name = level_0.name + str(F_0) + str(m_F_0) + "->" + level_1.name + str(F_1) + str(m_F_1)

        self.transition_table = self.data_table()

    def data_table(self):

        table_0 = self.level_0.data_table()
        table_1 = self.level_1.data_table()
        line_0 = table_0.loc[(table_0["m_F"] == self.m_F_0) & (table_0["F"] == self.F_0)]
        line_1 = table_1.loc[(table_1["m_F"] == self.m_F_1) & (table_1["F"] == self.F_1)]
        line_0 = line_0.drop(["J_frac", "F_frac", "m_F_frac"], axis=1)
        line_1 = line_1.drop(["J_frac", "F_frac", "m_F_frac"], axis=1)
        if line_0.empty or line_1.empty:
            raise ValueError("The selected quantum numbers don't yield a transition. Check that they exist in the specified levels")
        line_0.columns = [str(col) + '_0' for col in line_0.columns]
        line_1.columns = [str(col) + '_1' for col in line_1.columns]
        line_0 = line_0.reset_index(drop=True)
        line_1 = line_1.reset_index(drop=True)

        data_table = pd.concat([line_0, line_1], axis=1, ignore_index=False)

        delta_l = abs(data_table["level_0"] - data_table["level_1"])
        wavelength = 299792.458/delta_l

        data_table["delta_l"] = [delta_l]
        data_table["wavelength"] = [wavelength]
        data_table["name"] = [self.name]
        return data_table

    def get_type(self):
        # TODO: finish this method sometime when I have time to figure out what all the relevant selection rules are
        d_m = self.m_F_1 - self.m_F_0
        d_F = self.F_1 - self.F_0
        d_J = self.J_1 - self.J_0

        transition_type = "unknown"
        if abs(d_J) <= 1 and not (self.J_0 == 0  and self.J_1 == 0):
            if abs(d_F) <= 1:
                if self.parity_0 != self.parity_1:
                    transition_type = "E1"

        return "dm: {}, dF: {}, dJ:{}".format(d_m, d_F, d_J), transition_type


class Atom:
    def __init__(self, name, levels=(), transitions=()):
        self.name = name
        self.levels = {}
        for level in levels:
            if level.name in self.levels.keys():
                raise Exception("{} is already a level in this atom. ".format(level.name))
            self.levels[level.name] = level

        self.transitions = {}
        for transition in transitions:
            self.transitions[transition.name] = transition
        if len(levels) > 0:
            self.rezero()

    def rezero(self):
        """Takes the minimum hyperfine level in the atom and sets it to zero, shifting all others appropriately"""
        gs = self.levels.values()[0]
        gs_level = min(gs.hf_levels.values())
        for state in self.levels.values():
            if state.level < gs_level:
                gs = state
                gs_level = min(gs.hf_levels.values())
        for state in self.levels.values():
            state.level -= gs_level

    def add_level(self, levels):
        for level in levels:
            self.levels[level.name] = level
        self.rezero()

    def remove_level(self, mode="keys", *levels):
        if mode == "keys":
            for key in levels:
                self.levels.pop(key)
        elif mode == "levels":
            for level in levels:
                self.levels.pop(level.name)
        else:
            raise ValueError("Unrecognized mode")
        self.rezero()

    def add_transition(self, transitions):
        for transition in transitions:
            self.transitions[transition.name] = transition

    def remove_transition(self, mode="keys", *transitions):
        if mode == "keys":
            for key in transitions:
                self.transitions.pop(key)
        elif mode == "transitions":
            for transition in transitions:
                self.transitions.pop(transition.name)
        else:
            raise ValueError("Unrecognized mode")


if __name__ == "__main__":
    from parsers import parse_NIST_levels

    df = parse_NIST_levels("YbII_NIST_levels.csv")

    S12_171 = EnergyLevel(df, 0, I=0.5, A_coeff=12.645e-3)  # 2S1/2
    P12_171 = EnergyLevel(df, 8, I=0.5, A_coeff=2.1079e-3)  # 2P1/2
    cycling_171 = Transition(level_0=S12_171, F_0=0, m_F_0=0, level_1=P12_171, F_1=1, m_F_1=0)
    print cycling_171.name
