from pandas import DataFrame
from scipy import interpolate
import pandas as pd


def term_frac(term):
    num = int(term / 0.5)
    if num % 2 != 0:
        return str(num) + "/2"
    else:
        return str(num / 2)


class EnergyLevel:
    # TODO: isotope shifts! Frequencies are currently off by like terahertz
    def __init__(self, df, df_index, name="term", A_coeff=0.0, B_coeff=0.0, I=0.0):
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

        self.coupling = self.get_coupling()
        self.A_coeff = A_coeff
        self.B_coeff = B_coeff
        self.I = I
        self.Fs = self.get_Fs()
        self.hf_levels, self.hf_shifts = self.get_hyperfine_data()
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
        else:
            for F in self.Fs:
                K1 = (F * (F + 1) - self.I * (self.I + 1) - self.J * (self.J + 1))
                K2 = self.I * (2 * self.I - 1) * self.J * (2 * self.J - 1) / (
                            K1 * (K1 + 1) - 4 * self.I * (self.I + 1) * self.J * (self.J + 1))
                shift = (self.A_coeff * K1 / 2) + (self.B_coeff * 8 * K2 / 3)
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

    def rename(self, n):
        self.name = n

    # def physics_table(self):
    #     table = DataFrame(columns=["configuration", "term", "level",
    #                                "J", "F", "m_F", "y0", "hf", "z"])
    #     for F in self.Fs:
    #         for m_F in self.z_shifts[F].keys():
    #             z = self.z_shifts[F][m_F]
    #             hf = self.hf_shifts[F]
    #             y = y0 + self.hf_shifts[F] * scale + self.z_shifts[F][m_F] * scale
    #             level = self.level + self.hf_shifts[F] + self.z_shifts[F][m_F]
    #             line = DataFrame(
    #                 data={"configuration": [self.configuration], "term": [self.term], "level": level,
    #                       "J"            : [self.J], "F": [F], "m_F": [m_F],
    #                       "J_frac"       : [term_frac(self.J)], "F_frac": [term_frac(F)],
    #                       "m_F_frac"     : [term_frac(m_F)],
    #                       "color"        : [params["color"]], "y0": [y0], "hf": [hf], "z": [z],
    #                       "y"            : [y], "x0": [x], "x1": [x1]})
    #             table = table.append(line, ignore_index=True)

    def plottable_table(self, **kwargs):
        params = {"width"            : 1.0,
                  "sublevel_spacing" : 0.03,
                  "zeeman"           : True,
                  "compact"          : True,
                  "hyperfine"        : True,
                  "scale_splitting"  : 1.0,
                  "override_position": False,
                  "offset_position"  : (0.0, 0.0),
                  "color"            : "black"}
        for param in kwargs.keys():
            params[param] = kwargs[param]
        if not params["override_position"]:
            x0 = self.J + params["offset_position"][0] - params["width"] / 2
            y0 = self.level + params["offset_position"][1]
        else:
            x0, y0 = params["override_position"]

        table = DataFrame(columns=["configuration", "term", "level",
                                   "J", "F", "m_F", "J_frac", "F_frac", "m_F_frac",
                                   "color", "y0", "hf", "z",
                                   "y", "x0", "x1"])
        scale = params["scale_splitting"]

        if not params["hyperfine"]:
            x1 = x0 + params["width"]
            line = DataFrame(data={"configuration": [self.configuration], "term": [self.term], "level": self.level,
                                   "J"            : [self.J], "F": [None], "m_F": [None],
                                   "J_frac"       : [term_frac(self.J)], "F_frac": [None], "m_F_frac": [None],
                                   "color"        : [params["color"]], "y0": [y0], "hf": [0.0], "z": [0.0],
                                   "y"            : [y0], "x0": [x0], "x1": [x1]})
            table = table.append(line)
        else:
            for F in self.Fs:
                if not params["zeeman"]:
                    y = y0 + self.hf_shifts[F] * scale
                    hf = self.hf_shifts[F]
                    level = self.level + self.hf_shifts[F]
                    x1 = x0 + params["width"]
                    line = DataFrame(data={"configuration": [self.configuration], "term": [self.term], "level": level,
                                           "J"            : [self.J], "F": [F], "m_F": [None],
                                           "J_frac"       : [term_frac(self.J)], "F_frac": [term_frac(F)],
                                           "m_F_frac"     : [None],
                                           "color"        : [params["color"]], "y0": [y0], "hf": [hf], "z": [0.0],
                                           "y"            : [y], "x0": [x0], "x1": [x1]})
                    table = table.append(line, ignore_index=True)
                else:
                    delta = params["sublevel_spacing"]
                    Fmax = max(self.Fs)
                    wd = (params["width"] - delta * 2 * Fmax) / (2 * Fmax + 1)
                    for m_F in self.z_shifts[F].keys():
                        z = self.z_shifts[F][m_F]
                        hf = self.hf_shifts[F]
                        y = y0 + self.hf_shifts[F] * scale + self.z_shifts[F][m_F] * scale
                        level = self.level + self.hf_shifts[F] + self.z_shifts[F][m_F]
                        x = x0 + (m_F + Fmax) * (wd + delta)
                        x1 = x + wd
                        line = DataFrame(
                            data={"configuration": [self.configuration], "term": [self.term], "level": level,
                                  "J"            : [self.J], "F": [F], "m_F": [m_F],
                                  "J_frac"       : [term_frac(self.J)], "F_frac": [term_frac(F)],
                                  "m_F_frac"     : [term_frac(m_F)],
                                  "color"        : [params["color"]], "y0": [y0], "hf": [hf], "z": [z],
                                  "y"            : [y], "x0": [x], "x1": [x1]})
                        table = table.append(line, ignore_index=True)
        return table


class Transition:
    def __init__(self, level_0, F_0, m_F_0, level_1, F_1, m_F_1):
        self.level_0, self.level_1 = level_0, level_1
        self.F_0, self.m_F_0, self.F_1, self.m_F_1 = F_0, m_F_0, F_1, m_F_1
        self.J_0, self.J_1 = level_0.J, level_1.J
        self.L_0, self.L_1 = level_0.L, level_1.L
        self.parity_0, self.parity_1 = level_0.parity, level_1.parity

        self.transition_table = self.plottable_table()

        self.name = str(self.transition_table["term_0"]) + "->" + str(self.transition_table["term_1"])

    def plottable_table(self, x_off_0=0.5, x_off_1=0.5, **kwargs):
        params = {"scale_splitting"  : 1.0,
                  "color"            : None}
        for param in kwargs.keys():
            params[param] = kwargs[param]

        table_0 = self.level_0.plottable_table(scale_splitting = params["scale_splitting"])
        table_1 = self.level_1.plottable_table(scale_splitting = params["scale_splitting"])
        line_0 = table_0.loc[(table_0["m_F"] == self.m_F_0) & (table_0["F"] == self.F_0)]
        line_1 = table_1.loc[(table_1["m_F"] == self.m_F_1) & (table_1["F"] == self.F_1)]
        line_0 = line_0.drop(["color", "J_frac", "F_frac", "m_F_frac"], axis=1)
        line_1 = line_1.drop(["color", "J_frac", "F_frac", "m_F_frac"], axis=1)
        if line_0.empty or line_1.empty:
            raise ValueError("The selected quantum numbers don't yield a transition. Check that they exist in the specified levels")
        line_0.insert(9, "x", line_0["x0"]+x_off_0*(line_0["x1"]-line_0["x0"]))
        line_1.insert(9, "x", line_1["x0"]+x_off_1*(line_1["x1"]-line_1["x0"]))
        line_0.columns = [str(col) + '_0' for col in line_0.columns]
        line_1.columns = [str(col) + '_1' for col in line_1.columns]
        line_0 = line_0.reset_index(drop=True)
        line_1 = line_1.reset_index(drop=True)

        transition_table = pd.concat([line_0, line_1], axis=1, ignore_index=False)

        def frequency_to_color(freq):
            """
            Takes a frequency (in Hertz for now) and outputs a hex color value based on the frequency
            It does this by interpolating between set points where I know that the colors should be.
            """
            c = 3e8
            wl = c / abs(float(freq)) * 1e9

            wavelengths =   [0, 200, 280, 300, 350, 400, 445, 475, 493, 510, 570, 650, 780, 1000, 1500]
            rs =            [0, 0,   0,   50,  120, 80 , 110, 0  , 40 , 40 , 255, 235, 90 , 50 , 0]
            gs =            [0, 0,   0,   20,  120, 40 , 0  , 0  , 255, 255, 230, 30 , 30 , 20 , 0]
            bs =            [0, 0,   0,   100, 255, 120, 255, 255, 250, 0  , 0  , 30 , 30 , 20 , 0]

            rf = interpolate.interp1d(wavelengths, rs)
            gf = interpolate.interp1d(wavelengths, gs)
            bf = interpolate.interp1d(wavelengths, bs)

            if wl < 1500:
                r, g, b = rf(wl), gf(wl), bf(wl)
            else:
                r, g, b = 0, 0, 0

            R = "{:02x}".format(int(r))
            G = "{:02x}".format(int(g))
            B = "{:02x}".format(int(b))

            return "#" + R + G + B

        delta_l = abs(transition_table["level_0"] - transition_table["level_1"])
        if params["color"] is None:
            params["color"] = frequency_to_color(delta_l*1e12)
        wavelength = 299792.458/delta_l

        transition_table["delta_l"] = [delta_l]
        transition_table["color"] = [params["color"]]
        transition_table["wavelength"] = [wavelength]
        return transition_table

    def get_type(self):
        # TODO: finish this method sometime when I have time to actually fiure out what all the relevant selection rules are
        d_m = self.m_F_1 - self.m_F_0
        d_F = self.F_1 - self.F_0
        d_J = self.J_1 - self.J_0

        if abs(d_J) <= 1 and not (self.J_0 == 0  and self.J_1 == 0):
            if abs(d_F) <= 1:
                if self.parity_0 != self.parity_1:
                    transition_type = "E1"

        return "dm: {}, dF: {}, dJ:{}".format(d_m, d_F, d_J)


class Atom:
    def __init__(self, name, levels=(), transitions=()):
        self.levels = {}
        for level in levels:
            if level.name in self.levels.keys():
                raise Exception("{} is already a level in this atom. ".format(level.name))
            self.levels[level.name] = level

        self.transitions = {}
        for transition in transitions:
            self.transitions[transition.name] = transition
        self.name = name

    def add_level(self, *levels):
        for level in levels:
            self.levels[level.name] = level

    def remove_level(self, mode="keys", *levels):
        if mode == "keys":
            for key in levels:
                self.levels.pop(key)
        elif mode == "levels":
            for level in levels:
                self.levels.pop(level.name)
        else:
            raise ValueError("Unrecognized mode")

    def add_transition(self, *transitions):
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

    def list_levels(self):
        return self.levels

    def list_transitions(self):
        return self.transitions

    def line_table(self, **kwargs):
        table = DataFrame(columns=["configuration", "term", "level",
                                   "J", "F", "m_F", "J_frac", "F_frac", "m_F_frac",
                                   "color", "y0", "hf", "z",
                                   "y", "x0", "x1"])
        for level in self.levels.values():
            table = table.append(level.plottable_table(**kwargs), ignore_index=True)
        return table

    def transition_table(self, **kwargs):
        table = DataFrame(columns = [u'F_0', u'J_0', u'configuration_0', u'hf_0', u'level_0', u'm_F_0',
       u'term_0', u'x0_0', u'x1_0', u'y_0', u'y0_0', u'z_0', u'F_1', u'J_1',
       u'configuration_1', u'hf_1', u'level_1', u'm_F_1', u'term_1', u'x0_1',
       u'x1_1', u'y_1', u'y0_1', u'z_1', u'delta_l', u'color', u'wavelength'])
        for transition in self.transitions.values():
            table = table.append(transition.plottable_table(**kwargs), ignore_index=True)
        return table