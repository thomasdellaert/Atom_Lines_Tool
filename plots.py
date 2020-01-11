# TODO: Horizontal line plot given two levels (like in page 92 of Foot, or like in Barends/Maleki)
# TODO: Lorentzian plot
# TODO: King plot


from bokeh.io import show
from bokeh.plotting import figure, ColumnDataSource
import bokeh.models as models
from bokeh.layouts import row, column
from scipy import interpolate


class Grotrian:
    def __init__(self, atom, hf=True, zeeman=True):
        self.hf = hf
        self.zeeman = zeeman
        self.atom = atom
        self.plot_line_table = DataFrame(columns=["configuration", "term", "level",
                                                  "J", "F", "m_F", "J_frac", "F_frac", "m_F_frac",
                                                  "color", "y0", "hf", "z",
                                                  "y", "x0", "x1"])
        self.plot_transition_table = DataFrame(columns=["F_0", "J_0", "configuration_0", "hf_0", "level_0", "m_F_0",
                                                        "term_0", "x0_0", "x1_0", "y_0", "y0_0", "z_0", "F_1", "J_1",
                                                        "configuration_1", "hf_1", "level_1", "m_F_1", "term_1", "x0_1",
                                                        "x1_1", "y_1", "y0_1", "z_1", "delta_l", "color", "wavelength"])

    def level_table(self, level, width=1.0, sublevel_spacing=0.03, scale_splitting=1.0, override_position=False, 
                    offset_position=(0.0, 0.0), color="black"):
        table = DataFrame(columns=["configuration", "term", "level",
                                   "J", "F", "m_F", "J_frac", "F_frac", "m_F_frac",
                                   "color", "y0", "hf", "z",
                                   "y", "x0", "x1", "name"])
        if not override_position:
            x0 = level.J + offset_position[0] - width / 2
            y0 = level.level + offset_position[1]
        else:
            x0, y0 = override_position
            
        scale = scale_splitting
        
        if not self.hf:
            x1 = x0 + width
            line = DataFrame(data={"configuration": [level.configuration], "term": [level.term], "level": level.level,
                                   "J": [level.J], "F": [None], "m_F": [None],
                                   "J_frac": [term_frac(level.J)], "F_frac": [None], "m_F_frac": [None],
                                   "color": [color], "y0": [y0], "hf": [0.0], "z": [0.0],
                                   "y": [y0], "x0": [x0], "x1": [x1], "name": [level.name]})
            table = table.append(line, ignore_index=True)
        else:
            for F in level.Fs:
                if not self.zeeman:
                    hf = level.hf_shifts[F]
                    y = y0 + hf * scale
                    e_level = level.hf_levels[F]
                    x1 = x0 + width
                    line = DataFrame(
                        data={"configuration": [level.configuration], "term": [level.term], "level": e_level,
                              "J": [level.J], "F": [F], "m_F": [None],
                              "J_frac": [term_frac(level.J)], "F_frac": [term_frac(F)], "m_F_frac": [None],
                              "color": [color], "y0": [y0], "hf": [hf], "z": [0.0],
                              "y": [y], "x0": [x0], "x1": [x1], "name": [level.name]})
                    table = table.append(line, ignore_index=True)
                else:
                    delta = sublevel_spacing
                    F_max = max(level.Fs)
                    wd = (width - delta*2*F_max)/(2*F_max+1)
                    for m_F in level.z_shifts[F].keys():
                        z = level.z_shifts[F][m_F]
                        hf = level.hf_shifts[F]
                        y = y0 + hf*scale + z*scale
                        e_level = level.z_levels[F][m_F]
                        x = x0 + (m_F + F_max) * (wd + delta)
                        x1 = x + wd
                        line = DataFrame(
                            data={"configuration": [level.configuration], "term": [level.term], "level": e_level,
                                  "J": [level.J], "F": [F], "m_F": [m_F],
                                  "J_frac": [term_frac(level.J)], "F_frac": [term_frac(F)],
                                  "m_F_frac": [term_frac(m_F)],
                                  "color": [color], "y0": [y0], "hf": [hf], "z": [z],
                                  "y": [y], "x0": [x], "x1": [x1], "name": [level.name]})
                        table = table.append(line, ignore_index=True)
        return table
    
    def transition_table(self, transition, x_off_0=0.5, x_off_1=0.5, scale_splitting=1.0, color=None):
        # TODO: figure out a way to arbitrary kwargs into the table
        table_0 = self.level_table(transition.level_0, scale_splitting=scale_splitting)
        table_1 = self.level_table(transition.level_1, scale_splitting=scale_splitting)
        line_0 = table_0.loc[(table_0["m_F"] == transition.m_F_0) & (table_0["F"] == transition.F_0)]
        line_1 = table_1.loc[(table_1["m_F"] == transition.m_F_1) & (table_1["F"] == transition.F_1)]
        line_0 = line_0.drop(["color", "J_frac", "F_frac", "m_F_frac"], axis=1)
        line_1 = line_1.drop(["color", "J_frac", "F_frac", "m_F_frac"], axis=1)
        if line_0.empty or line_1.empty:
            raise ValueError("The selected quantum numbers don't yield a transition. "
                             "Check that they exist in the specified levels")
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

            wavelengths = [
                  280, 300, 350, 400, 445, 475, 493, 510, 570, 650, 780, 1000, 1500]
            rs = [0,   50,  120, 80,  110, 0,   40,  40,  255, 235, 90,  50,   0]
            gs = [0,   20,  120, 40,  0,   0,   255, 255, 230, 30,  30,  20,   0]
            bs = [0,   100, 255, 120, 255, 255, 250, 0,   0,   30,  30,  20,   0]

            rf = interpolate.interp1d(wavelengths, rs)
            gf = interpolate.interp1d(wavelengths, gs)
            bf = interpolate.interp1d(wavelengths, bs)

            if 280 < wl < 1500:
                r, g, b = rf(wl), gf(wl), bf(wl)
            else:
                r, g, b = 0, 0, 0

            return "#{:02x}{:02x}{:02x}".format(int(r), int(g), int(b))

        delta_l = abs(transition_table["level_0"] - transition_table["level_1"])
        if color is None:
            color = frequency_to_color(delta_l*1e12)
        wavelength = 299792.458/delta_l

        transition_table["delta_l"] = [delta_l]
        transition_table["color"] = [color]
        transition_table["wavelength"] = [wavelength]
        transition_table["name"] = [transition.name]

        return transition_table

    def add_level(self, levels, **kwargs):
        for level in levels:
            self.plot_line_table = self.plot_line_table.append(self.level_table(level, **kwargs))
        return self.plot_line_table

    def add_transition(self, transitions, **kwargs):
        for transition in transitions:
            self.plot_transition_table = self.plot_transition_table.append(self.transition_table(transition, **kwargs))
        return self.plot_transition_table

    # TODO: remove_level, remove_transition, level_style(), transition_style()

    def build_figure(self, dimensions=(800, 1000), y_range=(-1e2, 1.3e3), x_range=(0, 4), title=None, scale_splitting=1):
        p = figure(title=title, plot_width=dimensions[0], plot_height=dimensions[1], y_range=y_range, x_range=x_range)
        line_source = ColumnDataSource(self.plot_line_table)
        arrow_source = ColumnDataSource(self.plot_transition_table)
        lines = p.segment(x0="x0", y0="y", x1="x1", y1="y",
                          color="color", source=line_source)
        arrows = p.segment(x0="x_0", y0="y_0", x1="x_1", y1="y_1",
                           color="color", line_width=3, source=arrow_source)
        # TODO: Maybe make the arrows arrow-y? Might be more trouble than it's worth

        hover_lines = models.HoverTool(tooltips=[("Term", "@name F=@F_frac, m_F=@m_F_frac"),
                                                 ("Level", "@level{0.000000}")], renderers=[lines])
        hover_arrows = models.HoverTool(tooltips=[("Name", "@name"), ("Frequency", "@delta_l{0.000000} THz"),
                                                  ("Wavelength", "@wavelength{0.00} nm")], renderers=[arrows])
        p.add_tools(hover_lines)
        p.add_tools(hover_arrows)

        scale_slider = models.Slider(start=1, end=10000, value=scale_splitting, step=10, title="Scaling")
        b_field_slider = models.Slider(start=0, end=100, value=5, step=0.01, title="B-field (G)")
        line_callback = models.CustomJS(args=dict(source=line_source, scale=scale_slider, b_field=b_field_slider),
                                        code="""
                var data = source.data;
                var b_field = b_field.value;
                var scale = scale.value;

                hf = data['hf'];
                z = data['z'];
                y = data['y'];
                y0 = data['y0'];
                level = data['level'];

                for (i=0; i < y.length; i++) {
                    y[i] = hf[i]*scale + z[i]*b_field*scale + y0[i];
                    level[i] = y0[i] + hf[i] + z[i]*b_field;
                };       
                source.change.emit();
            """)
        arrow_callback = models.CustomJS(args=dict(source=arrow_source, scale=scale_slider, b_field=b_field_slider),
                                         code="""
                var data = source.data;
                var b_field = b_field.value;
                var scale = scale.value;

                hf0 = data['hf_0'];
                hf1 = data['hf_1'];
                z0 = data['z_0'];
                z1 = data['z_1'];
                y0 = data['y_0'];
                y1 = data['y_1'];
                y01 = data['y0_1'];
                y00 = data['y0_0'];
                level0 = data['level_0'];
                level1 = data['level_1'];
                delta_l = data['delta_l'];

                for (i=0; i < y0.length; i++) {
                    y0[i] = hf0[i]*scale + z0[i]*b_field*scale + y00[i];
                    y1[i] = hf1[i]*scale + z1[i]*b_field*scale + y01[i];
                    level0[i] = y00[i] + hf0[i] + z0[i]*b_field;
                    level1[i] = y01[i] + hf1[i] + z1[i]*b_field;
                    delta_l[i] = level1[i] - level0[i]
                };
                source.change.emit();
            """)
        scale_slider.js_on_change('value', line_callback)
        scale_slider.js_on_change('value', arrow_callback)
        b_field_slider.js_on_change('value', line_callback)
        b_field_slider.js_on_change('value', arrow_callback)

        show(row(p, column(scale_slider, b_field_slider)))

        return row(p, column(scale_slider, b_field_slider))


if __name__ == "__main__":
    from atom_library import *

    pd.set_option('display.max_rows', 500)
    pd.set_option('display.max_columns', 500)
    pd.set_option('display.width', 1000)

    atom = Yb_171

    g = Grotrian(atom)
    levels = atom.levels.values()
    mods = []
    defs = []
    for i in range(len(levels)):
        if levels[i].tamper:
            mods.append(levels[i])
        else:
            defs.append(levels[i])
    g.add_level(defs, color="lightgray")
    g.add_level(mods, color="black")
    g.add_transition(atom.transitions.values())

    g.build_figure()
