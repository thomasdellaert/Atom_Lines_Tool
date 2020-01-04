# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 13:57:39 2019

@author: Thomas Dellaert

"""
# TODO: Horizontal line plot given two levels (like in page 92 of Foot, or like in Barends/Maleki)
# TODO: Lorentzian plot
# TODO: King plot


from bokeh.io import output_file, show
from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import HoverTool, Label, LabelSet, Panel, Tabs, Slider, CustomJS
from bokeh.layouts import row, column
from scipy import interpolate
from pandas import DataFrame
import pandas as pd


def build_reference_plot(df, nstates):
    """
    Generate a reference plot used to pick out the levels that you care about.
    Those will receive their own EnergyLevel objects, with more versatile plotting options.
    Requires a datafile produced from NIST by the term parser.
    """

    p = figure(plot_width=600, plot_height=800, y_range=(-1e3, df["Level (cm-1)"][nstates] * 1.1))

    source = ColumnDataSource(data=dict(
        x0=df["J"][:nstates],
        x1=df["J"][:nstates] + 1,
        y0=df["Level (THz)"][:nstates],
        term=df["Term"][:nstates], ))

    labels = LabelSet(x="x0", y="y0", text="term", level="glyph", source=source,
                      x_offset=-5, y_offset=0, text_font_size="8pt")

    hover = HoverTool(tooltips=[("index", "$index"), ("Term", "@term"), ("Level", "@y0")])
    p.add_tools(hover)
    p.segment("x0", "y0", "x1", "y0", line_width=2, source=source)
    p.add_layout(labels)

    show(p)


# TODO: Make a Grotrian class with methods like g.transition, g.level, g.label_level, g.label_splitting. g.build_figure.
# TODO: migrate any references to plotting from the atoms/levels/transitions and move it into the plotting classes.

class Grotrian:
    def __init__(self, atom, hf=True, zeeman=True):
        self.hf=hf
        self.zeeman=zeeman
        self.atom = atom
        self.plot_line_table = DataFrame(columns=["configuration", "term", "level",
                                   "J", "F", "m_F", "J_frac", "F_frac", "m_F_frac",
                                   "color", "y0", "hf", "z",
                                   "y", "x0", "x1"])
        self.plot_transition_table = DataFrame(columns = ["F_0", "J_0", "configuration_0", "hf_0", "level_0", "m_F_0",
                                                     "term_0", "x0_0", "x1_0", "y_0", "y0_0", "z_0", "F_1", "J_1",
                                                     "configuration_1", "hf_1", "level_1", "m_F_1", "term_1", "x0_1",
                                                     "x1_1", "y_1", "y0_1", "z_1", "delta_l", "color", "wavelength"])

    def level_table(self, level, width=1.0, sublevel_spacing=0.03, scale_splitting=1.0, override_position=False, 
                    offset_position=(0.0, 0.0), color="black"):
        table = DataFrame(columns=["configuration", "term", "level",
                                   "J", "F", "m_F", "J_frac", "F_frac", "m_F_frac",
                                   "color", "y0", "hf", "z",
                                   "y", "x0", "x1"])
        if not override_position:
            x0 = level.J + offset_position[0] - width / 2
            y0 = level.level + offset_position[1]
        else:
            x0, y0 = override_position
            
        scale = scale_splitting
        
        if not self.hf:
            x1 = x0 + width
            line = DataFrame(data={"configuration": [level.configuration], "term": [level.term], "level": level.level,
                                   "J"            : [level.J], "F": [None], "m_F": [None],
                                   "J_frac"       : [level.J_frac], "F_frac": [None], "m_F_frac": [None],
                                   "color"        : [color], "y0": [y0], "hf": [0.0], "z": [0.0],
                                   "y"            : [y0], "x0": [x0], "x1": [x1]})
            table = table.append(line, ignore_index=True)
        else:
            for F in level.Fs:
                if not self.zeeman:
                    hf = level.hf_shifts[F]
                    y = y0 + hf * scale
                    level = level.hf_levels[F]
                    x1 = x0 + width
                    line = DataFrame(data={"configuration": [level.configuration], "term": [level.term], "level": level,
                                           "J"            : [level.J], "F": [F], "m_F": [None],
                                           "J_frac"       : [level.J_frac], "F_frac": [level.F_frac], "m_F_frac": [None],
                                           "color"        : [color], "y0": [y0], "hf": [hf], "z": [0.0],
                                           "y"            : [y], "x0": [x0], "x1": [x1]})
                    table = table.append(line, ignore_index=True)
                else:
                    delta = sublevel_spacing
                    F_max = max(level.Fs)
                    wd = (width - delta*2*F_max)/(2*F_max+1)
                    for m_F in level.z_shifts[F].keys():
                        z = level.z_shifts[F][m_F]
                        hf = level.hf_shifts[F]
                        y = y0 + hf*scale + z*scale
                        level = level.z_levels[F][m_F]
                        x = x0 + (m_F + F_max) * (wd + delta)
                        x1 = x + wd
                        line = DataFrame(
                            data={"configuration": [level.configuration], "term": [level.term], "level": level,
                                  "J"            : [level.J], "F": [F], "m_F": [m_F],
                                  "J_frac"       : [level.J_frac], "F_frac": [level.F_frac], "m_F_frac": [level.m_F_frac],
                                  "color"        : [color], "y0": [y0], "hf": [hf], "z": [z],
                                  "y"            : [y], "x0": [x], "x1": [x1]})
                        table = table.append(line, ignore_index=True)
        return table
    
    def transition_table(self, transition, x_off_0=0.5, x_off_1=0.5, scale_splitting=1.0, color=None):
        table_0 = self.level_table(transition.level_0, scale_splitting=scale_splitting)
        table_1 = self.level_table(transition.level_1, scale_splitting=scale_splitting)
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
        if color is None:
            color = frequency_to_color(delta_l*1e12)
        wavelength = 299792.458/delta_l

        transition_table["delta_l"] = [delta_l]
        transition_table["color"] = [color]
        transition_table["wavelength"] = [wavelength]
        return transition_table

    def add_level(self, levels, **kwargs):
        for level in levels:
            self.plot_line_table = self.plot_line_table.append(self.level_table(level, **kwargs))
        return self.plot_line_table()

    def add_transition(self, transitions, **kwargs):
        for transition in transitions:
            self.plot_transition_table = self.plot_transition_table(self.transition_table(transition, **kwargs))
        return self.plot_transition_table()

    def build_figure(self, dimensions=(800, 1000), y_range=(-1e2, 1e3), title=None, zeeman=True, hf=True, scale_splitting=3000):
        p = figure(title=title, plot_width=dimensions[0], plot_height=dimensions[1], y_range=y_range, x_range=(0, 4))
        line_source = ColumnDataSource(self.plot_line_table)
        arrow_source = ColumnDataSource(self.plot_transition_table)
        lines = p.segment(x0="x0", y0="y", x1="x1", y1="y",
                          color="color", source=line_source)
        arrows = p.segment(x0="x_0", y0="y_0", x1="x_1", y1="y_1",
                           color="color", line_width=3, source=arrow_source)
        # TODO: Maybe make the arrows arrow-y? Might be more trouble than it's worth

        hover_lines = HoverTool(tooltips=[("Term", "@term @J_frac F=@F_frac, m_F=@m_F_frac"),
                                          ("Level", "@level{0.000000}")], renderers=[lines])
        hover_arrows = HoverTool(tooltips=[("Frequency", "@delta_l{0.000000} THz"),
                                           ("Wavelength", "@wavelength{0.00} nm")], renderers=[arrows])
        p.add_tools(hover_lines)
        p.add_tools(hover_arrows)

        scale_slider = Slider(start=1, end=10000, value=scale_splitting, step=10, title="Scaling")
        b_field_slider = Slider(start=0, end=100, value=5, step=0.01, title="B-field (G)")

        line_callback = CustomJS(args=dict(source=line_source, scale=scale_slider, b_field=b_field_slider), code="""
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

        arrow_callback = CustomJS(args=dict(source=arrow_source, scale=scale_slider, b_field=b_field_slider), code="""
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

def build_grotrian_diagram(atom, dimensions=(800, 1000), y_range=(-1e2, 1e3),
                           title=None, zeeman=True, hf=True, scale=3000, show_plot=True):
    """
    Makes a grotrian diagram of the given atom

    dimensions (tuple) - dimensions of the plot
    y_range (tuple)    - the portion of the y-axis to display
    x_range (tuple)    - the portion of the x-axis to display
    name (string)      - the name of the plot. Will show up as a title
    zeeman (bool)      - whether to display zeeman sublevels
    hf (bool)          - whether to display hyperfine sublevels
    scale (int)        - initial scaling for the hyperfine structure display. Can change later w/ a slider
    arrows (list)      - a list of tuples: (start_index, end_index)
    show_plot (bool)   - Whether to show the plot when done. Set to false for passing to other things
    """

    p = figure(title=title, plot_width=dimensions[0], plot_height=dimensions[1], y_range=y_range, x_range=(0, 4))

    line_table = atom.line_table(hyperfine=hf, zeeman=zeeman, scale_splitting=scale)
    line_source = ColumnDataSource(line_table)

    arrow_table = atom.transition_table(scale_splitting=scale)
    arrow_source = ColumnDataSource(arrow_table)

    lines = p.segment(x0="x0", y0="y", x1="x1", y1="y",
                      color="color", source=line_source)
    arrows = p.segment(x0="x_0", y0="y_0", x1="x_1", y1="y_1",
                       color="color", line_width=3, source=arrow_source)
    # TODO: Maybe make the arrows arrow-y? Might be more trouble than it's worth

    hover_lines = HoverTool(tooltips=[("Term", "@term @J_frac F=@F_frac, m_F=@m_F_frac"),
                                      ("Level", "@level{0.000000}")], renderers=[lines])
    hover_arrows = HoverTool(tooltips=[("Frequency", "@delta_l{0.000000} THz"),
                                       ("Wavelength", "@wavelength{0.00} nm")], renderers=[arrows])
    p.add_tools(hover_lines)
    p.add_tools(hover_arrows)

    scale_slider = Slider(start=1, end=10000, value=scale, step=10, title="Scaling")
    b_field_slider = Slider(start=0, end=100, value=5, step=0.01, title="B-field (G)")

    line_callback = CustomJS(args=dict(source=line_source, scale=scale_slider, b_field=b_field_slider), code="""
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

    arrow_callback = CustomJS(args=dict(source=arrow_source, scale=scale_slider, b_field=b_field_slider), code="""
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

    if show_plot:
        show(row(p, column(scale_slider, b_field_slider)))

    return row(p, column(scale_slider, b_field_slider))


def build_final_plots(atoms, **kwargs):
    """
    Builds a mulltitab grotrian plot out of multiple atom objects. Not sure how useful this is.
    """
    panels = []
    for atom in atoms:
        panels.append(Panel(child=build_grotrian_diagram(atom, **kwargs), title=atom.name))
    tabs = Tabs(tabs=panels)
    show(tabs)

if __name__ == "__main__":
    from Atom_Library import Yb_171, Yb_173, Yb_174

    g = Grotrian(Yb_173)
    g.add_level(Yb_173.list_levels().values())
    g.add_transition(Yb_173.list_transitions().values())

    g.build_figure()

    build_grotrian_diagram(Yb_173, show_plot=True)

    output_file("Grotrian.html")
