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

    build_grotrian_diagram(Yb_173, show_plot=True)

    output_file("Grotrian.html")
