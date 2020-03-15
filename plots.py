# TODO: Horizontal line plot given two levels (like in page 92 of Foot, or like in Barends/Maleki)
# TODO: Lorentzian plot


from bokeh.io import show
from bokeh.plotting import figure, ColumnDataSource
import bokeh.models as models
from bokeh.layouts import row, column
from colors import default_lookup

class Grotrian:
    def __init__(self, atom):
        self.atom = atom
        self.plot_line_table = DataFrame(columns=["configuration", "term", "level",
                                                  "J", "F", "m_F", "J_frac", "F_frac", "m_F_frac",
                                                  "color", "y0", "hf", "z",
                                                  "y", "x0", "x1",
                                                  "hflabel", "hflx", "hfly",
                                                  "zlabel", "zlx", "zly",
                                                  "tlx", "tly"])
        self.plot_transition_table = DataFrame(columns=["F_0", "J_0", "configuration_0", "hf_0", "level_0", "m_F_0",
                                                        "term_0", "x0_0", "x1_0", "y_0", "y0_0", "z_0", "F_1", "J_1",
                                                        "configuration_1", "hf_1", "level_1", "m_F_1", "term_1", "x0_1",
                                                        "x1_1", "y_1", "y0_1", "z_1", "delta_l", "color", "wavelength"])

    def level_table(self, level, width=1.0, sublevel_spacing=0.03, scale_splitting=1.0, override_position=False,
                    offset_position=(0.0, 0.0), color="black", zeeman=True, hf=True):
        table = level.data_table(hf=hf, zeeman=zeeman)
        table['color'] = color
        table['name'] = level.name

        if not override_position:
            x0 = level.J + offset_position[0] - width / 2
            y0 = level.level + offset_position[1]
        else:
            x0, y0 = override_position

        table['y0'] = y0

        if not zeeman:
            table['x0'] = x0
            table['x1'] = table['x0'] + width
            table['y'] = table['y0'] + table['hf'] * scale_splitting
            table['hflabel'] = "F="+str(table["F"])
            table['hflx'] = table['x1']+0.05
            table['hfly'] = table['y']
        else:
            delta = sublevel_spacing
            F_max = max(table['F'])
            wd = (width - delta*2*F_max)/(2*F_max+1)
            table['y'] = table['y0'] + table['hf'] * scale_splitting + table['z'] * scale_splitting
            table['x0'] = x0 + (table['m_F'] + F_max) * (wd + delta)
            table['x1'] = table['x0'] + wd

            table.loc[table['m_F'] == 0, 'hflabel'] = "F="+table["F_frac"]
            table.loc[table['m_F'] != 0, 'hflabel'] = ""
            table['hflx'] = x0+width+0.05
            table['hfly'] = table['y']

            table['zlabel'] = table['m_F_frac']
            table['zlx'] = table['x0']
            table['zly'] = table['y']

        table['tlabel'] = ""
        table.loc[[0], 'tlabel'] = table['term'] + table['J_frac']
        table['tlx'] = x0
        table['tly'] = (max(table['y']) + min(table['y'])) / 2

        return table
    
    def transition_table(self, transition, x_off_0=0.5, x_off_1=0.5, color=default_lookup):
        table_0 = self.level_table(transition.level_0)
        table_1 = self.level_table(transition.level_1)
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

        delta_l = abs(transition_table["level_0"] - transition_table["level_1"])
        wavelength = 299792.458/delta_l
        if type(color) is dict:
            wl = int(wavelength)
            if wl < min(color.keys()):
                wl = min(color.keys())
            elif wl > max(color.keys()):
                wl = max(color.keys())
            mycolor = color[wl]
        else:
            mycolor = color

        transition_table["delta_l"] = [delta_l]
        transition_table["color"] = [mycolor]
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

    def remove_level(self, level_names):
        for name in level_names:
            self.plot_line_table = self.plot_line_table[self.plot_line_table.name == name]

    def remove_transition(self, transition_names):
        for name in transition_names:
            self.plot_transition_table = self.plot_transition_table[self.plot_transition_table.name == name]

    def build_figure(self, dimensions=(800, 1000), y_range=(-1e2, 1.3e3), x_range=(-0.5, 4.5), title=None, scale_splitting=1,
                     labels=[], display=False):
        p = figure(title=title, plot_width=dimensions[0], plot_height=dimensions[1], y_range=y_range, x_range=x_range)
        line_source = ColumnDataSource(self.plot_line_table)
        arrow_source = ColumnDataSource(self.plot_transition_table)
        print "plotting levels"
        lines = p.segment(x0="x0", y0="y", x1="x1", y1="y",
                          color="color", source=line_source)
        print "plotting transitions"
        arrows = p.segment(x0="x_0", y0="y_0", x1="x_1", y1="y_1",
                           color="color", line_width=3, source=arrow_source)
        # TODO: Maybe make the arrows arrow-y? Might be more trouble than it's worth. Bokeh arrows are insufficient.

        if labels is not []:
            print "drawing labels"
        if "hf" in labels:
            print "  hf labels"
            hflabels = models.LabelSet(x="hflx", y="y", text="hflabel", level="glyph", source=line_source,
                                       text_baseline='middle', text_font_size="10pt")
            p.add_layout(hflabels)
        if "zeeman" in labels:
            print "  zeeman labels"
            zlabels = models.LabelSet(x="zlx", y="y", text="zlabel", level="glyph", source=line_source,
                                      text_font_size="8pt")
            p.add_layout(zlabels)
        if "term" in labels:
            print "  term labels"
            tlabels = models.LabelSet(x="tlx", y="tly", text="tlabel", level="glyph", source=line_source,
                                      text_align='right', text_font_style="bold", text_font_size="12pt",
                                      text_baseline="middle")
            p.add_layout(tlabels)

        print "applying hovertext"
        hover_lines = models.HoverTool(tooltips=[("Term", "@name F=@F_frac, m_F=@m_F_frac"),
                                                 ("Level", "@level{0.000 000 000}")], renderers=[lines])
        hover_arrows = models.HoverTool(tooltips=[("Name", "@name"), ("Frequency", "@delta_l{0.000000} THz"),
                                                  ("Wavelength", "@wavelength{0.00} nm")], renderers=[arrows])
        p.add_tools(hover_lines)
        p.add_tools(hover_arrows)

        print "applying sliders"
        scale_slider = models.Slider(start=1, end=10000, value=scale_splitting, step=10, title="Scaling")
        b_field_slider = models.Slider(start=0, end=20, value=0, step=0.001, title="B-field (G)")
        line_callback = models.CustomJS(args=dict(source=line_source, scale=scale_slider, b_field=b_field_slider),
                                        code="""
                var data = source.data;
                var b_field = b_field.value;
                var scale = scale.value;

                hf = data['hf'];
                z = data['z'];
                y = data['y'];
                y0 = data['y0'];
                hfly = data['hfly'];
                level = data['level'];

                for (i=0; i < y.length; i++) {
                    y[i] = hf[i]*scale + z[i]*b_field*scale + y0[i];
                    hfly[i] = y[i]
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

        if display:
            print "displaying Grotrian diagram"
            show(row(p, column(scale_slider, b_field_slider)))

        return row(p, column(scale_slider, b_field_slider))


class HF_plot:
    def __init__(self, level):
        self.level = level
        self.plot_line_table = DataFrame(columns=["configuration", "term", "level",
                                                  "J", "F", "m_F", "J_frac", "F_frac", "m_F_frac",
                                                  "color", "y0", "hf", "z",
                                                  "y", "x0", "x1"])
        self.plot_transition_table = DataFrame(columns=["F_0", "J_0", "configuration_0", "hf_0", "level_0", "m_F_0",
                                                        "term_0", "x0_0", "x1_0", "y_0", "y0_0", "z_0", "F_1", "J_1",
                                                        "configuration_1", "hf_1", "level_1", "m_F_1", "term_1", "x0_1",
                                                        "x1_1", "y_1", "y0_1", "z_1", "delta_l", "wavelength"])
        self.plot_line_table = self.level_table(level)

    def level_table(self, level, width=1.0, scale_splitting_hf=1.0, scale_splitting_z=1.0, color="black", zeeman=True, hf=True):
        table = level.data_table(hf=hf, zeeman=zeeman)
        table['color'] = color
        table['name'] = level.name

        x0 = 0
        y0 = 0

        table['y0'] = 0
        table['x0'] = 0
        table['x1'] = table['x0'] + width

        if not zeeman:

            table['y'] = table['y0'] + table['hf'] * scale_splitting_hf
            table['hflabel'] = "F="+str(table["F"])
            table['hflx'] = table['x1']+0.05
            table['hfly'] = table['y']
        else:
            table['y'] = table['y0'] + table['hf'] * scale_splitting_hf + table['z'] * scale_splitting_z

            table.loc[table['m_F'] == 0, 'hflabel'] = "F="+table["F_frac"]
            table.loc[table['m_F'] != 0, 'hflabel'] = ""
            table['hflx'] = 0+width+0.05
            table['hfly'] = table['y']

            table['zlabel'] = table['m_F_frac']
            table['zlx'] = table['x0']
            table['zly'] = table['y']

        table['tlabel'] = ""
        table.loc[[0], 'tlabel'] = table['term'] + table['J_frac']
        table['tlx'] = x0
        table['tly'] = (max(table['y']) + min(table['y'])) / 2

        return table

    # TODO: make the HF plot
    def build_figure(self, dimensions=(800, 1000), y_range=(-1e2, 1.3e3), x_range=(-0.5, 4.5), title=None, scale_splitting_hf=1, scale_splitting_z=1, display=False, labels=[]):
        p = figure(title=title, plot_width=dimensions[0], plot_height=dimensions[1], y_range=y_range, x_range=x_range)
        line_source = ColumnDataSource(self.plot_line_table)
        # arrow_source = ColumnDataSource(self.plot_transition_table)
        print "plotting levels"
        lines = p.segment(x0="x0", y0="y", x1="x1", y1="y",
                          color="color", source=line_source)
        # print "plotting transitions"
        # arrows = p.segment(x0="x_0", y0="y_0", x1="x_1", y1="y_1",
        #                    color="color", line_width=3, source=arrow_source)
        # if labels is not []:
        #     print "drawing labels"
        # if "hf" in labels:
        #     print "  hf labels"
        #     hflabels = models.LabelSet(x="hflx", y="y", text="hflabel", level="glyph", source=line_source,
        #                                text_baseline='middle', text_font_size="10pt")
        #     p.add_layout(hflabels)
        # if "zeeman" in labels:
        #     print "  zeeman labels"
        #     zlabels = models.LabelSet(x="zlx", y="y", text="zlabel", level="glyph", source=line_source,
        #                               text_font_size="8pt")
        #     p.add_layout(zlabels)
        # if "term" in labels:
        #     print "  term labels"
        #     tlabels = models.LabelSet(x="tlx", y="tly", text="tlabel", level="glyph", source=line_source,
        #                               text_align='right', text_font_style="bold", text_font_size="12pt",
        #                               text_baseline="middle")
        #     p.add_layout(tlabels)

        print "applying hovertext"
        hover_lines = models.HoverTool(tooltips=[("Term", "@name F=@F_frac, m_F=@m_F_frac"),
                                                 ("Level", "@level{0.000 000 000}")], renderers=[lines])
        # hover_arrows = models.HoverTool(tooltips=[("Name", "@name"), ("Frequency", "@delta_l{0.000000} THz"),
        #                                           ("Wavelength", "@wavelength{0.00} nm")], renderers=[arrows])
        p.add_tools(hover_lines)
        # p.add_tools(hover_arrows)

        print "applying sliders"
        scale_hf_slider = models.Slider(start=1, end=10000, value=scale_splitting_hf, step=10, title="HF Scaling")
        scale_z_slider = models.Slider(start=1, end=10000, value=scale_splitting_z, step=10, title="Zeeman Scaling")
        b_field_slider = models.Slider(start=0, end=20, value=0, step=0.001, title="B-field (G)")
        line_callback = models.CustomJS(args=dict(source=line_source, hf_scale=scale_hf_slider, z_scale=scale_z_slider, b_field=b_field_slider),
                                        code="""
                       var data = source.data;
                       var b_field = b_field.value;
                       var hf_scale = hf_scale.value;
                       var z_scale = z_scale.value

                       hf = data['hf'];
                       z = data['z'];
                       y = data['y'];
                       y0 = data['y0'];
                       hfly = data['hfly'];
                       level = data['level'];

                       for (i=0; i < y.length; i++) {
                           y[i] = hf[i]*hf_scale + z[i]*b_field*z_scale + y0[i];
                           hfly[i] = y[i]
                           level[i] = y0[i] + hf[i] + z[i]*b_field;
                       };       
                       source.change.emit();
                   """)
        # arrow_callback = models.CustomJS(args=dict(source=arrow_source, scale=scale_slider, b_field=b_field_slider),
        #                                  code="""
        #                var data = source.data;
        #                var b_field = b_field.value;
        #                var scale = scale.value;
        #
        #                hf0 = data['hf_0'];
        #                hf1 = data['hf_1'];
        #                z0 = data['z_0'];
        #                z1 = data['z_1'];
        #                y0 = data['y_0'];
        #                y1 = data['y_1'];
        #                y01 = data['y0_1'];
        #                y00 = data['y0_0'];
        #                level0 = data['level_0'];
        #                level1 = data['level_1'];
        #                delta_l = data['delta_l'];
        #
        #                for (i=0; i < y0.length; i++) {
        #                    y0[i] = hf0[i]*scale + z0[i]*b_field*scale + y00[i];
        #                    y1[i] = hf1[i]*scale + z1[i]*b_field*scale + y01[i];
        #                    level0[i] = y00[i] + hf0[i] + z0[i]*b_field;
        #                    level1[i] = y01[i] + hf1[i] + z1[i]*b_field;
        #                    delta_l[i] = level1[i] - level0[i]
        #                };
        #                source.change.emit();
        #            """)
        scale_hf_slider.js_on_change('value', line_callback)
        scale_z_slider.js_on_change('value', line_callback)
        # scale_slider.js_on_change('value', arrow_callback)
        b_field_slider.js_on_change('value', line_callback)
        # b_field_slider.js_on_change('value', arrow_callback)

        if display:
            print "displaying Hyperfine diagram"
            show(row(p, column(scale_hf_slider, scale_z_slider, b_field_slider)))

        return row(p, column(scale_hf_slider, scale_z_slider, b_field_slider))




if __name__ == "__main__":
    from atom_library import *
    from colors import uv_ir_lookup

    pd.set_option('display.max_rows', 500)
    pd.set_option('display.max_columns', 500)
    pd.set_option('display.width', 1000)

    atom = Yb_173
    #
    # g = Grotrian(atom)
    # levels = atom.levels.values()
    # g.add_level(levels, color="black")
    # g.add_transition(atom.transitions.values(), color=uv_ir_lookup)
    #
    # g.build_figure(display=True, labels=["hf", "zeeman", "term"])

    h = HF_plot(atom.levels["2S1/2"])
    h.build_figure(display=True)
