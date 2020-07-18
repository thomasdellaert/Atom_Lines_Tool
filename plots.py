# TODO: Lorentzian plot
import copy

from bokeh.io import show
from bokeh.plotting import figure, ColumnDataSource
import bokeh.models as models
from bokeh.layouts import row, column, gridplot
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
    def __init__(self, *levels):
        self.levels = levels
        self.plot_line_table = DataFrame(columns=["configuration", "term", "level",
                                                  "J", "F", "m_F", "J_frac", "F_frac", "m_F_frac",
                                                  "color", "y0", "hf", "z",
                                                  "y", "x0", "x1"])
        self.plot_arrow_table = DataFrame(columns=["F_0", "hf_0", "level_0", "m_F_0", "term_0", "x0_0", "x1_0", "y_0", "y0_0", "z_0",
                                                   "F_1", "hf_1", "level_1", "m_F_1", "term_1", "x0_1", "x1_1", "y_1", "y0_1", "z_1",
                                                   "delta_l", "wavelength"])
        for level in self.levels:
            self.plot_line_table = self.plot_line_table.append(self.level_table(level[0], level[1]))

        self.plot_arrow_table = self.arrow_table(self.plot_line_table, spacing='physical')

    def level_table(self, level, F, width=0.01, b_field=1, scale_splitting_hf=1.0, scale_splitting_z=1.0, color="black", y0=0):
        table = level.data_table(hf=True, zeeman=True, b_field=b_field)
        table = table[table['F']==F].reset_index()
        table['color'] = color
        table['name'] = level.name

        table['y0'] = y0
        table['x0'] = 0
        table['x1'] = width

        table['y'] = table['y0'] + table['hf'] * scale_splitting_hf + table['z'] * scale_splitting_z
        table['level'] = table['y0'] + table['hf'] + table['z'] * b_field

        table.loc[table['m_F'] == 0, 'hflabel'] = "F="+table["F_frac"]
        table.loc[table['m_F'] != 0, 'hflabel'] = ""
        table['hflx'] = width+0.05
        table['hfly'] = table['y']

        table['zlabel'] = table['m_F_frac']
        table['zlx'] = table['x0']
        table['zly'] = table['y']

        table['tlabel'] = ""
        table.loc[[0], 'tlabel'] = table['term'] + table['J_frac']
        table['tlx'] = table['x0']
        table['tly'] = (max(table['y']) + min(table['y'])) / 2

        return table

    def arrow_table(self, level_table, spacing="physical"):
        table = DataFrame(columns=["F_0", "hf_0", "m_F_0", "y_0", "y0_0", "z_0",
                                   "F_1", "hf_1", "m_F_1", "y_1", "y0_1", "z_1",
                                   "level_0", "level_1", "delta_l", "x", "J0", "J1"])
        for index0, sl0 in level_table.iterrows():
            for index1, sl1 in level_table.iterrows():
                m_F_0, m_F_1, F_0, F_1 = sl0['m_F'], sl1['m_F'], sl0['F'], sl1['F']
                term_0, term_1 = sl0['term'], sl1['term']
                if F_0 != F_1 or term_0 != term_1:
                    if abs(m_F_0-m_F_1) <= 1 and index0 < index1:
                        line = DataFrame(data={'F_0': [F_0], 'hf_0': [sl0['hf']], 'm_F_0': [m_F_0],
                                               'y_0': [sl0['y']], 'y0_0': [sl0['y0']], 'z_0': [sl0['z']],
                                               'F_1': [F_1], 'hf_1': [sl1['hf']], 'm_F_1': [m_F_1],
                                               'y_1': [sl1['y']], 'y0_1': [sl1['y0']], 'z_1': [sl1['z']],
                                               'level_0': sl0['level'], 'level_1': sl1['level'],
                                               'delta_l': abs(sl0['level']-sl1['level']),
                                               'x': abs(sl0['level']-sl1['level']),
                                               'J_0': sl0['J'], 'J_1': sl1['J']})
                        table = table.append(line, ignore_index=True)

        #TODO: Color based on pi vs sigma transition

        def space_out_lines(series, spacing):
            s = copy.deepcopy(series)
            s.sort_values(ascending=True, inplace=True)
            for i in range(len(s) - 1):
                if s[i + 1] - s[i] < spacing:
                    s[i] -= spacing - (s[i + 1] - s[i]) / 2
                    s[i + 1] += spacing - (s[i + 1] - s[i]) / 2
                    s = space_out_lines(s, spacing)
            return s.sort_index()

        if spacing != 'physical':
            table['x'] = space_out_lines(table['x'], spacing)

        return table

    def build_figure(self, dimensions=(1600, 1000), title=None, scale_splitting_hf=1, scale_splitting_z=1, display=False, labels=[],
                     b_field_slider = models.Slider(start=0, end=5, value=1, step=0.0001, title="B-field (G)")):
        p = figure(title=title, plot_width=dimensions[0], plot_height=dimensions[1],
                   y_range=(min(self.plot_line_table['y']), max(self.plot_line_table['y'])),
                   x_range=(min(self.plot_arrow_table['delta_l'])-1e-6, max(self.plot_arrow_table['delta_l'])+1e-6))
        line_source = ColumnDataSource(self.plot_line_table)
        arrow_source = ColumnDataSource(self.plot_arrow_table)
        print "plotting levels"
        lines = p.segment(x0="x0", y0="y", x1="x1", y1="y",
                          color="color", source=line_source)
        print "plotting transitions"
        arrows = p.segment(x0="delta_l", y0="y_0", x1="delta_l", y1="y_1",
                           line_width=3, source=arrow_source)

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
        hover_arrows = models.HoverTool(tooltips=[("Name", "@name"), ("Frequency", "@delta_l{0.000000} THz")], renderers=[arrows])
        p.add_tools(hover_lines)
        p.add_tools(hover_arrows)

        print "applying sliders"
        scale_hf_slider = models.Slider(start=0.1, end=2, value=scale_splitting_hf, step=0.01, title="HF Scaling")
        scale_z_slider = models.Slider(start=1, end=10, value=scale_splitting_z, step=0.1, title="Zeeman Scaling")
        b_field_slider = b_field_slider
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
        arrow_callback = models.CustomJS(args=dict(source=arrow_source, hf_scale=scale_hf_slider, z_scale=scale_z_slider, b_field=b_field_slider),
                                         code="""
                       var data = source.data;
                       var b_field = b_field.value;
                       var hf_scale = hf_scale.value;
                       var z_scale = z_scale.value;

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
                           y0[i] = hf0[i]*hf_scale + z0[i]*b_field*z_scale + y00[i];
                           y1[i] = hf1[i]*hf_scale + z1[i]*b_field*z_scale + y01[i];
                           level0[i] = hf0[i] + z0[i]*b_field + y00[i];
                           level1[i] = hf1[i] + z1[i]*b_field + y01[i];
                           delta_l[i] = Math.abs(level1[i] - level0[i]);
                       };
                       source.change.emit();
                   """)
        scale_hf_slider.js_on_change('value', line_callback)
        scale_z_slider.js_on_change('value', line_callback)
        b_field_slider.js_on_change('value', line_callback)

        scale_hf_slider.js_on_change('value', arrow_callback)
        scale_z_slider.js_on_change('value', arrow_callback)
        b_field_slider.js_on_change('value', arrow_callback)

        if display:
            print "displaying Hyperfine diagram"
            show(row(p, column(scale_hf_slider, scale_z_slider, b_field_slider)))

        return p, column(scale_hf_slider, scale_z_slider, b_field_slider)


class Lorentzian_plot(HF_plot):
    def build_figure(self, linewidth=5e-9, dimensions=(1600, 400), title=None, scale_splitting_hf=1, scale_splitting_z=1, display=False, x_range=None, labels=[],
                     b_field_slider = models.Slider(start=0, end=5, value=1, step=0.0001, title="B-field (G)")):
        import numpy as np
        from math import pi
        from transition_strengths import rel_transiton_strength
        if x_range is None:
            x_range=(min(self.plot_arrow_table['delta_l'])-1e-6, max(self.plot_arrow_table['delta_l'])+1e-6)
        p = figure(title=title, plot_width=dimensions[0], plot_height=dimensions[1], x_range=x_range)

        xaxis = np.linspace(min(self.plot_arrow_table['delta_l'])-1e-6, max(self.plot_arrow_table['delta_l'])+1e-6, 10000)

        lines = []
        strengths = []
        for index, transition in self.plot_arrow_table.iterrows():
            m0, m1 = transition['m_F_0'], transition['m_F_1']
            F0, F1 = transition['F_0'], transition['F_1']
            J1= transition['J_1']
            I = self.levels[0][0].I
            q = m1-m0

            line = (1/(2*pi))*linewidth/((xaxis-transition['delta_l'])**2+linewidth**2/4)
            strength = rel_transiton_strength(I, q, J1, F0, m0, F1, m1)
            lines.append(line*strength)
            strengths.append(strength)

        self.plot_arrow_table['strength'] = strengths
        totalline = np.sum(lines)

        transition_data = ColumnDataSource(data={
        'hf_0': self.plot_arrow_table['hf_0'],
        'hf_1': self.plot_arrow_table['hf_1'],
        'z_0': self.plot_arrow_table['z_0'],
        'z_1': self.plot_arrow_table['z_1'],
        'y_0': self.plot_arrow_table['y_0'],
        'y_1': self.plot_arrow_table['y_1'],
        'y0_1': self.plot_arrow_table['y0_1'],
        'y0_0': self.plot_arrow_table['y0_0'],
        'level_0': self.plot_arrow_table['level_0'],
        'level_1': self.plot_arrow_table['level_1'],
        'delta_l': self.plot_arrow_table['delta_l'],
        'strength': self.plot_arrow_table['strength']
        })

        lorentz_table = DataFrame(data={'xaxis': xaxis, 'totalline': totalline})#, 'transition_data': transition_data})
        lorentz_source = ColumnDataSource(lorentz_table)

        print "plotting transitions"

        p.line(x='xaxis', y='totalline', source=lorentz_source)

        print "applying sliders"

        linewidth_slider = models.Slider(start=1e-9, end=1e-7, value=5e-9, step=1e-9, title="linewidth (THz)")
        lorentz_callback = models.CustomJS(
            args=dict(source=lorentz_source, b_field=b_field_slider, linewidth=linewidth_slider, transition_data=transition_data),
            code="""
                var linedata = source.data;
                var transition_data = transition_data.data;
                var b_field = b_field.value;
                var linewidth = linewidth.value;
                
                hf0 = transition_data['hf_0'];
                hf1 = transition_data['hf_1'];
                z0 = transition_data['z_0'];
                z1 = transition_data['z_1'];
                y01 = transition_data['y0_1'];
                y00 = transition_data['y0_0'];
                level0 = transition_data['level_0'];
                level1 = transition_data['level_1'];
                delta_l = transition_data['delta_l'];
                strength = transition_data['strength'];
                
                xaxis = linedata['xaxis'];
                totalline = linedata['totalline'];
                                
                for (i=0; i < hf0.length; i++) {
                    level0[i] = hf0[i] + z0[i]*b_field + y00[i];
                    level1[i] = hf1[i] + z1[i]*b_field + y01[i];
                    delta_l[i] = Math.abs(level1[i] - level0[i]);
                };
                
                for (i=0; i< xaxis.length; i++){
                    let tot = 0;
                    for (j=0; j<hf0.length; j++){
                        tot += strength[j]*(0.5/Math.PI)*linewidth/((xaxis[i]-delta_l[j])*(xaxis[i]-delta_l[j])+linewidth*linewidth/4);
                    };
                    totalline[i] = tot;
                };
                
                source.change.emit();
                """)

        b_field_slider.js_on_change('value', lorentz_callback)
        linewidth_slider.js_on_change('value', lorentz_callback)

        if display:
            print "displaying Lorentzian diagram"
            show(row(p, column(b_field_slider, linewidth_slider)))

        return p, column(b_field_slider, linewidth_slider)


if __name__ == "__main__":
    from atom_library import *
    from colors import uv_ir_lookup

    pd.set_option('display.max_rows', 500)
    pd.set_option('display.max_columns', 500)
    pd.set_option('display.width', 1000)

    atom = Yb_173
    def MakeGrotrian(atom):
        g = Grotrian(atom)
        levels = atom.levels.values()
        g.add_level(levels, color="black")
        g.add_transition(atom.transitions.values(), color=uv_ir_lookup)

        g.build_figure(display=True)#, labels=["hf", "zeeman", "term"])

    def MakeMixedPlot(level0, level1):
        b_field_slider = models.Slider(start=0, end=5, value=1, step=0.0001, title="B-field (G)")
        h = HF_plot(level0, level1)
        HFplot = h.build_figure(dimensions=(1600, 400), b_field_slider=b_field_slider)

        l = Lorentzian_plot(level0, level1)
        Lplot = l.build_figure(x_range=HFplot[0].x_range, b_field_slider=b_field_slider)

        show(row(gridplot([[Lplot[0]], [HFplot[0]]]), b_field_slider))

    def MakeLorentzPlot(level0, level1):
        l = Lorentzian_plot(level0, level1)
        l.build_figure(display=True, linewidth=1e-7)

    level0 = (atom.levels["2S1/2"], 2)
    level1 = (atom.levels["2S1/2"], 3)

    # MakeGrotrian(Yb_173)
    MakeMixedPlot(level0, level1)
    # MakeLorentzPlot(level0, level1)