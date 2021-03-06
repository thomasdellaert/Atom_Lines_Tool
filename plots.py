"""
Make Bokeh plot representations of Atoms and Energy Levels

Classes:
    Grotrian: Makes a Grotrian diagram for a given Atom object
    HFPlot: Makes a plot of hyperfine sublevels and transitions for a given level or set of levels
    LorentzianPlot: Plots what the spectrum would look like for a given level or set of levels
"""

import copy
from pandas import DataFrame
from bokeh.plotting import figure, ColumnDataSource, show
import bokeh.models as models
from bokeh.layouts import row, column, gridplot
from colors import default_lookup


class Grotrian:
    def __init__(self):
        self.plot_line_table = DataFrame(columns=['configuration', 'term', 'level',
                                                  'J', 'F', 'm_F', 'J_frac', 'F_frac', 'm_F_frac',
                                                  'color', 'y0', 'hf', 'z', 'y', 'x0', 'x1',
                                                  'hflabel', 'hflx', 'hfly', 'zlabel', 'zlx', 'zly',
                                                  'tlx', 'tly'])
        self.plot_transition_table = DataFrame(columns=['F_0', 'J_0', 'configuration_0', 'hf_0', 'level_0', 'm_F_0',
                                                        'term_0', 'x0_0', 'x1_0', 'y_0', 'y0_0', 'z_0', 'F_1', 'J_1',
                                                        'configuration_1', 'hf_1', 'level_1', 'm_F_1', 'term_1', 'x0_1',
                                                        'x1_1', 'y_1', 'y0_1', 'z_1', 'delta_l', 'color', 'wavelength'])

    @staticmethod
    def level_table(level, width=1.0, sublevel_spacing=0.03, scale_splitting=1.0, override_position=False,
                    offset_position=(0.0, 0.0), color='black', zeeman=True, hf=True):
        """
        Generates a level table containing the data needed to draw the diagram

        Args:
            level: the Energy Level being processed
            width: the total width the level should be drawn at
            sublevel_spacing: the spacing between zeeman sublevels
            scale_splitting: the initial scaling of the hyperfine and Zeeman splitting
            override_position: override coordinates
            offset_position: offset to the physical coordinates
            color: the color the level should be drawn in
            zeeman: whether to draw zeeman sublevels
            hf: whether to draw hyperfine sublevels

        Returns: None
        """
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
            table['hflabel'] = 'F=' + str(table['F'])
            table['hflx'] = table['x1'] + 0.05
            table['hfly'] = table['y']
        else:
            delta = sublevel_spacing
            F_max = max(table['F'])
            wd = (width - delta * 2 * F_max) / (2 * F_max + 1)
            table['y'] = table['y0'] + table['hf'] * scale_splitting + table['z'] * scale_splitting
            table['x0'] = x0 + (table['m_F'] + F_max) * (wd + delta)
            table['x1'] = table['x0'] + wd

            table.loc[table['m_F'] == 0, 'hflabel'] = 'F=' + table['F_frac']
            table.loc[table['m_F'] != 0, 'hflabel'] = ''
            table['hflx'] = x0 + width + 0.05
            table['hfly'] = table['y']

            table['zlabel'] = table['m_F_frac']
            table['zlx'] = table['x0']
            table['zly'] = table['y']

        table['tlabel'] = ''
        table.loc[[0], 'tlabel'] = table['term'] + table['J_frac']
        table['tlx'] = x0
        table['tly'] = (max(table['y']) + min(table['y'])) / 2

        return table

    def transition_table(self, transition, x_off_0=0.5, x_off_1=0.5, color=default_lookup):
        """
        Generates the table needed to draw a transition in the diagram

        Args:
            transition: the transition being processed
            x_off_0: the x-offset to apply on the first level
            x_off_1: the x-offset to apply on the second level
            color: the color of the transition

        Returns: None
        """
        table_0 = self.level_table(transition.level_0)
        table_1 = self.level_table(transition.level_1)
        line_0 = table_0.loc[(table_0['m_F'] == transition.m_F_0) & (table_0['F'] == transition.F_0)]
        line_1 = table_1.loc[(table_1['m_F'] == transition.m_F_1) & (table_1['F'] == transition.F_1)]
        line_0 = line_0.drop(['color', 'J_frac', 'F_frac', 'm_F_frac'], axis=1)
        line_1 = line_1.drop(['color', 'J_frac', 'F_frac', 'm_F_frac'], axis=1)
        if line_0.empty or line_1.empty:
            raise ValueError("The selected quantum numbers don't yield a transition. " +
                             "Check that they exist in the specified levels")
        line_0.insert(9, 'x', line_0['x0'] + x_off_0 * (line_0['x1'] - line_0['x0']))
        line_1.insert(9, 'x', line_1['x0'] + x_off_1 * (line_1['x1'] - line_1['x0']))
        line_0.columns = [str(col) + '_0' for col in line_0.columns]
        line_1.columns = [str(col) + '_1' for col in line_1.columns]
        line_0 = line_0.reset_index(drop=True)
        line_1 = line_1.reset_index(drop=True)

        transition_table = pd.concat([line_0, line_1], axis=1, ignore_index=False)

        delta_l = abs(transition_table['level_0'] - transition_table['level_1'])
        wavelength = 299792.458 / delta_l
        if type(color) is dict:
            wl = int(wavelength)
            if wl < min(color.keys()):
                wl = min(color.keys())
            elif wl > max(color.keys()):
                wl = max(color.keys())
            mycolor = color[wl]
        else:
            mycolor = color

        transition_table['delta_l'] = [delta_l]
        transition_table['color'] = [mycolor]
        transition_table['wavelength'] = [wavelength]
        transition_table['name'] = [transition.name]

        return transition_table

    def add_level(self, levels, **kwargs):
        """
        Adds levels to the diagram

        Args:
            levels: the levels to be added
            **kwargs: any settings with which the levels should be drawn

        Returns: None
        """
        for level in levels:
            self.plot_line_table = self.plot_line_table.append(self.level_table(level, **kwargs))
        return self.plot_line_table

    def add_transition(self, transitions, **kwargs):
        """
        Adds transitions to the diagram

        Args:
            transitions: the transitions to be added
            **kwargs: any settings with which the transitions should be drawn

        Returns: None
        """
        for transition in transitions:
            self.plot_transition_table = self.plot_transition_table.append(self.transition_table(transition, **kwargs))
        return self.plot_transition_table

    def remove_level(self, level_names):
        """
        removes the given levels from the table

        Args:
            level_names: the names of the levels to be deleted
        """
        for name in level_names:
            self.plot_line_table = self.plot_line_table[self.plot_line_table.name != name]

    def remove_transition(self, transition_names):
        """
        removes the given transitions from the table

        Args:
            transition_names: the names of the transitions to be deleted
        """
        for name in transition_names:
            self.plot_transition_table = self.plot_transition_table[self.plot_transition_table.name != name]

    def build_figure(self, dimensions=(800, 1000), y_range=(-1e2, 1.3e3), x_range=(-0.5, 4.5), title=None, scale_splitting=1,
                     labels=None, display=False):
        """
        Builds the Bokeh plot

        Args:
            dimensions: the dimensions of the plot
            y_range: the range of frequencies to plot (in THz)
            x_range: the range of J-values to plot
            title: the title of the plot
            scale_splitting: The initial scale of the splittings
            labels: The selection of labels to display. Can be any combination of 'hf', 'term', and 'zeeman', or None
            display: whether to display the table

        Returns:
            a complete Bokeh Grotrian diagram, with magnetic field and scale sliders

        """
        if labels is None:
            labels = []
        p = figure(title=title, plot_width=dimensions[0], plot_height=dimensions[1], y_range=y_range, x_range=x_range)
        line_source = ColumnDataSource(self.plot_line_table)
        arrow_source = ColumnDataSource(self.plot_transition_table)
        print('plotting levels')
        lines = p.segment(x0='x0', y0='y', x1='x1', y1='y',
                          color='color', source=line_source)
        print('plotting transitions')
        arrows = p.segment(x0='x_0', y0='y_0', x1='x_1', y1='y_1',
                           color='color', line_width=3, source=arrow_source)
        # TODO: Maybe make the arrows arrow-y? Might be more trouble than it's worth. Bokeh arrows are insufficient.

        if labels is not []:
            print('drawing labels')
        if 'hf' in labels:
            print('  hf labels')
            hflabels = models.LabelSet(x='hflx', y='y', text='hflabel', level='glyph', source=line_source,
                                       text_baseline='middle', text_font_size='10pt')
            p.add_layout(hflabels)
        if 'zeeman' in labels:
            print('  zeeman labels')
            zlabels = models.LabelSet(x='zlx', y='y', text='zlabel', level='glyph', source=line_source,
                                      text_font_size='8pt')
            p.add_layout(zlabels)
        if 'term' in labels:
            print('  term labels')
            tlabels = models.LabelSet(x='tlx', y='tly', text='tlabel', level='glyph', source=line_source,
                                      text_align='right', text_font_style='bold', text_font_size='12pt',
                                      text_baseline='middle')
            p.add_layout(tlabels)

        print('applying hovertext')
        hover_lines = models.HoverTool(tooltips=[('Term', '@name F=@F_frac, m_F=@m_F_frac'),
                                                 ('Level', '@level{0.000 000 000}')], renderers=[lines])
        hover_arrows = models.HoverTool(tooltips=[('Name', '@name'), ('Frequency', '@delta_l{0.000000} THz'),
                                                  ('Wavelength', '@wavelength{0.00} nm')], renderers=[arrows])
        p.add_tools(hover_lines)
        p.add_tools(hover_arrows)

        print('applying sliders')
        scale_slider = models.Slider(start=1, end=10000, value=scale_splitting, step=10, title='Scaling')
        b_field_slider = models.Slider(start=0, end=20, value=0, step=0.001, title='B-field (G)')
        line_callback = models.CustomJS(args=dict(source=line_source, scale=scale_slider, b_field=b_field_slider),
                                        code='''
                var data = source.data;
                var b_field = b_field.value;
                var scale = scale.value;

                const hf = data['hf'];
                const z = data['z'];
                const y = data['y'];
                const y0 = data['y0'];
                const hfly = data['hfly'];
                const level = data['level'];

                for (var i=0; i < y.length; i++) {
                    y[i] = hf[i]*scale + z[i]*b_field*scale + y0[i];
                    hfly[i] = y[i]
                    level[i] = y0[i] + hf[i] + z[i]*b_field;
                };       
                source.change.emit();
            ''')
        arrow_callback = models.CustomJS(args=dict(source=arrow_source, scale=scale_slider, b_field=b_field_slider),
                                         code='''
                var data = source.data;
                var b_field = b_field.value;
                var scale = scale.value;

                const hf0 = data['hf_0'];
                const hf1 = data['hf_1'];
                const z0 = data['z_0'];
                const z1 = data['z_1'];
                const y0 = data['y_0'];
                const y1 = data['y_1'];
                const y01 = data['y0_1'];
                const y00 = data['y0_0'];
                const level0 = data['level_0'];
                const level1 = data['level_1'];
                const delta_l = data['delta_l'];

                for (var i=0; i < y0.length; i++) {
                    y0[i] = hf0[i]*scale + z0[i]*b_field*scale + y00[i];
                    y1[i] = hf1[i]*scale + z1[i]*b_field*scale + y01[i];
                    level0[i] = y00[i] + hf0[i] + z0[i]*b_field;
                    level1[i] = y01[i] + hf1[i] + z1[i]*b_field;
                    delta_l[i] = level1[i] - level0[i]
                };
                source.change.emit();
            ''')
        scale_slider.js_on_change('value', line_callback, arrow_callback)
        b_field_slider.js_on_change('value', line_callback, arrow_callback)

        if display:
            print('displaying Grotrian diagram')
            show(row(p, column(scale_slider, b_field_slider)))

        return row(p, column(scale_slider, b_field_slider))


class HFPlot:
    def __init__(self, *levels):
        """
        Generates a plot of the hyperfine sublevels and transitions in a given level or levels

        Args:
            *levels: the levels to consider
        """
        self.levels = levels
        self.plot_line_table = DataFrame(columns=[
            'configuration', 'term', 'level',
            'J', 'F', 'm_F', 'J_frac', 'F_frac', 'm_F_frac',
            'color', 'y0', 'hf', 'z', 'y', 'x0', 'x1'])
        self.plot_arrow_table = DataFrame(columns=[
            'F_0', 'hf_0', 'level_0', 'm_F_0', 'term_0', 'x0_0', 'x1_0', 'y_0', 'y0_0', 'z_0',
            'F_1', 'hf_1', 'level_1', 'm_F_1', 'term_1', 'x0_1', 'x1_1', 'y_1', 'y0_1', 'z_1',
            'delta_l', 'wavelength'])
        if self.levels is None:
            self.levels = []
        if len(self.levels) == 1:
            self.internal = True
        else:
            self.internal = False

        for level in self.levels:
            self.plot_line_table = self.plot_line_table.append(self.level_table(level[0], level[1]), ignore_index=True)

        self.plot_arrow_table = self.arrow_table(self.plot_line_table, spacing='physical')

    @staticmethod
    def level_table(level, F, width=0.01, b_field=1, scale_splitting_hf=1.0, scale_splitting_z=1.0, color='black', y0=0):
        """

        Args:
            level: the level to be processed
            F: the F-value of the level
            width: the total width of the level
            b_field: the magnetic field
            scale_splitting_hf: the hyperfine scaling
            scale_splitting_z: the extra scaling to be applied for the Zeeman effect
            color: the color of the level
            y0: the y-offset at which to draw the level

        Returns: None
        """
        table = level.data_table(hf=True, zeeman=True, b_field=b_field)
        table = table[table['F'] == F].reset_index()
        table['color'] = color
        table['name'] = level.name

        table['y0'] = y0

        table['y'] = table['y0'] + table['hf'] * scale_splitting_hf + table['z'] * scale_splitting_z
        table['level'] = table['y0'] + table['hf'] + table['z'] * b_field

        table['x0'] = max(table['level']) - min(table['level'])
        table['x1'] = min(table['level']) - max(table['level'])

        table.loc[table['m_F'] == 0, 'hflabel'] = 'F=' + table['F_frac']
        table.loc[table['m_F'] != 0, 'hflabel'] = ''
        table['hflx'] = width + 0.05
        table['hfly'] = table['y']

        table['zlabel'] = table['m_F_frac']
        table['zlx'] = table['x0']
        table['zly'] = table['y']

        table['tlabel'] = ''
        table.loc[[0], 'tlabel'] = table['term'] + table['J_frac']
        table['tlx'] = table['x0']
        table['tly'] = (max(table['y']) + min(table['y'])) / 2

        return table

    def arrow_table(self, level_table, spacing='physical'):
        """
        Takes in the plot's level table and popultes a new table containing all the allowed transitions

        Args:
            level_table: the level table being processed
            spacing: can be 'physical' or 'spaced'

        Returns:
            a table containing all the data for the transitions

        """
        from transition_strengths import M1_transition_strength_avg, E1_transition_strength_avg, E2_transition_strength_avg

        table = DataFrame(columns=['F_0', 'hf_0', 'm_F_0', 'y_0', 'y0_0', 'z_0',
                                   'F_1', 'hf_1', 'm_F_1', 'y_1', 'y0_1', 'z_1',
                                   'level_0', 'level_1', 'delta_l', 'x', 'J0', 'J1',
                                   'strength', 'color'])
        for index0, sl0 in level_table.iterrows():
            for index1, sl1 in level_table.iterrows():
                # TODO: Make this iteration more efficient
                I = sl0['I']
                L_0, S_0, J_0, F_0, m_F_0 = sl0['L'], sl0['S'], sl0['J'], sl0['F'], sl0['m_F']
                L_1, S_1, J_1, F_1, m_F_1 = sl1['L'], sl1['S'], sl1['J'], sl1['F'], sl1['m_F']
                term_0, term_1 = sl0['term'], sl1['term']

                if term_0 == term_1:
                    strength = M1_transition_strength_avg(I, L_0, S_0, J_0, F_0, m_F_0, L_1, S_1, J_1, F_1, m_F_1)
                elif abs(J_0 - J_1) > 1:
                    strength = E2_transition_strength_avg(I, L_0, S_0, J_0, F_0, m_F_0, L_1, S_1, J_1, F_1, m_F_1)
                else:
                    strength = E1_transition_strength_avg(I, L_0, S_0, J_0, F_0, m_F_0, L_1, S_1, J_1, F_1, m_F_1, depth=1)

                if m_F_0 == m_F_1:
                    color = 'black'
                elif m_F_0 - m_F_1 == 1:
                    color = 'orange'
                elif m_F_1 - m_F_0 == 1:
                    color = 'blue'
                else:
                    color = 'red'

                # if E1_str == 0 and F_0 != F_1:
                #     print 'strength for F={} m={} to F={} m={} was 0'.format(F_0, m_F_0, F_1, m_F_1)
                if index0 <= index1 and strength != 0:
                    if (term_0 != term_1) or (self.internal) or (not self.internal and F_0 != F_1):
                        line = DataFrame(data={'F_0': [F_0], 'hf_0': [sl0['hf']], 'm_F_0': [m_F_0],
                                               'y_0': [sl0['y']], 'y0_0': [sl0['y0']], 'z_0': [sl0['z']],
                                               'F_1': [F_1], 'hf_1': [sl1['hf']], 'm_F_1': [m_F_1],
                                               'y_1': [sl1['y']], 'y0_1': [sl1['y0']], 'z_1': [sl1['z']],
                                               'level_0': sl0['level'], 'level_1': sl1['level'],
                                               'delta_l': abs(sl0['level'] - sl1['level']),
                                               'x': abs(sl0['level'] - sl1['level']),
                                               'J_0': [J_0], 'J_1': [J_1],
                                               'strength': [strength], 'color': [color]})
                        table = table.append(line, ignore_index=True)

        def space_out_lines(series, new_spacing):
            # TODO: Fix this
            s = copy.deepcopy(series)
            s.sort_values(ascending=True, inplace=True)
            for i in range(len(s) - 1):
                if s[i + 1] - s[i] < new_spacing:
                    s[i] -= new_spacing - (s[i + 1] - s[i]) / 2
                    s[i + 1] += new_spacing - (s[i + 1] - s[i]) / 2
                    s = space_out_lines(s, new_spacing)
            return s.sort_index()

        if spacing != 'physical':
            table['x'] = space_out_lines(table['x'], spacing)

        return table

    def add_level(self, *levels):
        self.levels = self.levels.append(levels)
        if len(self.levels) == 1:
            self.internal = True
        else:
            self.internal = False
        for level in levels:
            self.plot_line_table = self.plot_line_table.append(self.level_table(level[0], level[1]))
        self.plot_arrow_table = self.arrow_table(self.plot_line_table)
        return self.plot_line_table, self.plot_arrow_table

    def remove_level(self, *level_names):
        for name in level_names:
            if type(name) is tuple:
                self.plot_line_table = self.plot_line_table[
                    self.plot_line_table.F != name[1] and self.plot_line_table.name != name[0]]
                self.plot_arrow_table = self.plot_arrow_table[
                    self.plot_arrow_table.F_0 != name[1] and self.plot_arrow_table.name_0 != name[0]]
                self.plot_arrow_table = self.plot_arrow_table[
                    self.plot_arrow_table.F_1 != name[1] and self.plot_arrow_table.name_1 != name[0]]
            else:
                self.plot_line_table = self.plot_line_table[self.plot_line_table.name != name]
                self.plot_arrow_table = self.plot_arrow_table[self.plot_arrow_table.name_0 != name]
                self.plot_arrow_table = self.plot_arrow_table[self.plot_arrow_table.name_0 != name]

    def build_figure(self, dimensions=(1600, 1000), title=None, display=False, labels=None, sliders=None, x_range=None):
        """

        Args:
            dimensions: the dimensions of the figure
            title: the title
            display: whether to display the final figure
            labels: which labels to include. Can be any of 'zeeman', 'level', 'term'
            sliders: which sliders to include. HF plot listens for 'b_field_slider', 'scale_hf_slider', and 'scale_z_slider'
            x_range: the range of the plot

        Returns:
            a bokeh figure, and a column of the relevant sliders

        """
        import sliders as sli

        if labels is None:
            labels = []
        if sliders is None:
            sliders = {}
        if 'b_field_slider' not in sliders.keys():
            sliders['b_field_slider'] = sli.b_field_slider
        if 'scale_hf_slider' not in sliders.keys():
            sliders['scale_hf_slider'] = sli.scale_hf_slider
        if 'scale_z_slider' not in sliders.keys():
            sliders['scale_z_slider'] = sli.scale_z_slider

        scale_hf_slider = sliders['scale_hf_slider']
        scale_z_slider = sliders['scale_z_slider']
        b_field_slider = sliders['b_field_slider']

        p = figure(title=title, plot_width=dimensions[0], plot_height=dimensions[1],
                   y_range=(min(self.plot_line_table['y']), max(self.plot_line_table['y'])),
                   x_range=(min(self.plot_arrow_table['delta_l']) - 1e-6, max(self.plot_arrow_table['delta_l']) + 1e-6))
        line_source = ColumnDataSource(self.plot_line_table)
        arrow_source = ColumnDataSource(self.plot_arrow_table)
        print('plotting levels')
        if x_range is not None:
            x_min = min(self.plot_arrow_table['delta_l'])
            x_max = max(self.plot_arrow_table['delta_l'])
            lines = p.segment(x0=x_min, y0='y', x1=x_max, y1='y',
                              color='color', source=line_source)
        else:
            lines = p.segment(x0='x0', y0='y', x1='x1', y1='y',
                              color='color', source=line_source)
        print('plotting transitions')
        arrows = p.segment(x0='delta_l', y0='y_0', x1='delta_l', y1='y_1',
                           line_width=3, color='color', source=arrow_source)

        if labels is not []:
            print('drawing labels')
        if 'hf' in labels:
            print('  hf labels')
            hflabels = models.LabelSet(x='hflx', y='y', text='hflabel', level='glyph', source=line_source,
                                       text_baseline='middle', text_font_size='10pt')
            p.add_layout(hflabels)
        if 'zeeman' in labels:
            print('  zeeman labels')
            zlabels = models.LabelSet(x='zlx', y='y', text='zlabel', level='glyph', source=line_source,
                                      text_font_size='8pt')
            p.add_layout(zlabels)
        if 'term' in labels:
            print('  term labels')
            tlabels = models.LabelSet(x='tlx', y='tly', text='tlabel', level='glyph', source=line_source,
                                      text_align='right', text_font_style='bold', text_font_size='12pt',
                                      text_baseline='middle')
            p.add_layout(tlabels)

        print('applying hovertext')
        hover_lines = models.HoverTool(tooltips=[('Term', '@name F=@F_frac, m_F=@m_F_frac'),
                                                 ('Level', '@level{0.000 000 000}')], renderers=[lines])
        hover_arrows = models.HoverTool(tooltips=[('Name', '@name'), ('Frequency', '@delta_l{0.000000} THz')], renderers=[arrows])
        p.add_tools(hover_lines)
        p.add_tools(hover_arrows)

        print('applying sliders')

        line_callback = models.CustomJS(
            args=dict(source=line_source, hf_scale=scale_hf_slider, z_scale=scale_z_slider, b_field=b_field_slider),
            code='''
                   var data = source.data;
                   var b_field = b_field.value;
                   var hf_scale = hf_scale.value;
                   var z_scale = z_scale.value

                   const hf = data['hf'];
                   const z = data['z'];
                   const y = data['y'];
                   const y0 = data['y0'];
                   const hfly = data['hfly'];
                   const level = data['level'];

                   for (var i=0; i < y.length; i++) {
                       y[i] = hf[i]*hf_scale + z[i]*b_field*z_scale + y0[i];
                       hfly[i] = y[i];
                       level[i] = y0[i] + hf[i] + z[i]*b_field;
                   };       
                   source.change.emit();
                   ''')
        arrow_callback = models.CustomJS(
            args=dict(source=arrow_source, hf_scale=scale_hf_slider, z_scale=scale_z_slider, b_field=b_field_slider),
            code='''
                   var data = source.data;
                   var b_field = b_field.value;
                   var hf_scale = hf_scale.value;
                   var z_scale = z_scale.value;

                   const hf0 = data['hf_0'];
                   const hf1 = data['hf_1'];
                   const z0 = data['z_0'];
                   const z1 = data['z_1'];
                   const y0 = data['y_0'];
                   const y1 = data['y_1'];
                   const y01 = data['y0_1'];
                   const y00 = data['y0_0'];
                   const level0 = data['level_0'];
                   const level1 = data['level_1'];
                   const delta_l = data['delta_l'];

                   for (var i=0; i < y0.length; i++) {
                       y0[i] = hf0[i]*hf_scale + z0[i]*b_field*z_scale + y00[i];
                       y1[i] = hf1[i]*hf_scale + z1[i]*b_field*z_scale + y01[i];
                       level0[i] = hf0[i] + z0[i]*b_field + y00[i];
                       level1[i] = hf1[i] + z1[i]*b_field + y01[i];
                       delta_l[i] = Math.abs(level1[i] - level0[i]);
                   };
                   source.change.emit();
                   ''')
        scale_hf_slider.js_on_change('value', line_callback, arrow_callback)
        scale_z_slider.js_on_change('value', line_callback, arrow_callback)
        b_field_slider.js_on_change('value', line_callback, arrow_callback)

        if display:
            print('displaying Hyperfine diagram')
            show(row(p, column(scale_hf_slider, scale_z_slider, b_field_slider)))

        return p, column(scale_hf_slider, scale_z_slider, b_field_slider)


class LorentzianPlot(HFPlot):

    # TODO: Autoscale point density by minimum linewidth
    # TODO: Maintain a vertical axis, at least optionally

    def build_figure(self, linewidth=5e-9, dimensions=(1600, 400), title=None, display=False, x_range=None, labels=None,
                     sliders=None):
        """

       Args:
           linewidth: the linewidth of the transition (in THz)
           dimensions: the dimensions of the figure
           title: the title
           display: whether to display the final figure
           x_range: the range of frequencies the plot should care about
           labels: which labels to include.
           sliders: which sliders to include. Lorentz plot listens for 'b_field_slider', 'linewidth_slider'

       Returns:
           a bokeh figure, and a column of the relevant sliders

        """
        import numpy as np
        from math import pi
        import sliders as sli

        # if labels is None:
        #     labels = []
        # TODO: implement labels
        if sliders is None:
            sliders = {}
        if 'linewidth_slider' not in sliders.keys():
            sliders['linewidth_slider'] = sli.linewidth_slider
        if 'b_field_slider' not in sliders.keys():
            sliders['b_field_slider'] = sli.b_field_slider

        b_field_slider = sliders['b_field_slider']
        linewidth_slider = sliders['linewidth_slider']

        if x_range is None:
            x_range = (min(self.plot_arrow_table['delta_l']) - 1e-6, max(self.plot_arrow_table['delta_l']) + 1e-6)
        p = figure(title=title, plot_width=dimensions[0], plot_height=dimensions[1], x_range=x_range)

        x_axis = np.linspace(min(self.plot_arrow_table['delta_l']) - 1e-6, max(self.plot_arrow_table['delta_l']) + 1e-6, 30000)

        lines = []
        for index, transition in self.plot_arrow_table.iterrows():
            line = (1 / (2 * pi)) * linewidth / ((x_axis - transition['delta_l']) ** 2 + linewidth ** 2 / 4)
            strength = transition['strength']
            lines.append(line * strength)

        total_line = np.sum(lines)

        transition_data = ColumnDataSource(data={
            'hf_0': self.plot_arrow_table['hf_0'],
            'hf_1': self.plot_arrow_table['hf_1'],
            'z_0': self.plot_arrow_table['z_0'],
            'z_1': self.plot_arrow_table['z_1'],
            'y0_1': self.plot_arrow_table['y0_1'],
            'y0_0': self.plot_arrow_table['y0_0'],
            'level_0': self.plot_arrow_table['level_0'],
            'level_1': self.plot_arrow_table['level_1'],
            'delta_l': self.plot_arrow_table['delta_l'],
            'strength': self.plot_arrow_table['strength']
        })

        lorentz_table = DataFrame(data={'x_axis': x_axis, 'total_line': total_line})  # , 'transition_data': transition_data})
        lorentz_source = ColumnDataSource(lorentz_table)

        print('plotting transitions')

        p.line(x='x_axis', y='total_line', source=lorentz_source)

        print('applying sliders')

        lorentz_callback = models.CustomJS(
            args=dict(source=lorentz_source, b_field=b_field_slider, linewidth=linewidth_slider, transition_data=transition_data),
            code='''
                var line_data = source.data;
                var transition_data = transition_data.data;
                var b_field = b_field.value;
                var linewidth = linewidth.value;
                
                const hf0 = transition_data['hf_0'];
                const hf1 = transition_data['hf_1'];
                const z0 = transition_data['z_0'];
                const z1 = transition_data['z_1'];
                const y01 = transition_data['y0_1'];
                const y00 = transition_data['y0_0'];
                const level0 = transition_data['level_0'];
                const level1 = transition_data['level_1'];
                const delta_l = transition_data['delta_l'];
                const strength = transition_data['strength'];
                
                const x_axis = line_data['x_axis'];
                const total_line = line_data['total_line'];
                                
                for (i=0; i < hf0.length; i++) {
                    level0[i] = hf0[i] + z0[i]*b_field + y00[i];
                    level1[i] = hf1[i] + z1[i]*b_field + y01[i];
                    delta_l[i] = Math.abs(level1[i] - level0[i]);
                };
                
                for (var i=0; i< x_axis.length; i++){
                    let tot = 0;
                    for (var j=0; j<hf0.length; j++){
                        tot += strength[j]*(0.5/Math.PI)*linewidth/((x_axis[i]-delta_l[j])*
                        (x_axis[i]-delta_l[j])+linewidth*linewidth/4);
                    };
                    total_line[i] = tot;
                };
                
                source.change.emit();
                ''')

        b_field_slider.js_on_change('value', lorentz_callback)
        linewidth_slider.js_on_change('value', lorentz_callback)

        if display:
            print('displaying Lorentzian diagram')
            show(row(p, column(b_field_slider, linewidth_slider)))

        return p, column(b_field_slider, linewidth_slider)


if __name__ == '__main__':
    import pandas as pd
    from atom_import import load_from_pickle
    from colors import uv_ir_lookup
    from sliders import default_sliders

    Yb_171 = load_from_pickle("171Yb.atom")
    Yb_173 = load_from_pickle("173Yb.atom")
    Yb_174 = load_from_pickle("174Yb.atom")

    pd.set_option('display.max_rows', 500)
    pd.set_option('display.max_columns', 500)
    pd.set_option('display.width', 1000)


    def MakeGrotrian(atom):
        g = Grotrian()
        levels = atom.levels.values()
        g.add_level(levels, color='black')
        g.add_transition(atom.transitions.values(), color=uv_ir_lookup)

        g.build_figure(display=True)  # , labels=['hf', 'zeeman', 'term'])


    def MakeMixedPlot(l0, l1):
        l = LorentzianPlot(l0, l1)
        L_plot = l.build_figure(dimensions=(1600, 400), sliders=default_sliders)

        h = HFPlot(l0, l1)
        HF_plot = h.build_figure(dimensions=(1600, 400), x_range=L_plot[0].x_range, sliders=default_sliders)

        show(row(gridplot([[L_plot[0]], [HF_plot[0]]]), column(list(default_sliders.values()))))


    def MakeLorentzPlot(l0, l1):
        l = LorentzianPlot(l0, l1)
        l.build_figure(display=True, linewidth=1e-7)


    def MakeHFPlot(l0, l1):
        h = HFPlot(l0, l1)
        h.build_figure(display=True)


    # level0 = (Yb_171.levels['2F*7/2'], 3)
    # level1 = (Yb_171.levels['2F*7/2'], 4)

    # level0 = (Yb_171.levels['2S1/2'], 1)
    # level1 = (Yb_171.levels['2S1/2'], 0)

    # level0 = (Yb_171.levels['2D5/2'], 3)
    # level1 = (Yb_171.levels['2P*3/2'], 2)

    level0 = (Yb_171.levels['2D3/2'], 2)
    level1 = (Yb_171.levels['1[3/2]*3/2'], 2)

    # MakeGrotrian(Yb_171)
    MakeMixedPlot(level0, level1)
    # MakeLorentzPlot(level0, level1)
    # MakeHFPlot(level0, level1)
