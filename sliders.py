from bokeh.models import Slider

b_field_slider = Slider(start=0, end=5, value=1, step=0.0001, title='B-field (G)')
linewidth_slider = Slider(start=1e-9, end=1e-7, value=5e-9, step=1e-9, title='linewidth (THz)')
scale_hf_slider = Slider(start=0.1, end=2, value=1, step=0.01, title='HF Scaling')
scale_z_slider = Slider(start=1, end=1000, value=1, step=0.5, title='Zeeman Scaling')

default_sliders = {'b_field_slider': b_field_slider,
                   'linewidth_slider': linewidth_slider,
                   'scale_hf_slider': scale_hf_slider,
                   'scale_z_slider': scale_z_slider}
