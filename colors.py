from scipy import interpolate

def wavelength_to_color(wl, colormap):
    """
    Takes a frequency (in Hertz for now) and outputs a hex color value based on the frequency
    It does this by interpolating between set points where I know that the colors should be.
    """

    wavelengths, rs, gs, bs = colormap

    rf = interpolate.interp1d(wavelengths, rs)
    gf = interpolate.interp1d(wavelengths, gs)
    bf = interpolate.interp1d(wavelengths, bs)

    if wl >= max(wavelengths):
        wl = max(wavelengths)
    elif wl <= min(wavelengths):
        wl = min(wavelengths)
    r, g, b = rf(wl), gf(wl), bf(wl)

    return "#{:02x}{:02x}{:02x}".format(int(r), int(g), int(b))

def make_lookup(colormap):
    lookup = {}
    for wl in range(colormap[0][0], colormap[0][-1]+1):
        lookup[wl] = wavelength_to_color(wl, colormap)
    return lookup

# colormap format: wavelengths, rs, gs, bs

default_map = ([280, 300, 350, 400, 445, 475, 493, 510, 570, 650, 780, 1000, 1500],
               [0, 50, 120, 80, 110, 0, 40, 40, 255, 235, 90, 50, 0],
               [0, 20, 120, 40, 0, 0, 255, 255, 230, 30, 30, 20, 0],
               [0, 100, 255, 120, 255, 255, 250, 0, 0, 30, 30, 20, 0])

rich_ir_map = ([280, 300, 350, 400, 445, 475, 493, 510, 570, 650, 780, 1000, 1500, 3000],
               [0, 50, 120, 80, 110, 0, 40, 40, 255, 235, 90, 80, 50, 0],
               [0, 20, 120, 40, 0, 0, 255, 255, 230, 30, 30, 20, 0, 0],
               [0, 100, 255, 120, 255, 255, 250, 0, 0, 30, 30, 20, 0, 0])

uv_ir_map = ([280, 300, 350, 400, 445, 475, 493, 510, 570, 650, 780, 1000, 1500, 3000],
                [0,   50,  120, 80,  110, 0,   40,  40,  255, 235, 90,  80,   50,   0],
                [255, 255, 120, 40,  0,   0,   255, 255, 230, 30,  30,  20,   0,    0],
                [0,   100, 255, 120, 255, 255, 250, 0,   0,   30,  30,  20,   0,    0])

green_uv_map = ([280, 300, 350, 400, 445, 475, 493, 510, 570, 650, 780, 1000, 1500],
                [0, 50, 120, 80, 110, 0, 40, 40, 255, 235, 90, 50, 0],
                [255, 255, 120, 40, 0, 0, 255, 255, 230, 30, 30, 20, 0],
                [0, 100, 255, 120, 255, 255, 250, 0, 0, 30, 30, 20, 0])

#  lookup tables are what should be imported

default_lookup = make_lookup(default_map)
rich_ir_lookup = make_lookup(rich_ir_map)
uv_ir_lookup = make_lookup(uv_ir_map)
green_uv_lookup = make_lookup(green_uv_map)

if __name__ == "__main__":
    pass
