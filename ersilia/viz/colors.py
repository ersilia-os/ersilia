import numpy as np
import colorsys
import matplotlib.colors as mc

def rgb2hex(r, g, b):
    """RGB to hexadecimal."""
    return '#%02x%02x%02x' % (r, g, b)

# Ersilia colors
eos = {
	"dark": "#50285a",
	"gray": "#d2d2d0",
	"green": "#bee6b4",
	"white": "#ffffff",
	"purple": "#aa96fa",
	"pink": "#dca0dc",
	"yellow": "#fad782",
	"blue": "#8cc8fa",
	"red": "#faa08c"
}

# Chemical Checker colors
cc = {
    "red": rgb2hex(250, 100, 80),
	"purple": rgb2hex(200, 100, 225),
    "blue": rgb2hex(80, 120, 220),
	"green": rgb2hex(120, 180, 60),
	"orange": rgb2hex(250, 100, 80)
}


def lighten_color(color, amount=0.5):
    """Lighthen a color."""
    try:
        c = mc.cnames[color]
    except Exception:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])


def cc_colors(coord, lighness=0):
    """Predefined CC colors."""
    colors = {
        'A': ['#EA5A49', '#EE7B6D', '#F7BDB6'],
        'B': ['#B16BA8', '#C189B9', '#D0A6CB'],
        'C': ['#5A72B5', '#7B8EC4', '#9CAAD3'],
        'D': ['#7CAF2A', '#96BF55', '#B0CF7F'],
        'E': ['#F39426', '#F5A951', '#F8BF7D'],
        'Z': ['#000000', '#666666', '#999999']}
    return colors[coord[:1]][lighness]


def make_cmap(colors, position=None, bit=False):
    bit_rgb = np.linspace(0, 1, 256)
    if position is None:
        position = np.linspace(0, 1, len(colors))
    else:
        if len(position) != len(colors):
            raise Exception("position length must be the same as colors")
        elif position[0] != 0 or position[-1] != 1:
            raise Exception("position must start with 0 and end with 1")
    if bit:
        for i in range(len(colors)):
            colors[i] = (bit_rgb[colors[i][0]],
                         bit_rgb[colors[i][1]],
                         bit_rgb[colors[i][2]])
    cdict = {'red': [], 'green': [], 'blue': []}
    for pos, color in zip(position, colors):
        cdict['red'].append((pos, color[0], color[0]))
        cdict['green'].append((pos, color[1], color[1]))
        cdict['blue'].append((pos, color[2], color[2]))

    cmap = matplotlib.colors.LinearSegmentedColormap(
        'my_colormap', cdict, 256)
    return cmap
