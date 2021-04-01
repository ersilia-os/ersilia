import seaborn as sns
import matplotlib


def set_style(style=None):
    """Set basic plotting style and fonts."""
    try:
        matplotlib.font_manager._rebuild()
    except Exception as ex:
        print(str(ex))
    if style is None:
        style = ('ticks', {
            'font.family': 'sans-serif',
            'font.serif': ['Arial'],
            'font.size': 16,
            'axes.grid': True})
    else:
        style = style
    sns.set_style(*style)
