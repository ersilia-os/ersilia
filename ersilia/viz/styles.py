import seaborn as sns
import matplotlib


def set_style(style=None):
    """Set basic plotting style and fonts."""
    matplotlib.font_manager._rebuild()
    if style is None:
        style = (
            "ticks",
            {
                "font.family": "sans-serif",
                "font.serif": ["Arial"],
                "font.size": 16,
                "axes.grid": True,
            },
        )
    else:
        style = style
    sns.set_style(*style)
