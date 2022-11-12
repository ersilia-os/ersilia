---
description: Sytlia is a small Python library for styling plots
---

# Scientific figures with Stylia

[Stylia](https://github.com/ersilia-os/stylia) is a small package to stylize [Matplotlib](https://matplotlib.org/) plots in Python so that they are publication-ready. Stylia provides modified axes (`ax`) that can be used as drop-in replacements for Matplotlib axes.

## Getting started

### Installation

Stylia is constantly evolving, so we recommend that you install it directly from the [GitHub repository](https://github.com/ersilia-os/stylia).

```python
git clone https://github.com/ersilia-os/stylia.git
cd stylia
pip install -e . 
```

### Create a single panel figure

```python
import stylia as st
import numpy as np

fig, axs = st.create_figure(1,1)
ax = axs[0]

x = np.random.normal(size=100)
y = np.random.normal(size=100)

ax.scatter(x, y)

st.save_figure("my_single_plot.png")
```

### Create a multipanel figure

```python
import stylia as st
import numpy as np

# create a figure to be used in a slide
fig, axs = st.create_figure(nrows=2, ncols=2, width_ratios=[2, 1])

# get data
x = np.random.normal(size=100)
y = np.random.normal(size=100)

# first plot, access with flat subplots coordinates (0)
ax = axs[0]
# a default color is used
ax.scatter(x, y)
# write labels to axis, title and numbering of the subplot
st.label(
    ax,
    title="My first plot",
    xlabel="This is the X axis",
    ylabel="This is the Y axis",
    abc="A",
)

# second plot, acces with subplots 2D coordinates (0,1)
ax = axs[0, 1]
# use a named color
named_colors = st.NamedColors()

def my_scatterplot(ax, x, y):
    ax.scatter(x, y, color=named_colors.get("red"))

my_scatterplot(ax, x, y)
# write only a new title (the rest are defaults)
st.label(ax, title="My second plot")

# third plot, access with next() method
ax = axs.next()
cmap = st.ContinuousColorMap()
cmap.fit(x)
colors = cmap.get(x)
ax.scatter(x, y, color=colors)

# fourth plot
ax = axs[1, 1]
# add transparency
ax.scatter(x, y, color=named_colors.get("blue", alpha=0.2))

# save figure
st.save_figure("my_grid_plot.png")
```

## Sizes

### Figure size

We follow the [Nature Figure Guidelines](https://www.nature.com/nature/for-authors/formatting-guide). Please read those style guidelines carefully. In brief, the entire figure should be have the following sizes:

* `SINGLE_COLUMN_WIDTH`: 90 mm or 3.54 in
* `TWO_COLOUMNS_WIDTH`: 180 mm or 7.09 in

These variables are built-in within Stylia. You can access them as follows:

```python
from stylia import TWO_COLUMNS_WIDTH
```

### Font size

* `FONTSIZE_SMALL`: 5
* `FONTSIZE`: 6
* `FONTSIZE_BIG` : 8

### Marker sizes

* `MARKERSIZE_SMALL`: 5
* `MARKERSIZE`: 10
* `MARKERSIZE_BIG`: 30

### Line widths

* `LINEWIDTH`: 0.5
* `LINEWIDTH_THICK`: 1.0

## Colors

### Named colors

You can use predefined colors:

```python
from stylia import NamedColors

named_colors = NamedColors()
color = named_colors.get('blue')
```

Available color names are:

* `'red'`
* `'blue'`
* `'green'`
* `'orange'`
* `'purple'`
* `'yellow'`
* `'gray'`
* `'white'`
* `'black'`

### Color maps

#### Continuous color maps

Color maps can be created with the `fit` method.

```python
from stylia import ContinuousColorMap
import numpy as np

cmap = ContinuousColorMap("spectral")
x = np.random.normal(size=100)
y = np.random.normal(size=200) / 2
cmap.fit(x)
# get colors of x
colors_x = cmap.get(x)
# get colors of y based on the x scale
colors_y = cmap.get(y)
```

Available color maps are:

* `'spectral'`
* `'viridis'`
* `'coolwarm'`

#### Discrete colormaps

{% hint style="warning" %}
Discrete colormaps are work in progress
{% endhint %}

{% hint style="info" %}
Please note that, by default, we use [Scientific Color Maps](https://www.fabiocrameri.ch/colourmaps/). Non-scientific color maps look brighter, though. If you want to use non-scientific color maps, simply specify `ContinuousColorMap("spectral", scientific=False)`.
{% endhint %}
