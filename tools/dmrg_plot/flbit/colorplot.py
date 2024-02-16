import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

cmaps = [('Perceptually Uniform Sequential', [
            'viridis', 'plasma', 'inferno', 'magma', 'cividis']),
         ]

gradient = np.linspace(0, 1, 256)
gradient = np.vstack((gradient, gradient))


def plot_color_gradients(cmap_category, cmap_list):
    # Create figure and adjust figure height to number of colormaps
    nrows = len(cmap_list)
    figh = 0.35 + 0.15 + (nrows + (nrows-1)*0.1)*0.22
    fig, axs = plt.subplots(nrows=nrows, figsize=(6.4, figh))
    fig.subplots_adjust(top=1-.35/figh, bottom=.15/figh, left=0.2, right=0.99)

    axs[0].set_title(f"{cmap_category} colormaps", fontsize=14)

    palette = sns.color_palette("viridis", n_colors=4, as_cmap=False)
    for ax, cmap_name in zip(axs, cmap_list):
        cmap = plt.colormaps[cmap_name]
        palette = [cmap(i) for i in range(cmap.N)]
        print(palette)
        for idx, color in enumerate(palette):
            ax.vlines(x=idx, ymin=0, ymax=1, color=color)
        # ax.imshow(gradient, aspect='auto', cmap=palette)
        # ax.text(-.01, .5, cmap_name, va='center', ha='right', fontsize=10,
        #         transform=ax.transAxes)

    # Turn off *all* ticks & spines, not just the ones with colormaps.
    for ax in axs:
        ax.set_axis_off()


for cmap_category, cmap_list in cmaps:
    plot_color_gradients(cmap_category, cmap_list)

plt.show()