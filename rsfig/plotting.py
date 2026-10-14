"""Panel drawing and page layout of the summary figure."""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

from . import COLORS

MUTED = "#8a8984"


def style(p):
    f = p["font"]
    return {
        "font.size": f, "axes.labelsize": f, "axes.titlesize": f,
        "xtick.labelsize": f, "ytick.labelsize": f, "legend.fontsize": f,
        "font.family": "serif", "mathtext.fontset": "dejavusans",
        "axes.edgecolor": "black", "axes.labelcolor": "black", "axes.linewidth": 0.8,
        "xtick.color": "black", "ytick.color": "black",
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.size": 4, "ytick.major.size": 4,
        "xtick.major.width": 0.8, "ytick.major.width": 0.8,
        "xtick.minor.visible": False, "ytick.minor.visible": False,
        "axes.spines.top": True, "axes.spines.right": True,
        "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
        "legend.frameon": True, "legend.fancybox": True, "legend.framealpha": 1,
        "legend.edgecolor": "#dcdcdc", "legend.borderpad": 0.4, "legend.handletextpad": 0.3,
        "legend.labelspacing": 0.3, "legend.borderaxespad": 0.4,
    }


def spread_positions(x, min_sep):
    """
    Nudge x positions apart so neighbours are at least min_sep apart, keeping
    each group of crowded points centred on its original mean (display only).
    """
    order = np.argsort(x)
    xs = x[order]
    groups = [[i] for i in range(len(xs))]

    def placed(g):
        return xs[g].mean() + min_sep * (np.arange(len(g)) - (len(g) - 1) / 2)

    merged = True
    while merged:
        merged = False
        for k in range(len(groups) - 1):
            if placed(groups[k + 1])[0] - placed(groups[k])[-1] < min_sep:
                groups[k] += groups.pop(k + 1)
                merged = True
                break

    shifted = np.empty_like(x)
    shifted[order] = np.concatenate([placed(g) for g in groups])
    return shifted


def _arr(col):
    return np.ma.filled(np.ma.asarray(col, dtype=float), np.nan)


def _legend(ax, **kw):
    leg = ax.legend(**kw)
    leg.get_frame().set_linewidth(0.6)
    return leg


def _panel_label(ax, letter, font, x=0.02, y=0.97, **kw):
    ax.text(x, y, f"({letter})", transform=ax.transAxes, ha="left", va="top",
            fontsize=font, fontweight="bold", zorder=10, **kw)


def draw_cmd(ax, an, p):
    ax.scatter(_arr(an.m_b), _arr(an.color_index), s=2, alpha=0.08, color="#393939", lw=0,
               label="All galaxies", rasterized=True)
    mm = np.linspace(an.m_lower_limit, an.m_upper_limit, 100)
    line = an.a * mm + an.b
    ax.fill_between(mm, line - an.tolerance, line + an.tolerance, color="firebrick",
                    alpha=0.18, edgecolor="firebrick", linewidth=0.6,
                    label="Red-sequence selection")
    ax.scatter(_arr(an.m_b_mid), _arr(an.color_index_mid), s=5, edgecolors="k",
               linewidths=0.3, facecolors="#13FF02", label="Known members")
    ax.scatter(np.concatenate([_arr(an.m_b_low), _arr(an.m_b_high)]),
               np.concatenate([_arr(an.color_index_low), _arr(an.color_index_high)]),
               s=5, edgecolors="k", linewidths=0.3, facecolors="#FF1A7D",
               label="Known non-members")
    ax.set_xlabel(r"$b$ [mag]")
    ax.set_ylabel(r"$b-g$ [mag]")
    ax.set_xlim(17.0, 25.5)
    ax.set_ylim(-0.62, 1.5)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=6))
    leg = _legend(ax, loc="lower left", ncol=2, markerscale=2.0, columnspacing=0.8)
    leg.legend_handles[0].set_alpha(0.6)
    _panel_label(ax, "a", p["font"])


def draw_density(fig, ax, an, p):
    im = ax.imshow(an.smoothed_density.T[:, ::-1], origin="lower", aspect="equal",
                   cmap="turbo", interpolation="bicubic",
                   extent=[an.ra_max, an.ra_min, an.dec_min, an.dec_max])
    ax.xaxis.set_major_locator(MaxNLocator(nbins=3))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
    ax.set_xlabel("RA [deg]")
    ax.set_ylabel("Dec [deg]")
    cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.03)
    label = ("Luminosity density [counts arcmin$^{-2}$]" if an.density_mode == "luminosity"
             else "Number density [galaxies arcmin$^{-2}$]")
    cb.set_label(label)
    cb.outline.set_linewidth(0.8)
    cb.ax.tick_params(which="both", width=0.8)
    _panel_label(ax, "b", p["font"], color="white")


def draw_violin(ax, data, p):
    parts = ax.violinplot(data, showmedians=True, showextrema=True, widths=0.8)
    for body in parts["bodies"]:
        body.set_facecolor("#d9d9d9")
        body.set_edgecolor("black")
        body.set_linewidth(0.8)
        body.set_alpha(1)
    for key in ("cbars", "cmins", "cmaxes", "cmedians"):
        parts[key].set_color("black")
        parts[key].set_linewidth(1.0)
    ax.set_xticks([1, 2, 3], [lab for _, lab, *_ in COLORS])
    ax.set_xlim(0.4, 3.6)
    ax.set_ylabel("Member color [mag]")
    ax.set_ylim(min(d.min() for d in data) - 0.15, max(d.max() for d in data) + 0.2)
    _panel_label(ax, "c", p["font"])


def draw_scatter(ax, tab, z, p):
    vals = np.vstack([np.asarray(tab[f"{c}_STD"], float) for c, *_ in COLORS])
    top = np.nanmax(vals, axis=0)
    ax.vlines(z, np.nanmin(vals, axis=0), top, color=MUTED, lw=0.9, zorder=1)
    for (col, lab, color, marker), v in zip(COLORS, vals):
        ax.scatter(z, v, s=p["marker_size"], color=color, marker=marker, label=lab,
                   edgecolor="white", linewidth=0.4, zorder=3)
    for zi, yi, n in zip(z, top, tab["N_MEMBERS"]):
        ax.annotate(str(n), (zi, yi), xytext=(0, 7), textcoords="offset points",
                    ha="center", va="bottom", rotation=90, fontsize=p["font_count"])
    ax.set_ylim(0, 1.3 * np.nanmax(top))
    ax.set_ylabel("Color scatter [mag]")
    _legend(ax, loc="upper left", bbox_to_anchor=(0.3, 1.0), ncol=3, columnspacing=0.8,
            handlelength=1.0)
    _panel_label(ax, "d", p["font"])


def draw_purity(ax, tab, z, p, show_counts=True):
    vals = np.vstack([np.asarray(tab[f"{c}_PURITY"], float) for c, *_ in COLORS])
    low = np.nanmin(vals, axis=0)
    ax.vlines(z, low, np.nanmax(vals, axis=0), color=MUTED, lw=0.9, zorder=1)
    for (col, lab, color, marker), v in zip(COLORS, vals):
        ax.scatter(z, v, s=p["marker_size"], color=color, marker=marker, label=lab,
                   edgecolor="white", linewidth=0.4, zorder=3)
    if show_counts:
        for zi, yi, nm, nn in zip(z, low, tab["N_MEMBERS"], tab["N_NONMEMBERS"]):
            ax.annotate(f"{nm}/{nn}", (zi, yi), xytext=(0, -7), textcoords="offset points",
                        ha="center", va="top", rotation=90, fontsize=p["font_count"])
        ax.set_ylim(np.nanmin(vals) - 0.45, 1.05)  # room for the labels below the points
    else:
        ax.set_ylim(p["purity_ymin"], 1.05)        # points below this are not shown
    ax.set_ylabel("Red-sequence purity")
    _panel_label(ax, "e", p["font"])


def draw_median(ax, tab, z, p):
    order = np.argsort(z)
    for col, lab, color, marker in COLORS:
        ax.errorbar(z[order], np.asarray(tab[f"{col}_MEDIAN"], float)[order],
                    yerr=np.asarray(tab[f"{col}_MEDIAN_ERR"], float)[order],
                    fmt=marker + "-", color=color, ms=3.5, lw=0.8, elinewidth=0.8,
                    capsize=1.5, mec=color, label=lab)
    ax.set_ylabel("Median color [mag]")
    ax.set_xlabel("Redshift")
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    _panel_label(ax, "f", p["font"])


def make_figure(an, violin_data, scatter_tab, purity_tab, p, out, show_counts=True):
    """Assemble the six panels on one page and save `out`.pdf / `out`.png."""
    plt.rcParams.update(style(p))

    sc = scatter_tab[~np.asarray(scatter_tab["FLAG_LOW_N"])]
    # one set of display positions for all redshift panels, keyed by cluster
    zpos = dict(zip(sc["CLUSTER"],
                    spread_positions(np.asarray(sc["Z_MEMBERS_MEAN"], float), p["min_sep"])))
    pu = purity_tab[[c in zpos for c in purity_tab["CLUSTER"]]]

    fig = plt.figure(figsize=(p["fig_width"], p["fig_height"]), layout="constrained")
    fig.get_layout_engine().set(w_pad=0.03, h_pad=0.03, wspace=0.02, hspace=0.02)
    outer = fig.add_gridspec(2, 1, height_ratios=[1, 2.05])
    top = outer[0].subgridspec(1, 2, width_ratios=[1.45, 1])
    bot = outer[1].subgridspec(1, 3, width_ratios=[1, 0.18, 3.3])  # middle column = gap
    stack = bot[2].subgridspec(3, 1)

    ax_cmd = fig.add_subplot(top[0])
    ax_map = fig.add_subplot(top[1])
    ax_vio = fig.add_subplot(bot[0])
    ax_sc = fig.add_subplot(stack[0])
    ax_pu = fig.add_subplot(stack[1], sharex=ax_sc)
    ax_md = fig.add_subplot(stack[2], sharex=ax_sc)

    draw_cmd(ax_cmd, an, p)
    draw_density(fig, ax_map, an, p)
    draw_violin(ax_vio, violin_data, p)
    z_sc = np.array([zpos[c] for c in sc["CLUSTER"]])
    draw_scatter(ax_sc, sc, z_sc, p)
    draw_purity(ax_pu, pu, np.array([zpos[c] for c in pu["CLUSTER"]]), p, show_counts)
    draw_median(ax_md, sc, z_sc, p)
    for ax in (ax_sc, ax_pu):
        ax.tick_params(labelbottom=False)
    ax_sc.set_xlim(-0.01, p["zmax"])

    # Let constrained layout place everything once, then freeze it and stretch
    # (c) so its top/bottom edges line up with the top of (d) and bottom of (f).
    fig.canvas.draw()
    fig.set_layout_engine("none")
    top_d, bot_f, pv = ax_sc.get_position(), ax_md.get_position(), ax_vio.get_position()
    ax_vio.set_position([pv.x0, bot_f.y0, pv.width, top_d.y1 - bot_f.y0])

    fig.savefig(out + ".pdf")
    fig.savefig(out + ".png", dpi=p["dpi"])
    plt.close(fig)
    print(f"Saved {out}.pdf / .png")
