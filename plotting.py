"""Plot originali MeteoD estratti dallo script Leonardo, senza accesso SSH/R2.
Alcune forme della grafica (tercili 40/20/40 riscalati) replicano esattamente l'originale.
"""
import os
from datetime import datetime,timedelta
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker
from matplotlib.patches import Patch,Rectangle
from matplotlib.colors import to_rgba
try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
except ImportError:
    ccrs = cfeature = None  # Only city charts can run without Cartopy.
TARGET_DATE = datetime(2026,10,22)
def infer_year_from_mmdd(mmdd_str, ref_date):
    m = int(mmdd_str[:2])
    d = int(mmdd_str[2:4])

    try:
        candidate = datetime(ref_date.year, m, d)
    except ValueError:
        return ref_date

    delta_days = (candidate - ref_date).days
    if abs(delta_days) > 180:
        if delta_days > 0:
            candidate = candidate.replace(year=ref_date.year - 1)
        else:
            candidate = candidate.replace(year=ref_date.year + 1)

    return candidate

def get_discrete_color(probability, tercile_index, varname):
    if tercile_index == 0:
        thresholds = [35, 50, 65, 80, 90]
        if varname == "2tp":
            colors = ["#cce5ff", "#99ccff", "#66b2ff", "#3399ff", "#0059b3"]
        else:
            colors = ["#ffe0cc", "#ffb380", "#ff944d", "#ff751a", "#994d00"]
    elif tercile_index == 1:
        thresholds = [35, 50, 65]
        colors = ["#e5e5e5", "#f5f5f5", "#ffffff"]
    elif tercile_index == 2:
        thresholds = [35, 50, 65, 80, 90]
        if varname == "2tp":
            colors = ["#ffcccc", "#ff9999", "#ff6666", "#ff3333", "#990000"]
        else:
            colors = ["#d9f2d9", "#a6d8a8", "#73bf78", "#409f49", "#1f5c33"]
    else:
        return "#ffffff"

    if probability >= 100:
        return colors[-1]

    for t, c in zip(thresholds[::-1], colors[::-1]):
        if probability >= t:
            return c
    return "#ffffff"

def add_colored_hatching(ax_, mask, lon2d_, lat2d_, color_):
    cs = ax_.contourf(
        lon2d_, lat2d_, mask,
        levels=[0.5, 1.5],
        colors="none",
        hatches=["////"],
        transform=ccrs.PlateCarree()
    )
    for coll in getattr(cs, "collections", [cs]):
        coll.set_edgecolor(color_)
        coll.set_linewidth(0.15)

def plot_terciles_week(week_key, array, varname, output_dir, lat2d, lon2d):
    tercile_data = array[10:13, :, :].copy()
    tercile_data[1, :, :] *= 2

    sum_terciles = np.sum(tercile_data, axis=0)
    sum_terciles[sum_terciles == 0] = 1
    tercile_data = tercile_data / sum_terciles * 100

    max_idx = np.argmax(tercile_data, axis=0)

    color_map = np.full(max_idx.shape, "#ffffff", dtype=object)
    for i in range(max_idx.shape[0]):
        for j in range(max_idx.shape[1]):
            idx = max_idx[i, j]
            prob_normalized = tercile_data[idx, i, j]
            color_map[i, j] = get_discrete_color(prob_normalized, idx, varname)

    fig = plt.figure(figsize=(12, 8))
    gs = gridspec.GridSpec(2, 1, height_ratios=[0.9, 0.3], hspace=0.05)

    map_gs = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs[0])
    ax = fig.add_subplot(map_gs[0], projection=ccrs.PlateCarree())

    ax.set_extent([-12, 42, 29.2, 52], crs=ccrs.PlateCarree())
    ax.coastlines(resolution="10m", linewidth=0.5)
    ax.add_feature(cfeature.BORDERS, linewidth=0.8)

    gl = ax.gridlines(draw_labels=True, linestyle="--", linewidth=0.5, color="gray")
    gl.top_labels = False
    gl.right_labels = False
    gl.xlocator = mticker.FixedLocator(range(-180, 181, 10))
    gl.ylocator = mticker.FixedLocator(range(-90, 91, 5))

    rgba_map = np.array([[to_rgba(c) for c in row] for row in color_map])

    ax.imshow(
        rgba_map,
        origin="upper",
        extent=(lon2d.min(), lon2d.max(), lat2d.min(), lat2d.max()),
        transform=ccrs.PlateCarree()
    )

    p10 = array[0, :, :]
    p90 = array[4, :, :]

    thresh_10 = 80
    thresh_90 = 80

    mask_10 = p10 >= thresh_10
    mask_90 = p90 >= thresh_90

    contour_10 = mask_10.astype(float)
    contour_90 = mask_90.astype(float)
    levels = [0.5, 1.5]

    if varname == "2tp":
        color_10 = "#66b2ff"
        color_90 = "#ff57a5"
    else:
        color_10 = "#a30021"
        color_90 = "#009300"

    ax.contour(
        lon2d, lat2d, contour_10,
        levels=levels, colors=color_10,
        linewidths=0.5,
        transform=ccrs.PlateCarree()
    )
    ax.contour(
        lon2d, lat2d, contour_90,
        levels=levels, colors=color_90,
        linewidths=0.5,
        transform=ccrs.PlateCarree()
    )

    add_colored_hatching(ax, contour_10, lon2d, lat2d, color_10)
    add_colored_hatching(ax, contour_90, lon2d, lat2d, color_90)

    try:
        ref = globals().get("TARGET_DATE", datetime.today())

        datecode = week_key.split("____")[0][-4:]
        temp_date = infer_year_from_mmdd(datecode, ref)
        end_date = temp_date - timedelta(days=1)
        start_date = end_date - timedelta(days=6)

        run_date_string = week_key[3:7]
        run_date = infer_year_from_mmdd(run_date_string, ref)

        subtitle_date = f"{start_date.strftime('%a %d %b %Y')} - {end_date.strftime('%a %d %b %Y')}"
        run_info = f"(run {run_date.strftime('%d/%m/%Y')})"
        subtitle = f"{subtitle_date} {run_info}"
    except Exception:
        subtitle = ""
        start_date = datetime.today()

    main_title_text = (
        "Weekly Anomaly Probability - Temperature 2m"
        if varname == "2tp"
        else "Weekly Anomaly Probability - Total Precipitation"
    )

    full_title = r"$\bf{" + main_title_text.replace(" ", r"\ ") + r"}$" + f"\n{subtitle}"
    ax.set_title(full_title, fontsize=12, loc="center", pad=10, linespacing=1.5)

    bottom_gs = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs[1])
    legend_container_ax = fig.add_subplot(bottom_gs[0])
    legend_container_ax.set_axis_off()

    legend_outer = gridspec.GridSpecFromSubplotSpec(
        2, 3, subplot_spec=bottom_gs[0],
        height_ratios=[0.25, 0.75], hspace=0.4, wspace=0.2
    )

    legend_title_ax = fig.add_subplot(legend_outer[0, :])
    legend_title_ax.text(
        0.5, 0.5, "Most likely scenario",
        fontsize=10, fontweight="bold", ha="center", va="center"
    )
    legend_title_ax.axis("off")

    if varname == "2tp":
        colors_1st = ["#cce5ff", "#99ccff", "#66b2ff", "#3399ff", "#0059b3"]
        colors_3rd = ["#ffcccc", "#ff9999", "#ff6666", "#ff3333", "#990000"]
    else:
        colors_1st = ["#ffe0cc", "#ffb380", "#ff944d", "#ff751a", "#994d00"]
        colors_3rd = ["#d9f2d9", "#a6d8a8", "#73bf78", "#409f49", "#1f5c33"]

    colors_mid = ["#e5e5e5", "#f5f5f5", "#ffffff"]
    probs_mid = [35, 50, 65]

    legend_cols_info = {
        0: {"title": "Below average", "colors": colors_1st, "probs": [35, 50, 65, 80, 90]},
        1: {"title": "Near average", "colors": colors_mid, "probs": probs_mid},
        2: {"title": "Above average", "colors": colors_3rd, "probs": [35, 50, 65, 80, 90]},
    }

    legend_axes = [legend_title_ax]

    for i, info in legend_cols_info.items():
        ax_col = fig.add_subplot(legend_outer[1, i])
        patches = [
            Patch(facecolor=c, edgecolor="k", label=f"> {p}%")
            for c, p in zip(info["colors"], info["probs"])
        ]
        ax_col.legend(
            handles=patches,
            loc="upper center",
            fontsize=8,
            frameon=True,
            labelspacing=0.8
        )
        ax_col.set_title(info["title"], fontsize=9)
        ax_col.axis("off")
        legend_axes.append(ax_col)

        if i == 1:
            ax_col.text(
                0.5, 0.18, "Extreme events",
                ha="center", va="top",
                fontsize=9, fontweight="bold",
                transform=ax_col.transAxes
            )
            contour_patches = [
                Patch(facecolor="none", edgecolor=color_10, hatch="////", label=f"20th percentile > {thresh_10}%"),
                Patch(facecolor="none", edgecolor=color_90, hatch="////", label=f"80th percentile > {thresh_90}%"),
            ]
            inset_ax = ax_col.inset_axes([0.15, -0.40, 0.7, 0.42], transform=ax_col.transAxes)
            inset_ax.legend(handles=contour_patches, loc="center", fontsize=8, frameon=True)
            inset_ax.axis("off")
            legend_axes.append(inset_ax)

    p = legend_container_ax.get_position(fig)
    target_width = 0.60
    cx = 0.5 * (p.x0 + p.x1)
    new_half_w = 0.5 * target_width
    new_x0 = cx - new_half_w
    legend_container_ax.set_position([new_x0, p.y0, 2 * new_half_w, p.height])

    bboxes = [ax_.get_position(fig) for ax_ in legend_axes]
    x0_raw = min(b.x0 for b in bboxes)
    y0_raw = min(b.y0 for b in bboxes)
    x1_raw = max(b.x1 for b in bboxes)
    y1_raw = max(b.y1 for b in bboxes)

    padding_x = 0.01
    padding_y = 0.03
    offset_y = 0.02

    legend_width = (x1_raw - x0_raw) + 2 * padding_x
    legend_height = (y1_raw - y0_raw) + 2 * padding_y

    x0_centered = 0.5 - legend_width / 2.0
    y0 = y0_raw - padding_y - offset_y

    legend_border = Rectangle(
        (x0_centered, y0),
        legend_width,
        legend_height,
        linewidth=0.5,
        edgecolor="black",
        facecolor="none",
        transform=fig.transFigure,
        zorder=1000
    )
    fig.patches.append(legend_border)

    os.makedirs(output_dir, exist_ok=True)
    filename = f"subseasonal_{start_date.strftime('%Y%m%d')}_{varname}.webp"
    plt.savefig(
        os.path.join(output_dir, filename),
        format="webp",
        dpi=150,
        bbox_inches="tight",
        pil_kwargs={"quality": 70}
    )
    plt.close(fig)

def get_dominant_quintile_simple(probs):
    probs = np.array(probs)
    max_idx = int(np.argmax(probs))
    max_val = probs[max_idx]
    return max_idx if max_val > 30 else None

def plot_city_combined(city, weeks, q_data, s_data, out_path, run_date):
    n = len(weeks)
    x = np.arange(n)

    colors_T = ["#990000", "#ff9999", "#ffffff", "#99ccff", "#004487"][::-1]
    colors_P = ["#1f5c33", "#a6d8a8", "#ffffff", "#ffb380", "#994d00"][::-1]
    colors_T_leg = ["#990000", "#ff9999", "#ffffff", "#99ccff", "#004487"]
    colors_P_leg = ["#1f5c33", "#a6d8a8", "#ffffff", "#ffb380", "#994d00"]

    labels_T = ["Più caldo", "Sopra media", "In media", "Sotto media", "Più freddo"]
    labels_P = ["Più piovoso", "Sopra media", "In media", "Sotto media", "Più secco"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 9))

    week_labels = [datetime.strptime(w, "%Y%m%d").strftime("%d/%m/%Y") for w in weeks]

    q_vals_T = np.array([[q_data[w]["T_1"], q_data[w]["T_2"], q_data[w]["T_3"],
                          q_data[w]["T_4"], q_data[w]["T_5"]] for w in weeks])
    q_vals_P = np.array([[q_data[w]["P_1"], q_data[w]["P_2"], q_data[w]["P_3"],
                          q_data[w]["P_4"], q_data[w]["P_5"]] for w in weeks])

    dominant_T = [get_dominant_quintile_simple(q_vals_T[i, :]) for i in range(n)]
    dominant_P = [get_dominant_quintile_simple(q_vals_P[i, :]) for i in range(n)]

    bottoms1 = np.zeros(n)
    for i in range(5):
        linewidths = [3 if dominant_T[j] == i else 1 for j in range(n)]
        bars = ax1.bar(x, q_vals_T[:, i], bottom=bottoms1, color=colors_T[i],
                       edgecolor="black", linewidth=linewidths)

        for j, bar in enumerate(bars):
            height = q_vals_T[j, i]
            if height > 3:
                center_y = bottoms1[j] + height / 2
                ax1.text(bar.get_x() + bar.get_width()/2, center_y, f"{height:.0f}%",
                         ha="center", va="center", fontsize=8, fontweight="normal", color="black")
        bottoms1 += q_vals_T[:, i]

    thr_prefix_T = "T_p"
    thr_levels = [0, 20, 40, 60, 80, 100]

    for i, w in enumerate(weeks):
        xi = x[i]
        thr_vals = [s_data[w][f"{thr_prefix_T}{lvl}"] for lvl in thr_levels]

        cum_probs = np.cumsum(q_vals_T[i, :])
        cum_probs = np.insert(cum_probs, 0, 0)

        max_quintile_idx = dominant_T[i]

        for j, (y_thr, v_thr) in enumerate(zip(cum_probs, thr_vals)):
            is_highlighted = False
            if max_quintile_idx is not None:
                if max_quintile_idx == 0:
                    is_highlighted = (j == 0) or (j == 1)
                elif max_quintile_idx == 4:
                    is_highlighted = (j == 4) or (j == 5)
                else:
                    is_highlighted = (j == max_quintile_idx) or (j == max_quintile_idx + 1)

            ax1.text(
                xi, y_thr, f"{v_thr:.1f}°C",
                ha="center", va="bottom", fontsize=9,
                fontweight="bold" if is_highlighted else "normal",
                color="black",
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.8)
            )

    ax1.set_xticks(x)
    ax1.set_xticklabels(week_labels, rotation=0, fontweight="bold")
    ax1.set_ylim(0, 100)
    ax1.set_ylabel("Probabilità (%)")
    ax1.set_title("Temperatura media", fontsize=12, fontweight="bold", pad=20)
    ax1.set_xlabel("Settimane")

    bottoms2 = np.zeros(n)
    for i in range(5):
        linewidths = [3 if dominant_P[j] == i else 1 for j in range(n)]
        bars = ax2.bar(x, q_vals_P[:, i], bottom=bottoms2, color=colors_P[i],
                       edgecolor="black", linewidth=linewidths)

        for j, bar in enumerate(bars):
            height = q_vals_P[j, i]
            if height > 3:
                center_y = bottoms2[j] + height / 2
                ax2.text(bar.get_x() + bar.get_width()/2, center_y, f"{height:.0f}%",
                         ha="center", va="center", fontsize=8, fontweight="normal", color="black")
        bottoms2 += q_vals_P[:, i]

    thr_prefix_P = "P_p"
    for i, w in enumerate(weeks):
        xi = x[i]
        thr_vals = [s_data[w][f"{thr_prefix_P}{lvl}"] for lvl in thr_levels]

        cum_probs = np.cumsum(q_vals_P[i, :])
        cum_probs = np.insert(cum_probs, 0, 0)

        max_quintile_idx = dominant_P[i]

        for j, (y_thr, v_thr) in enumerate(zip(cum_probs, thr_vals)):
            is_highlighted = False
            if max_quintile_idx is not None:
                if max_quintile_idx == 0:
                    is_highlighted = (j == 0) or (j == 1)
                elif max_quintile_idx == 4:
                    is_highlighted = (j == 4) or (j == 5)
                else:
                    is_highlighted = (j == max_quintile_idx) or (j == max_quintile_idx + 1)

            ax2.text(
                xi, y_thr, f"{v_thr:.1f}mm",
                ha="center", va="bottom", fontsize=9,
                fontweight="bold" if is_highlighted else "normal",
                color="black",
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.8)
            )

    ax2.set_xticks(x)
    ax2.set_xticklabels(week_labels, rotation=0, fontweight="bold")
    ax2.set_ylim(0, 100)
    ax2.set_ylabel("Probabilità (%)")
    ax2.set_title("Precipitazione totale", fontsize=12, fontweight="bold", pad=20)
    ax2.set_xlabel("Settimane")

    plt.suptitle(f"{city.upper()} - Quintili settimanali - RUN: {run_date}",
                 fontsize=14, fontweight="bold")

    legend_elements_T = [Patch(facecolor=colors_T_leg[i], edgecolor="black", linewidth=1.2,
                               label=labels_T[i]) for i in range(5)]
    legend_elements_P = [Patch(facecolor=colors_P_leg[i], edgecolor="black", linewidth=1.2,
                               label=labels_P[i]) for i in range(5)]

    fig.legend(legend_elements_T, labels_T, loc="lower center",
               bbox_to_anchor=(0.25, -0.05), ncol=1, fontsize=9)

    fig.legend(legend_elements_P, labels_P, loc="lower center",
               bbox_to_anchor=(0.75, -0.05), ncol=1, fontsize=9)

    plt.subplots_adjust(bottom=0.15, top=0.88)
    plt.savefig(out_path, format="webp", dpi=150, bbox_inches="tight",
                pil_kwargs={"quality": 70})
    plt.close(fig)
