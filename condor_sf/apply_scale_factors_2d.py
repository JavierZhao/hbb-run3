import warnings
warnings.filterwarnings('ignore')
import os
import numpy as np
from coffea import util

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import mplhep as hep
plt.style.use(hep.style.ROOT)

import matplotlib.pylab as pylab
params = {'legend.fontsize': 'medium',
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'medium',
         'ytick.labelsize':'medium'}
pylab.rcParams.update(params)

import matplotlib as mpl
mpl.rcParams['lines.linewidth'] = 5
import json
import argparse


def load_coffea_outputs(year, prod_mode, base_dir=None):
    """
    Load the coffea output files for data and MC scale factors.

    Parameters:
    - year: str - Year of the data
    - prod_mode: str - Production mode (VBF or ggF)
    - base_dir: str - Base directory for outputs (defaults to current working directory)

    Returns:
    - output_data: dict - Data output
    - output_mc: dict - MC output
    """
    if base_dir is None:
        base_dir = os.getcwd()

    output_dir = os.path.join(base_dir, "output")

    data_file = os.path.join(output_dir, f"sf_2d_{year}_{prod_mode}_data.coffea")
    mc_file = os.path.join(output_dir, f"sf_2d_{year}_{prod_mode}_mc.coffea")

    print(f"Loading data output from: {data_file}")
    output_data = util.load(data_file)

    print(f"Loading MC output from: {mc_file}")
    output_mc = util.load(mc_file)

    return output_data, output_mc


def calculate_2d_efficiencies_and_sf(output_data, output_mc, var_name, bins_x, bins_y):
    """
    Calculate 2D efficiencies for data and MC, and the scale factor.

    Parameters:
    - output_data: dict - Data output
    - output_mc: dict - MC output
    - var_name: str - Variable name (e.g., 'pt_vs_msd')
    - bins_x: array - Bin edges for x-axis
    - bins_y: array - Bin edges for y-axis

    Returns:
    - eff_data: 2D array - Data efficiency
    - eff_mc: 2D array - MC efficiency
    - sf: 2D array - Scale factor
    """
    # Get baseline (total) data for both data and MC
    baseline_x_data = np.array(output_data['Baseline'][var_name]['x'])
    baseline_y_data = np.array(output_data['Baseline'][var_name]['y'])
    baseline_x_mc = np.array(output_mc['Baseline'][var_name]['x'])
    baseline_y_mc = np.array(output_mc['Baseline'][var_name]['y'])

    # Get numerator data for both data and MC
    pass_x_data = np.array(output_data['Numerator'][var_name]['pass_x'])
    pass_y_data = np.array(output_data['Numerator'][var_name]['pass_y'])
    pass_x_mc = np.array(output_mc['Numerator'][var_name]['pass_x'])
    pass_y_mc = np.array(output_mc['Numerator'][var_name]['pass_y'])

    # Create 2D histograms for baseline (denominator)
    hist_baseline_data, _, _ = np.histogram2d(
        baseline_x_data, baseline_y_data,
        bins=[bins_x, bins_y]
    )
    hist_baseline_mc, _, _ = np.histogram2d(
        baseline_x_mc, baseline_y_mc,
        bins=[bins_x, bins_y]
    )

    # Create 2D histograms for passing events (numerator)
    hist_pass_data, _, _ = np.histogram2d(
        pass_x_data, pass_y_data,
        bins=[bins_x, bins_y]
    )
    hist_pass_mc, _, _ = np.histogram2d(
        pass_x_mc, pass_y_mc,
        bins=[bins_x, bins_y]
    )

    # Calculate efficiencies
    with np.errstate(divide='ignore', invalid='ignore'):
        eff_data = np.divide(
            hist_pass_data, hist_baseline_data,
            out=np.zeros_like(hist_pass_data, dtype=float),
            where=hist_baseline_data > 0
        )
        eff_mc = np.divide(
            hist_pass_mc, hist_baseline_mc,
            out=np.zeros_like(hist_pass_mc, dtype=float),
            where=hist_baseline_mc > 0
        )

    # Calculate scale factor
    valid_sf = (hist_baseline_data > 0) & (hist_baseline_mc > 0) & (eff_mc > 0)
    sf = np.divide(
        eff_data, eff_mc,
        out=np.ones_like(eff_data, dtype=float),
        where=valid_sf
    )

    # Mask invalid SF bins
    sf[~valid_sf] = np.nan

    return eff_data, eff_mc, sf


def plot_scaled_efficiency_comparison(eff_data, eff_mc, sf, bins_x, bins_y,
                                      label_x, label_y, save_path, year, prod_mode):
    """
    Plot comparison of data efficiency, MC efficiency, and MC*SF efficiency.

    Parameters:
    - eff_data: 2D array - Data efficiency
    - eff_mc: 2D array - MC efficiency
    - sf: 2D array - Scale factor
    - bins_x, bins_y: arrays - Bin edges
    - label_x, label_y: str - Axis labels
    - save_path: str - Path to save the figure
    - year: str - Year label
    - prod_mode: str - Production mode
    """
    # Calculate scaled MC efficiency
    eff_mc_scaled = eff_mc * sf

    # Calculate overall efficiencies
    total_eff_data = np.nanmean(eff_data[eff_data > 0]) if np.any(eff_data > 0) else 0.0
    total_eff_mc = np.nanmean(eff_mc[eff_mc > 0]) if np.any(eff_mc > 0) else 0.0
    total_eff_mc_scaled = np.nanmean(eff_mc_scaled[eff_mc_scaled > 0]) if np.any(eff_mc_scaled > 0) else 0.0

    # Create figure with 4 subplots (2x2)
    fig, axes = plt.subplots(2, 2, figsize=(20, 16))

    # Plot 1: Data Efficiency
    ax1 = axes[0, 0]
    im1 = ax1.pcolormesh(
        bins_x, bins_y, eff_data.T,
        cmap='RdYlGn',
        vmin=0, vmax=1,
        shading='flat'
    )
    cbar1 = plt.colorbar(im1, ax=ax1)
    cbar1.set_label('Efficiency', fontsize=14)
    ax1.set_xlabel(label_x, fontsize=14)
    ax1.set_ylabel(label_y, fontsize=14)
    ax1.set_title(f"Data Efficiency\nOverall: {total_eff_data:.2%}", fontsize=14)
    hep.cms.label(ax=ax1, data=True, year=year, com="13.6", fontsize=12)

    # Plot 2: MC Efficiency (Uncorrected)
    ax2 = axes[0, 1]
    im2 = ax2.pcolormesh(
        bins_x, bins_y, eff_mc.T,
        cmap='RdYlGn',
        vmin=0, vmax=1,
        shading='flat'
    )
    cbar2 = plt.colorbar(im2, ax=ax2)
    cbar2.set_label('Efficiency', fontsize=14)
    ax2.set_xlabel(label_x, fontsize=14)
    ax2.set_ylabel(label_y, fontsize=14)
    ax2.set_title(f"MC Efficiency (Uncorrected)\nOverall: {total_eff_mc:.2%}", fontsize=14)
    hep.cms.label(ax=ax2, data=False, year=year, com="13.6", label="Simulation", fontsize=12)

    # Plot 3: MC Efficiency (Scaled by SF)
    ax3 = axes[1, 0]
    im3 = ax3.pcolormesh(
        bins_x, bins_y, eff_mc_scaled.T,
        cmap='RdYlGn',
        vmin=0, vmax=1,
        shading='flat'
    )
    cbar3 = plt.colorbar(im3, ax=ax3)
    cbar3.set_label('Efficiency', fontsize=14)
    ax3.set_xlabel(label_x, fontsize=14)
    ax3.set_ylabel(label_y, fontsize=14)
    ax3.set_title(f"MC Efficiency × SF\nOverall: {total_eff_mc_scaled:.2%}", fontsize=14)
    hep.cms.label(ax=ax3, data=False, year=year, com="13.6", label="Simulation (Corrected)", fontsize=12)

    # Plot 4: Scale Factor
    ax4 = axes[1, 1]
    im4 = ax4.pcolormesh(
        bins_x, bins_y, sf.T,
        cmap='RdBu_r',
        vmin=0.8, vmax=1.2,
        shading='flat'
    )
    cbar4 = plt.colorbar(im4, ax=ax4)
    cbar4.set_label('Scale Factor (Data/MC)', fontsize=14)
    ax4.set_xlabel(label_x, fontsize=14)
    ax4.set_ylabel(label_y, fontsize=14)
    mean_sf = np.nanmean(sf[np.isfinite(sf)])
    ax4.set_title(f"Scale Factor\nMean: {mean_sf:.3f}", fontsize=14)
    hep.cms.label(ax=ax4, data=True, year=year, com="13.6", fontsize=12)

    plt.suptitle(f"Trigger Efficiency Comparison - {prod_mode} - {year}", fontsize=16, y=0.995)
    plt.tight_layout()

    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.savefig(save_path, dpi=200, bbox_inches='tight')
    print(f"Saved scaled efficiency comparison plot to: {save_path}")
    plt.close(fig)


def plot_1d_projections(eff_data, eff_mc, sf, bins_x, bins_y,
                        label_x, label_y, save_dir, year, prod_mode):
    """
    Plot 1D projections of efficiencies along x and y axes.

    Parameters:
    - eff_data: 2D array - Data efficiency
    - eff_mc: 2D array - MC efficiency
    - sf: 2D array - Scale factor
    - bins_x, bins_y: arrays - Bin edges
    - label_x, label_y: str - Axis labels
    - save_dir: str - Directory to save figures
    - year: str - Year label
    - prod_mode: str - Production mode
    """
    # Calculate scaled MC efficiency
    eff_mc_scaled = eff_mc * sf

    # Project onto x-axis (average over y)
    eff_data_x = np.nanmean(eff_data, axis=1)
    eff_mc_x = np.nanmean(eff_mc, axis=1)
    eff_mc_scaled_x = np.nanmean(eff_mc_scaled, axis=1)

    # Project onto y-axis (average over x)
    eff_data_y = np.nanmean(eff_data, axis=0)
    eff_mc_y = np.nanmean(eff_mc, axis=0)
    eff_mc_scaled_y = np.nanmean(eff_mc_scaled, axis=0)

    # Plot 1: X-axis projection
    fig, (ax, rax) = plt.subplots(
        2, 1, figsize=(10, 10),
        gridspec_kw={"height_ratios": (3, 1)},
        sharex=True
    )
    fig.subplots_adjust(hspace=0.05)

    # Top panel: efficiencies
    bin_centers_x = (bins_x[:-1] + bins_x[1:]) / 2
    ax.step(bins_x[:-1], eff_data_x, where='post', color='black', linewidth=2, label='Data')
    ax.step(bins_x[:-1], eff_mc_x, where='post', color='blue', linewidth=2, linestyle='--', label='MC (Uncorrected)')
    ax.step(bins_x[:-1], eff_mc_scaled_x, where='post', color='red', linewidth=2, linestyle=':', label='MC × SF')
    ax.set_ylabel("Efficiency", fontsize=14)
    ax.set_ylim(0.0, 1.1)
    ax.legend(fontsize=12, loc='best')
    ax.set_title(f"Efficiency vs {label_x} - {prod_mode} - {year}", fontsize=14)
    hep.cms.label("Preliminary", data=True, year=year, ax=ax)
    ax.grid(True, alpha=0.3)

    # Bottom panel: ratio to data
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio_mc = np.divide(eff_mc_x, eff_data_x, out=np.full_like(eff_mc_x, np.nan), where=(eff_data_x > 0))
        ratio_mc_scaled = np.divide(eff_mc_scaled_x, eff_data_x, out=np.full_like(eff_mc_scaled_x, np.nan), where=(eff_data_x > 0))

    rax.step(bins_x[:-1], ratio_mc, where='post', color='blue', linewidth=2, linestyle='--', label='MC / Data')
    rax.step(bins_x[:-1], ratio_mc_scaled, where='post', color='red', linewidth=2, linestyle=':', label='(MC × SF) / Data')
    rax.axhline(1.0, color='black', linestyle='-', linewidth=1)
    rax.set_ylabel("Ratio to Data", fontsize=14)
    rax.set_xlabel(label_x, fontsize=14)
    rax.set_ylim(0.8, 1.2)
    rax.legend(fontsize=10, loc='best')
    rax.grid(True, alpha=0.3)

    save_path = os.path.join(save_dir, f"scaled_eff_projection_x_{prod_mode}.png")
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.savefig(save_path, dpi=200, bbox_inches='tight')
    print(f"Saved x-projection plot to: {save_path}")
    plt.close(fig)

    # Plot 2: Y-axis projection
    fig, (ax, rax) = plt.subplots(
        2, 1, figsize=(10, 10),
        gridspec_kw={"height_ratios": (3, 1)},
        sharex=True
    )
    fig.subplots_adjust(hspace=0.05)

    # Top panel: efficiencies
    bin_centers_y = (bins_y[:-1] + bins_y[1:]) / 2
    ax.step(bins_y[:-1], eff_data_y, where='post', color='black', linewidth=2, label='Data')
    ax.step(bins_y[:-1], eff_mc_y, where='post', color='blue', linewidth=2, linestyle='--', label='MC (Uncorrected)')
    ax.step(bins_y[:-1], eff_mc_scaled_y, where='post', color='red', linewidth=2, linestyle=':', label='MC × SF')
    ax.set_ylabel("Efficiency", fontsize=14)
    ax.set_ylim(0.0, 1.1)
    ax.legend(fontsize=12, loc='best')
    ax.set_title(f"Efficiency vs {label_y} - {prod_mode} - {year}", fontsize=14)
    hep.cms.label("Preliminary", data=True, year=year, ax=ax)
    ax.grid(True, alpha=0.3)

    # Bottom panel: ratio to data
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio_mc = np.divide(eff_mc_y, eff_data_y, out=np.full_like(eff_mc_y, np.nan), where=(eff_data_y > 0))
        ratio_mc_scaled = np.divide(eff_mc_scaled_y, eff_data_y, out=np.full_like(eff_mc_scaled_y, np.nan), where=(eff_data_y > 0))

    rax.step(bins_y[:-1], ratio_mc, where='post', color='blue', linewidth=2, linestyle='--', label='MC / Data')
    rax.step(bins_y[:-1], ratio_mc_scaled, where='post', color='red', linewidth=2, linestyle=':', label='(MC × SF) / Data')
    rax.axhline(1.0, color='black', linestyle='-', linewidth=1)
    rax.set_ylabel("Ratio to Data", fontsize=14)
    rax.set_xlabel(label_y, fontsize=14)
    rax.set_ylim(0.8, 1.2)
    rax.legend(fontsize=10, loc='best')
    rax.grid(True, alpha=0.3)

    save_path = os.path.join(save_dir, f"scaled_eff_projection_y_{prod_mode}.png")
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.savefig(save_path, dpi=200, bbox_inches='tight')
    print(f"Saved y-projection plot to: {save_path}")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Apply 2D Scale Factors to MC Efficiencies")
    parser.add_argument('--year', required=False, choices=['2022', '2022EE', '2023', '2023BPix'], default='2023',
                        help="Year of the dataset (e.g., 2022, 2023, etc.).")
    parser.add_argument('--prod_mode', required=False, choices=['VBF', 'ggF', 'both'], default='both',
                        help="Production mode (VBF, ggF, or both).")
    parser.add_argument('--input_dir', required=False, default=None,
                        help="Base directory containing input coffea files (default: current directory)")
    args = parser.parse_args()

    year = args.year
    input_dir = args.input_dir if args.input_dir else os.getcwd()

    # Determine which production modes to process
    if args.prod_mode == 'both':
        prod_modes = ['VBF', 'ggF']
    else:
        prod_modes = [args.prod_mode]

    print(f"\n{'='*60}")
    print(f"Applying Scale Factors - Year: {year}")
    print(f"{'='*60}\n")

    # Variable configuration (matching scale_factor_2d.py)
    var_name = 'pt_vs_msd'
    label_x = "Leading Jet $m_{SD}$ [GeV]"
    label_y = "Leading Jet $p_{T}$ [GeV]"
    bins_x = np.linspace(40, 200, 16)  # 15 bins
    bins_y = np.linspace(200, 1000, 21)  # 20 bins

    for prod_mode in prod_modes:
        print(f"\n{'='*60}")
        print(f"Processing production mode: {prod_mode}")
        print(f"{'='*60}\n")

        try:
            # Load coffea outputs
            output_data, output_mc = load_coffea_outputs(year, prod_mode, base_dir=input_dir)

            # Calculate efficiencies and scale factors
            print("Calculating efficiencies and scale factors...")
            eff_data, eff_mc, sf = calculate_2d_efficiencies_and_sf(
                output_data, output_mc, var_name, bins_x, bins_y
            )

            # Create output directory
            save_dir = os.path.join(input_dir, "figures_scaled_sf_2d", year, prod_mode)
            os.makedirs(save_dir, exist_ok=True)

            # Plot 2D comparison
            print("Creating 2D comparison plot...")
            save_path_2d = os.path.join(save_dir, f"scaled_efficiency_comparison_{var_name}_2d.png")
            plot_scaled_efficiency_comparison(
                eff_data, eff_mc, sf, bins_x, bins_y,
                label_x, label_y, save_path_2d, year, prod_mode
            )

            # Plot 1D projections
            print("Creating 1D projection plots...")
            plot_1d_projections(
                eff_data, eff_mc, sf, bins_x, bins_y,
                label_x, label_y, save_dir, year, prod_mode
            )

            print(f"\nAll plots for {prod_mode} saved to {save_dir}")

        except FileNotFoundError as e:
            print(f"Error: Could not find input files for {prod_mode}")
            print(f"  {e}")
            print(f"  Skipping {prod_mode}...")
            continue
        except Exception as e:
            print(f"Error processing {prod_mode}: {e}")
            import traceback
            traceback.print_exc()
            continue

    print(f"\n{'='*60}")
    print(f"Completed applying scale factors for all production modes")
    print(f"{'='*60}")
