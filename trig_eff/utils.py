import warnings
warnings.filterwarnings('ignore')
import os
import awkward as ak
import uproot
import hist
import numpy as np
from coffea import processor
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema, BaseSchema
from coffea.processor import dict_accumulator, list_accumulator

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

#line thickness
import matplotlib as mpl
mpl.rcParams['lines.linewidth'] = 5
import itertools
import json

# for making fancy 2d histograms
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def plot_1d_trigger_efficiencies(trig_name, prod_mode, output, triggers_set, trig_vars, baseline_key, max_triggers_per_plot=5):
    # Create the directory to save figures
    save_dir = f"/srv/figures/{baseline_key}_baseline/{prod_mode}/1d_trigger_efficiencies"
    os.makedirs(save_dir, exist_ok=True)

    # Function to chunk the triggers_set
    def chunk_triggers(triggers, chunk_size):
        for i in range(0, len(triggers), chunk_size):
            yield triggers[i:i + chunk_size]

    # Define colors for plotting
    colors = ['blue', 'green', 'red', 'purple', 'orange', 'cyan', 'magenta', 'yellow', 'black', 'brown']

    # Chunk the triggers_set into parts of size max_triggers_per_plot
    trigger_chunks = list(chunk_triggers(triggers_set, max_triggers_per_plot))

    # Loop over each chunk of triggers
    for chunk_index, trigger_chunk in enumerate(trigger_chunks):
        num_vars = len(trig_vars)
        # Calculate the number of rows and columns
        nrows = 2
        ncols = (num_vars + 1) // 2  # Ensure all variables are covered
        fig_width = 8 * ncols  # Increase width per subplot
        fig_height = 6 * nrows  # Increase height per subplot
        fig, axs = plt.subplots(nrows, ncols, figsize=(fig_width, fig_height), squeeze=False)

        # Flatten axs array for easy indexing
        axs = axs.flatten()

        # Loop over variables to create subplots
        for var_idx, var_name in enumerate(trig_vars.keys()):
            ax = axs[var_idx]

            # Prepare bin edges for the variable
            bin_edges = trig_vars[var_name]['axis'].edges

            # Loop over triggers to plot efficiencies
            for i, trigger in enumerate(trigger_chunk):
                # Retrieve per-event data
                total_values = np.array(output[trigger][var_name]['total'])
                pass_values = np.array(output[trigger][var_name]['pass'])

                # Create histograms from per-event data
                total_hist, _ = np.histogram(total_values, bins=bin_edges)
                pass_hist, _ = np.histogram(pass_values, bins=bin_edges)

                # Compute efficiency
                with np.errstate(divide='ignore', invalid='ignore'):
                    efficiency = np.nan_to_num(pass_hist / total_hist, nan=0.0, posinf=0.0, neginf=0.0)

                # Compute bin centers for plotting
                bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

                # Plot efficiency
                ax.step(bin_centers, efficiency, where='mid', label=f'{trigger}', color=colors[i % len(colors)], linewidth=2)

            # Set axis labels and title with smaller font size
            ax.set_xlabel(trig_vars[var_name]['label'], fontsize=12)
            ax.set_ylabel('Efficiency', fontsize=12)
            ax.set_ylim(0, 1)
            ax.grid(True)
            ax.tick_params(axis='both', which='major', labelsize=10)

        # Hide any unused subplots
        for idx in range(len(trig_vars), len(axs)):
            fig.delaxes(axs[idx])

        # Place a shared legend outside the plots on the right side
        handles, labels = axs[0].get_legend_handles_labels()
        fig.legend(handles, labels, title='Triggers', loc='upper left', bbox_to_anchor=(1.05, 1), borderaxespad=0., fontsize=12)

        # Adjust layout to accommodate legends outside of plots
        plt.tight_layout(rect=[0, 0, 0.95, 1])

        # Save and show the figure
        plt.savefig(f"{save_dir}/{trig_name}_TriggerEfficiencies_Set{chunk_index + 1}.png", dpi=200, bbox_inches='tight')
        plt.show()
        plt.close(fig)


def plot_1d_trigger_efficiencies_improvement(trig_name, prod_mode, output, triggers_set, trig_vars, baseline_key, or_trig_name, max_triggers_per_plot=1):
    """
    Plot the improvement in trigger efficiency after adding a logical OR with one of the VBF triggers, e.g.   QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1
    """
    # Create the directory to save figures
    save_dir = f"/srv/figures/{baseline_key}_baseline/{prod_mode}/1d_trigger_efficiencies_improvement"
    os.makedirs(save_dir, exist_ok=True)

    # Function to chunk the triggers_set
    def chunk_triggers(triggers, chunk_size):
        for i in range(0, len(triggers), chunk_size):
            yield triggers[i:i + chunk_size]

    # Define colors for plotting
    colors = ['blue', 'green', 'red', 'purple', 'orange', 'cyan', 'magenta', 'yellow', 'black', 'brown']

    # Chunk the triggers_set into parts of size max_triggers_per_plot
    trigger_chunks = list(chunk_triggers(triggers_set, max_triggers_per_plot))

    # Loop over each chunk of triggers
    for chunk_index, trigger_chunk in enumerate(trigger_chunks):
        num_vars = len(trig_vars)
        # Calculate the number of rows and columns
        nrows = 2
        ncols = (num_vars + 1) // 2  # Ensure all variables are covered
        fig_width = 8 * ncols  # Increase width per subplot
        fig_height = 6 * nrows  # Increase height per subplot
        fig, axs = plt.subplots(nrows, ncols, figsize=(fig_width, fig_height), squeeze=False)

        # Flatten axs array for easy indexing
        axs = axs.flatten()

        # Loop over variables to create subplots
        for var_idx, var_name in enumerate(trig_vars.keys()):
            ax = axs[var_idx]

            # Prepare bin edges for the variable
            bin_edges = trig_vars[var_name]['axis'].edges

            # Loop over triggers to plot efficiencies
            for i, trigger in enumerate(trigger_chunk):
                # Retrieve per-event data
                total_values = np.array(output[trigger][var_name]['total'])
                pass_values = np.array(output[trigger][var_name]['pass'])
                pass_or_values = np.array(output[trigger][var_name]['pass_or'])

                # Create histograms from per-event data
                total_hist, _ = np.histogram(total_values, bins=bin_edges)
                pass_hist, _ = np.histogram(pass_values, bins=bin_edges)
                pass_or_hist, _ = np.histogram(pass_or_values, bins=bin_edges)

                # Compute efficiency
                with np.errstate(divide='ignore', invalid='ignore'):
                    efficiency = np.nan_to_num(pass_hist / total_hist, nan=0.0, posinf=0.0, neginf=0.0)
                    efficiency_or = np.nan_to_num(pass_or_hist / total_hist, nan=0.0, posinf=0.0, neginf=0.0)

                # Compute bin centers for plotting
                bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

                # Plot efficiency
                ax.step(bin_centers, efficiency, where='mid', label=f'{trigger}', color=colors[i % len(colors)], linewidth=2)
                ax.step(bin_centers, efficiency_or, where='mid', label=f'{trigger} or {or_trig_name}', color=colors[(i+1) % len(colors)], linewidth=2)

            # Set axis labels and title with smaller font size
            ax.set_xlabel(trig_vars[var_name]['label'], fontsize=12)
            ax.set_ylabel('Efficiency', fontsize=12)
            ax.set_ylim(0, 1)
            ax.grid(True)
            ax.tick_params(axis='both', which='major', labelsize=10)

        # Hide any unused subplots
        for idx in range(len(trig_vars), len(axs)):
            fig.delaxes(axs[idx])

        # Place a shared legend outside the plots on the right side
        handles, labels = axs[0].get_legend_handles_labels()
        fig.legend(handles, labels, title='Triggers', loc='upper left', bbox_to_anchor=(1.05, 1), borderaxespad=0., fontsize=12)

        # Adjust layout to accommodate legends outside of plots
        plt.tight_layout(rect=[0, 0, 0.95, 1])

        # Save and show the figure
        plt.savefig(f"{save_dir}/{trigger}_or_{or_trig_name}.png", dpi=200, bbox_inches='tight')
        plt.show()
        plt.close(fig)

def plot_2d_trigger_efficiencies(trig_name, prod_mode, output, triggers_set, baseline_key, trig_vars, max_triggers_per_plot=5):
    # Create the directory to save figures
    save_dir = f"/srv/figures/{prod_mode}/2d_trigger_efficiencies"
    os.makedirs(save_dir, exist_ok=True)

    # Function to chunk the triggers_set
    def chunk_triggers(triggers, chunk_size):
        for i in range(0, len(triggers), chunk_size):
            yield triggers[i:i + chunk_size]

    # Generate all pairs of variables
    variable_names = list(trig_vars.keys())
    variable_pairs = list(itertools.combinations(variable_names, 2))

    # Chunk the triggers_set into parts of size max_triggers_per_plot
    trigger_chunks = list(chunk_triggers(triggers_set, max_triggers_per_plot))

    # Loop over each chunk of triggers
    for chunk_index, trigger_chunk in enumerate(trigger_chunks):
        # Loop over variable pairs
        for var_pair in variable_pairs:
            var_x, var_y = var_pair
            num_triggers = len(trigger_chunk)
            nrows = 2
            ncols = (num_triggers + 1) // 2  # Ensure all triggers are covered
            fig_width = 8 * ncols
            fig_height = 6 * nrows
            fig, axs = plt.subplots(nrows, ncols, figsize=(fig_width, fig_height), squeeze=False)
            axs = axs.flatten()

            # Loop over triggers in the chunk
            for i, trigger in enumerate(trigger_chunk):
                ax = axs[i]

                # Get per-event variable arrays
                var_x_total = output[trigger][var_x]['total']
                var_x_pass = output[trigger][var_x]['pass']
                var_y_total = output[trigger][var_y]['total']
                var_y_pass = output[trigger][var_y]['pass']

                # Create 2D histograms for total and passing events
                H_total, xedges, yedges = np.histogram2d(
                    var_x_total, var_y_total,
                    bins=[trig_vars[var_x]['axis'].edges, trig_vars[var_y]['axis'].edges]
                )
                H_pass, _, _ = np.histogram2d(
                    var_x_pass, var_y_pass,
                    bins=[trig_vars[var_x]['axis'].edges, trig_vars[var_y]['axis'].edges]
                )

                # Calculate efficiency
                with np.errstate(divide='ignore', invalid='ignore'):
                    efficiency = np.nan_to_num(H_pass / H_total, nan=0.0, posinf=0.0, neginf=0.0)

                # Plot 2D efficiency
                X, Y = np.meshgrid(trig_vars[var_x]['axis'].edges, trig_vars[var_y]['axis'].edges)
                pcm = ax.pcolormesh(X, Y, efficiency.T, cmap='viridis', shading='auto')
                fig.colorbar(pcm, ax=ax, label='Efficiency')

                # Set axis labels with smaller font size
                ax.set_xlabel(trig_vars[var_x]['label'], fontsize=12)
                ax.set_ylabel(trig_vars[var_y]['label'], fontsize=12)
                ax.set_title(f'{trigger}', fontsize=14)
                ax.tick_params(axis='both', which='major', labelsize=10)

            # Hide any unused subplots
            for idx in range(len(trigger_chunk), len(axs)):
                fig.delaxes(axs[idx])

            # Adjust layout to accommodate larger subplots
            plt.tight_layout(rect=[0, 0, 1, 0.95])

            # Set a super title for the figure
            fig.suptitle(f'2D Efficiency: {var_x} vs {var_y}', fontsize=16)

            # Save and show the figure
            plt.savefig(f"{save_dir}/{var_x}_vs_{var_y}_Set{chunk_index + 1}.png", dpi=200, bbox_inches='tight')
            plt.show()
            plt.close(fig)

def plot_fancy_2d_trigger_efficiencies(trig_name, prod_mode, output, triggers_set, trig_vars, baseline_key, max_triggers_per_plot=5):
    # Create the directory to save figures
    save_dir = f"/srv/figures/{baseline_key}_baseline/{prod_mode}/fancy_2d_trigger_efficiencies"
    os.makedirs(save_dir, exist_ok=True)

    # Function to chunk the triggers_set
    def chunk_triggers(triggers, chunk_size):
        for i in range(0, len(triggers), chunk_size):
            yield triggers[i:i + chunk_size]

    # Generate all pairs of variables
    variable_names = list(trig_vars.keys())
    variable_pairs = list(itertools.combinations(variable_names, 2))

    # Chunk the triggers_set into parts of size max_triggers_per_plot
    trigger_chunks = list(chunk_triggers(triggers_set, max_triggers_per_plot))

    # Loop over each chunk of triggers
    for chunk_index, trigger_chunk in enumerate(trigger_chunks):
        # Loop over variable pairs
        for var_pair in variable_pairs:
            var_x, var_y = var_pair

            # Loop over triggers in the chunk
            for trigger in trigger_chunk:
                # Get per-event variable arrays
                var_x_total = np.array(output[trigger][var_x]['total'])
                var_x_pass = np.array(output[trigger][var_x]['pass'])
                var_y_total = np.array(output[trigger][var_y]['total'])
                var_y_pass = np.array(output[trigger][var_y]['pass'])

                # Create 2D histograms for total and passing events
                H_total, xedges, yedges = np.histogram2d(
                    var_x_total, var_y_total,
                    bins=[trig_vars[var_x]['axis'].edges, trig_vars[var_y]['axis'].edges]
                )
                H_pass, _, _ = np.histogram2d(
                    var_x_pass, var_y_pass,
                    bins=[trig_vars[var_x]['axis'].edges, trig_vars[var_y]['axis'].edges]
                )

                # Calculate efficiency
                with np.errstate(divide='ignore', invalid='ignore'):
                    efficiency = np.nan_to_num(H_pass / H_total, nan=0.0, posinf=0.0, neginf=0.0)

                # Prepare data for plotting
                x_centers = (xedges[:-1] + xedges[1:]) / 2
                y_centers = (yedges[:-1] + yedges[1:]) / 2
                X, Y = np.meshgrid(x_centers, y_centers)
                Z = efficiency.T  # Transpose to match X and Y

                # Create marginal histograms
                # For total events
                total_hist_x, _ = np.histogram(var_x_total, bins=xedges)
                total_hist_y, _ = np.histogram(var_y_total, bins=yedges)
                # For passing events
                pass_hist_x, _ = np.histogram(var_x_pass, bins=xedges)
                pass_hist_y, _ = np.histogram(var_y_pass, bins=yedges)
                # Compute efficiencies
                with np.errstate(divide='ignore', invalid='ignore'):
                    efficiency_x = np.nan_to_num(pass_hist_x / total_hist_x, nan=0.0, posinf=0.0, neginf=0.0)
                    efficiency_y = np.nan_to_num(pass_hist_y / total_hist_y, nan=0.0, posinf=0.0, neginf=0.0)
                # Compute bin centers
                bin_centers_x = x_centers
                bin_centers_y = y_centers

                # Create figure with subplots for marginal histograms
                fig = make_subplots(
                    rows=2, cols=2,
                    shared_xaxes=True,
                    shared_yaxes=True,
                    row_heights=[0.2, 0.8],
                    column_widths=[0.8, 0.2],
                    specs=[[{"type": "xy"}, {"type": "xy"}],
                           [{"type": "xy"}, {"type": "xy"}]],
                    horizontal_spacing=0.02,
                    vertical_spacing=0.02
                )

                # Add heatmap to the bottom-left cell (row=2, col=1)
                heatmap = go.Heatmap(
                    x=x_centers,
                    y=y_centers,
                    z=Z,
                    colorscale='Viridis',
                    colorbar=dict(title='Efficiency'),
                    showscale=True
                )
                fig.add_trace(heatmap, row=2, col=1)

                # Add x efficiency histogram to the top-left cell (row=1, col=1)
                hist_x = go.Bar(
                    x=bin_centers_x,
                    y=efficiency_x,
                    marker=dict(color='blue'),
                    showlegend=False
                )
                fig.add_trace(hist_x, row=1, col=1)

                # Add y efficiency histogram to the bottom-right cell (row=2, col=2)
                hist_y = go.Bar(
                    x=efficiency_y,
                    y=bin_centers_y,
                    orientation='h',
                    marker=dict(color='blue'),
                    showlegend=False
                )
                fig.add_trace(hist_y, row=2, col=2)

                # Hide empty subplot (top-right cell)
                fig.update_xaxes(visible=False, row=1, col=2)
                fig.update_yaxes(visible=False, row=1, col=2)

                # Update axes labels
                fig.update_xaxes(title_text='', row=1, col=1)
                fig.update_xaxes(title_text=trig_vars[var_x]['label'], row=2, col=1)
                fig.update_yaxes(title_text='Efficiency', row=1, col=1)
                fig.update_yaxes(title_text=trig_vars[var_y]['label'], row=2, col=1)
                fig.update_xaxes(title_text='Efficiency', row=2, col=2)
                fig.update_yaxes(title_text='', row=2, col=2)

                # Adjust layout
                fig.update_layout(
                    title_text=f'{var_x} vs {var_y} for {trigger}',
                    width=800,
                    height=800,
                    showlegend=False,
                )

                # Save figure as an HTML file
                new_save_dir = save_dir + f"/{var_x}_vs_{var_y}"
                os.makedirs(new_save_dir, exist_ok=True)
                fig.write_html(f"{new_save_dir}/{trigger}.html")
                # Save figure as a JPG file
                fig.write_image(f"{new_save_dir}/{trigger}.jpg")

def compare_trigger_efficiencies(outputs, triggers_dict, trig_vars, baseline_key):
    # Create a directory to save figures
    save_dir = f"./figures/{baseline_key}_baseline/trigger_efficiency_comparisons"
    os.makedirs(save_dir, exist_ok=True)
    
    # Define colors for plotting
    colors = {'ggF': 'blue', 'VBF': 'red'}
    
    # Loop over trigger categories
    for name, triggers_list in triggers_dict.items():
        # Loop over triggers in each category
        for trigger in triggers_list:
            # Loop over variables to create efficiency plots
            for var_name in trig_vars.keys():
                plt.figure(figsize=(10, 6))
                
                # Prepare bin edges for the variable
                bin_edges = trig_vars[var_name]['axis'].edges
                
                # Plot efficiencies for each production mode
                for prod_mode in outputs.keys():  # 'ggF' and 'VBF'
                    # Check if outputs[prod_mode][name] exists
                    if name in outputs[prod_mode]:
                        output = outputs[prod_mode][name]
                        # Check if the trigger exists in output
                        if trigger in output:
                            # Get the per-event data arrays
                            total_values = np.array(output[trigger][var_name]['total'])
                            pass_values = np.array(output[trigger][var_name]['pass'])
                            
                            # Create histograms from per-event data
                            total_hist, _ = np.histogram(total_values, bins=bin_edges)
                            pass_hist, _ = np.histogram(pass_values, bins=bin_edges)
                            
                            # Compute efficiency
                            with np.errstate(divide='ignore', invalid='ignore'):
                                efficiency = np.nan_to_num(pass_hist / total_hist, nan=0.0, posinf=0.0, neginf=0.0)
                            
                            # Compute bin centers
                            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
                            
                            # Plot efficiency
                            plt.step(
                                bin_centers,
                                efficiency,
                                where='mid',
                                label=f'{prod_mode}',
                                color=colors[prod_mode],
                                linewidth=2
                            )
                        else:
                            print(f"Trigger '{trigger}' not found in outputs for production mode '{prod_mode}' and category '{name}'.")
                    else:
                        print(f"Category '{name}' not found in outputs for production mode '{prod_mode}'.")
                
                # Set axis labels and title
                plt.xlabel(trig_vars[var_name]['label'], fontsize=12)
                plt.ylabel('Efficiency', fontsize=12)
                plt.title(f'Trigger Efficiency Comparison for {trigger} - {trig_vars[var_name]["label"]}', fontsize=14)
                plt.legend()
                plt.grid(True)
                plt.ylim(0, 1)
                
                # Save and show the figure
                filename = f"{trigger}_{var_name}_efficiency_comparison.png".replace(' ', '_').replace('/', '_')
                plt.tight_layout()
                plt.savefig(f"{save_dir}/{filename}", dpi=200)
                plt.show()
                plt.close()
def plot_cutflow(output, save_dir, small_font_size=12):
    """
    Creates and saves a cutflow diagram from the provided output dictionary,
    displaying the number of events passing each selection step on top of the bars.

    Parameters:
    output (list): A list containing dictionaries, where the first dictionary should have a 'cutflow' key.
    save_dir (str): The directory where the plot will be saved.
    small_font_size (int, optional): Font size for labels and title. Default is 12.
    """
    # Access the cutflow
    cutflow = output[0]['cutflow']

    # Print the cutflow
    for cut, count in cutflow.items():
        print(f"{cut}: {count}")

    # Prepare data
    cuts = list(cutflow.keys())
    counts = [cutflow[cut] for cut in cuts]

    # Plot
    plt.figure(figsize=(10, 6))
    plt.yscale("log")
    bars = plt.bar(cuts, counts, color='skyblue')

    # Add counts on top of each bar using bar_label (Matplotlib ≥ 3.4)
    plt.bar_label(bars, labels=[str(count) for count in counts], padding=3, fontsize=small_font_size)

    # Set smaller font sizes for labels and title
    plt.xlabel('Selection Steps', fontsize=small_font_size)
    plt.ylabel('Number of Events', fontsize=small_font_size)
    plt.title('Cutflow Diagram', fontsize=small_font_size)

    # Set smaller font sizes for tick labels
    plt.xticks(rotation=45, fontsize=small_font_size)
    plt.yticks(fontsize=small_font_size)

    plt.grid(axis='y')
    plt.tight_layout()

    # Ensure save directory exists
    os.makedirs(save_dir, exist_ok=True)

    # Save the plot
    plt.savefig(os.path.join(save_dir, "cutFlow.png"), bbox_inches='tight', dpi=300)
    plt.show()
    plt.close()

# Example usage:
# plot_cutflow(output, "/srv/figures/VBF_baseline/cutFlow/")

def plot_1d_trigger_soup(output, trig_vars, tags=[], save_dir=None):
    """
    Plot the trigger efficiencies for each cumulative trigger combination.

    Parameters:
    - output: dict
        The output dictionary from the TriggerEfficiencyImprovementProcessor.
    - trig_vars: dict
        Dictionary of trigger variables with their processing functions and axis information.
    - save_dir: str, optional
        Directory to save the figures. If None, defaults to "./figures".
    """
    if save_dir is None:
        save_dir = "./figures"
    os.makedirs(save_dir, exist_ok=True)

    # Define colors for plotting
    colors = ['blue', 'green', 'red', 'purple', 'orange', 'cyan', 'magenta', 'yellow', 'black', 'brown']

    # Extract the combinations from the output
    combinations = list(output.keys())

    num_vars = len(trig_vars)
    ncols = 2
    nrows = (num_vars + 1) // 2
    fig_width = 8 * ncols
    fig_height = 6 * nrows
    fig, axs = plt.subplots(nrows, ncols, figsize=(fig_width, fig_height), squeeze=False)
    axs = axs.flatten()

    for var_idx, var_name in enumerate(trig_vars.keys()):
        ax = axs[var_idx]
        # Prepare bin edges for the variable
        bin_edges = trig_vars[var_name]['axis'].edges

        for i, combination in enumerate(combinations):
            # Retrieve per-event data
            total_values = np.array(output[combination][var_name]['total'])
            pass_values = np.array(output[combination][var_name]['pass'])

            # Create histograms from per-event data
            total_hist, _ = np.histogram(total_values, bins=bin_edges)
            pass_hist, _ = np.histogram(pass_values, bins=bin_edges)

            # Compute efficiency
            with np.errstate(divide='ignore', invalid='ignore'):
                efficiency = np.nan_to_num(pass_hist / total_hist, nan=0.0, posinf=0.0, neginf=0.0)

            # Compute bin centers for plotting
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

            # Plot efficiency
            ax.step(
                bin_centers,
                efficiency,
                where='mid',
                label=f'{combination}',
                color=colors[i % len(colors)],
                linewidth=2
            )

        # Set axis labels and title
        ax.set_xlabel(trig_vars[var_name]['label'], fontsize=12)
        ax.set_ylabel('Efficiency', fontsize=12)
        ax.set_ylim(0, 1)
        ax.grid(True)
        ax.tick_params(axis='both', which='major', labelsize=10)

    # Hide any unused subplots
    for idx in range(len(trig_vars), len(axs)):
        fig.delaxes(axs[idx])

    # Place a shared legend outside the plots on the right side
    handles, labels = axs[0].get_legend_handles_labels()
    fig.legend(
        handles,
        tags,
        title='Trigger Combinations',
        loc='upper left',
        bbox_to_anchor=(1.05, 1),
        borderaxespad=0.,
        fontsize=12
    )

    # Adjust layout to accommodate legends outside of plots
    plt.tight_layout(rect=[0, 0, 0.95, 1])

    # Save and show the figure
    plt.savefig(f"{save_dir}/trigger_soup.png", dpi=200, bbox_inches='tight')
    plt.show()
    plt.close(fig)

