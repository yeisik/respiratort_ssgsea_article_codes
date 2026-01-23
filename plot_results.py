import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import glob

# Configuration
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_PP_DIR = os.path.join(BASE_DIR, "Results_PP")
OUTPUT_PLOT_DIR = os.path.join(BASE_DIR, "Plots")
os.makedirs(OUTPUT_PLOT_DIR, exist_ok=True)

def load_summaries():
    """Loads and concatenates all summary CSV files."""
    pattern = os.path.join(RESULTS_PP_DIR, "Summary_*.csv")
    files = glob.glob(pattern)
    
    if not files:
        print("No summary files found in", RESULTS_PP_DIR)
        return None
        
    dfs = []
    for f in files:
        try:
            df = pd.read_csv(f)
            dfs.append(df)
        except Exception as e:
            print(f"Error reading {f}: {e}")
            
    if not dfs:
        return None
        
    combined_df = pd.concat(dfs, ignore_index=True)
    return combined_df

def plot_performance(df, metric="AUPRC", target="SC1"):
    """
    Generates a faceted bar plot mimicking the user's example.
    Filters: Only 'Combined' and 'Probe'.
    Style: Stacked/Overlay bars (Combined behind Probe).
    """
    subset = df[df["Target"] == target].copy()
    if subset.empty:
        print(f"No data for target {target}")
        return

    # Filter Features
    desired_features = ["Combined", "Probe"]
    subset = subset[subset["Feature"].isin(desired_features)]
    
    if subset.empty:
        print(f"No data for features {desired_features} in target {target}")
        return

    subset = subset.sort_values("TimePoint")

    # Set style
    sns.set_style("ticks")
    sns.set_context("talk")
    
    # Custom Palette: Combined (Dark Blue), Probe (Light Blue)
    # We map specifically to ensure consistency
    palette = {"Combined": "#006FA6", "Probe": "#89CFF0"}  # Example hex codes matching descriptions

    # We want to plot Combined first (so it's in the background if it's larger), 
    # then Probe. But standard seaborn hue order might vary.
    # Logic: If we use dodge=False, they are plotted on top of each other.
    # We need to sort subset so that the Larger value is likely plotted first?
    # Actually, seaborn plots hue levels in order.
    # If we set hue_order=["Combined", "Probe"], Combined is plotted, then Probe.
    # If Combined > Probe (usually true), Probe will be visible on top.
    # If Probe > Combined, Probe covers Combined.
    # This matches the "Stacked" visual where the top is the max value.

    # FacetGrid
    g = sns.FacetGrid(
        data=subset,
        col="TimePoint",
        height=5,
        aspect=0.4,
        sharey=True
    )

    # Map Barplot
    # We wrap sns.barplot to pass specific args
    def overlaid_barplot(*args, **kwargs):
        data = kwargs.pop('data')
        # Ensure Combined is plotted first (bottom layer), Probe second (top layer)
        # Using hue_order to control plotting order
        sns.barplot(
            data=data,
            x="Classifier",
            y=metric,
            hue="Feature",
            hue_order=["Combined", "Probe"], # Plot Combined background, Probe foreground
            palette=palette,
            dodge=False, # Overlay
            edgecolor="white",
            linewidth=1,
            **kwargs
        )

    g.map_dataframe(overlaid_barplot)

    # Adjust layout
    g.fig.subplots_adjust(top=0.9, wspace=0.1)
    
    # Titles & Axes
    g.set_titles("{col_name}") 
    g.set_axis_labels("", metric)
    
    # Ticks & Grid
    for ax in g.axes.flat:
        ax.tick_params(axis='x', rotation=90)
        ax.grid(False, axis='x')
        ax.grid(True, axis='y', linestyle='--', alpha=0.5)

    # Legend
    # We need to manually add legend because map_dataframe with hue inside doesn't automatically add it to Grid 
    # properly when using custom wrapper sometimes, or it duplicates. 
    # Let's handle legend manually to ensure correct handles.
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=palette["Combined"], edgecolor='white', label='Combined'),
        Patch(facecolor=palette["Probe"], edgecolor='white', label='Probe')
    ]
    g.fig.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1, 1), title="Feature", frameon=False)

    output_path = os.path.join(OUTPUT_PLOT_DIR, f"BarPlot_{target}_{metric}.png")
    g.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved plot to {output_path}")

def main():
    df = load_summaries()
    if df is not None:
        targets = df["Target"].unique()
        for t in targets:
            plot_performance(df, metric="AUPRC", target=t)

if __name__ == "__main__":
    main()
