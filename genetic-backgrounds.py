"""
Simulates how polygenic score (PGS), log odds ratios, and genetic background affect
disease prevalence and risk. Generates a 4-panel figure illustrating these relationships.
"""
import matplotlib.pyplot as plt 
import numpy as np

def main():
    BETA = 1
    PGS_RANGE = (-4.4, 4.4)
    N_PGS_POINTS = 500

    LOW_RISK_PGS = -2
    LOW_RISK_VARIANT_PGS = -1
    HIGH_RISK_PGS = 3
    HIGH_RISK_VARIANT_PGS = 4
    
    # Define the range of PGS to plot
    pgs = np.linspace(PGS_RANGE[0], PGS_RANGE[1], N_PGS_POINTS)

    # Get data for logOR and OR
    log_odds_ratio = BETA * pgs  # Log odds ratio
    odds_ratio = np.exp2(log_odds_ratio)  # Odds ratio

    # Plotting parameters
    SMALL_SIZE = 12
    BIGGER_SIZE = 16

    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    fig, axs = plt.subplots(2,2, figsize=(10,8))

    # Figure labels
    plt.gcf().text(0.01, 1.0, 'A)', fontsize=14, fontweight='bold', va='top', ha='left')
    plt.gcf().text(0.51, 1.0, 'B)', fontsize=14, fontweight='bold', va='top', ha='left')
    plt.gcf().text(0.01, 0.5, 'C)', fontsize=14, fontweight='bold', va='top', ha='left')
    plt.gcf().text(0.51, 0.5, 'D)', fontsize=14, fontweight='bold', va='top', ha='left')
    
    # Colours and shading
    colors = {'low': 'green', 'high': 'red'}
    alphas = {'shade': 0.1,'fill': 0.4}

    # --- Plot A: log odds ratio vs PGS --- 
    plt.subplot(2, 2, 1)
    plt.xlim(-4.4, 4.4)
    plt.ylim(-4.4, 4.4)
    plt.plot(pgs, log_odds_ratio, color="black")
    plt.scatter([LOW_RISK_VARIANT_PGS], BETA*LOW_RISK_VARIANT_PGS, color = colors['low'])
    plt.scatter([LOW_RISK_PGS], [BETA*LOW_RISK_PGS], color = colors['low'])

    # delta-variant in high risk background 
    plt.scatter([HIGH_RISK_VARIANT_PGS], BETA*HIGH_RISK_VARIANT_PGS, color = colors['high'])
    plt.scatter([HIGH_RISK_PGS], [BETA*HIGH_RISK_PGS], color = colors['high'])

    # Create the y-range for the shaded area
    y_range = np.linspace(BETA*LOW_RISK_PGS, BETA*LOW_RISK_VARIANT_PGS, N_PGS_POINTS)
    x_right_bound = (y_range / BETA)

    # Fill the shaded area with four bounds
    plt.fill_betweenx(
        y_range,  # Vertical range
        plt.xlim()[0],  # Left edge of the plot
        x_right_bound,  # Right edge (odds_ratio curve)
        color=colors['low'], alpha=alphas['fill'], label="Low PGS background",
        edgecolor = None
    )
    plt.fill_between(y_range, PGS_RANGE[0], x_right_bound, color=colors['low'], alpha=alphas['shade'], edgecolor=None)

    y_range = np.linspace(BETA*HIGH_RISK_PGS, BETA*HIGH_RISK_VARIANT_PGS, N_PGS_POINTS)
    x_right_bound = (y_range) / BETA

    plt.fill_betweenx(
        y_range,  # Vertical range
        plt.xlim()[0],  # Left edge of the plot
        x_right_bound,  # Right edge (odds_ratio curve)
        color=colors['high'], alpha=alphas['fill'], label="High PGS background",
        edgecolor = None
    )

    plt.fill_between(y_range, PGS_RANGE[0], x_right_bound, color=colors['high'], alpha=alphas['shade'], edgecolor=None)

    plt.xlabel("Polygenic Score (PGS)")
    plt.ylabel("Log Odds Ratio")

    # --- Plot B: odds ratio ---
    plt.subplot(2, 2, 2)
    plt.plot(pgs, odds_ratio, color="black")
    plt.xlim(-4.4, 4.4)
    plt.ylim(bottom=-0.4, top=16.7)

    # delta-variant in low-risk background
    plt.scatter([LOW_RISK_VARIANT_PGS], np.exp2(BETA*LOW_RISK_VARIANT_PGS), color = 'green')
    plt.scatter([LOW_RISK_PGS], [np.exp2(BETA*LOW_RISK_PGS)], color = 'green')

    # delta-variant in high-risk background 
    plt.scatter([HIGH_RISK_VARIANT_PGS], np.exp2(BETA*HIGH_RISK_VARIANT_PGS), color = 'red')
    plt.scatter([HIGH_RISK_PGS], [np.exp2(BETA*HIGH_RISK_PGS)], color = 'red')

    # Create the y-range for the shaded area
    y_range = np.linspace(np.exp2(BETA*LOW_RISK_PGS), np.exp2(BETA*LOW_RISK_VARIANT_PGS), 500)
    x_right_bound = (np.log2(y_range) / BETA)

    # Fill the shaded area with four bounds
    plt.fill_betweenx(
        y_range,  # Vertical range
        plt.xlim()[0],  # Left edge of the plot
        x_right_bound,  # Right edge (odds_ratio curve)
        color=colors['low'], alpha=alphas['fill'], label="Low PGS background",
        edgecolor = None
    )
    plt.fill_between(y_range, -0.4, x_right_bound, color=colors['low'], alpha=alphas['shade'], edgecolor=None)
    plt.fill_between(x_right_bound, -0.4, y_range, color=colors['low'], alpha=alphas['shade'], edgecolor=None)


    y_range = np.linspace(np.exp2(BETA*HIGH_RISK_PGS), np.exp2(BETA*HIGH_RISK_VARIANT_PGS), N_PGS_POINTS)
    x_right_bound = (np.log2(y_range) / BETA)

    plt.fill_betweenx(
        y_range,  # Vertical range
        plt.xlim()[0],  # Left edge of the plot
        x_right_bound,  # Right edge (odds_ratio curve)
        color=colors['high'], alpha=alphas['fill'], label="High PGS background",
        edgecolor = None
    )
    plt.fill_between(y_range, -0.4, x_right_bound, color=colors['high'], alpha=alphas['shade'], edgecolor=None)
    plt.fill_between(x_right_bound, -0.4, y_range, color=colors['high'], alpha=alphas['shade'], edgecolor=None)

    plt.xlabel("Polygenic Score (PGS)")
    plt.ylabel("Odds Ratio")

    # --- Plot 3: Differing background disease rates ---
    ax3 = plt.subplot(2,2,3)
    plt.xlim(-4.4, 4.4)
    plt.ylabel("Disease prevalence (%)")
    plt.xlabel("Polygenic Score (PGS)")

    # Disease prevalence in reference group 
    prevalences = [0.001, 0.005, 0.01, 0.05, 0.1]

    for prevalence in prevalences:
        odds = prevalence / (1 - prevalence)  # Calculate initial odds
        new_prevalence = (odds * odds_ratio) / (1 + (odds * odds_ratio))  # Adjust odds, convert back to prevalence
        plt.plot(pgs, new_prevalence * 100, label=f"{prevalence*100:.1f}%")  # Convert to percentage
        
    plt.legend(title="Disease prevalence in reference group", ncols=2)

    # Shade high- and low-risk ranges
    y_range_bottom, y_range_top = plt.ylim()
    plt.ylim(plt.ylim())
    y = np.arange(y_range_bottom, y_range_top, 0.1)
    plt.fill_betweenx(y, LOW_RISK_PGS, LOW_RISK_VARIANT_PGS, 
                    color=colors['low'], alpha=alphas['shade'], edgecolor = None)

    plt.fill_betweenx(y, HIGH_RISK_PGS, HIGH_RISK_VARIANT_PGS, 
                    color=colors['high'], alpha=alphas['shade'], edgecolor = None)

    # --- Plot 4: Barplots of background-dependent change in disease prevalence ---
    ax = plt.subplot(2,2,4)
    high_risk_changes = []
    low_risk_changes = [] 
    for prevalence in prevalences:
        odds = prevalence / (1 - prevalence)
        
        # Calculate low-risk changes
        low_risk_changes.append(
            (np.exp2(BETA * LOW_RISK_VARIANT_PGS) * odds /
            (1 + odds * np.exp2(BETA * LOW_RISK_VARIANT_PGS))) * 100 -
            (np.exp2(BETA * LOW_RISK_PGS) * odds /
            (1 + odds * np.exp2(BETA * LOW_RISK_VARIANT_PGS))) * 100
        )
        
        # Calculate high-risk changes
        high_risk_changes.append(
            (np.exp2(BETA * HIGH_RISK_VARIANT_PGS) * odds /
            (1 + odds * np.exp2(BETA * HIGH_RISK_VARIANT_PGS))) * 100 -
            (np.exp2(BETA * HIGH_RISK_PGS) * odds /
            (1 + odds * np.exp2(BETA * HIGH_RISK_PGS))) * 100
        )


    x = np.arange(len(prevalences))  # x positions for bars
    width = 0.4  # Width of  bars

    bars_low = ax.bar(x-width/2, low_risk_changes, width, color = colors['low'], alpha=alphas['fill'])
    bars_high = ax.bar(x+width/2, high_risk_changes, width, color=colors['high'], alpha=alphas['fill'])
    plt.xticks(x)

    for bar in bars_low:
        ax.text(
            bar.get_x() + bar.get_width() / 2,  # Center of the bar
            bar.get_height(),                  # Height of the bar
            f"{bar.get_height():.2f}",         # Label text
            ha='center', va='bottom', fontsize=10
        )
    for bar in bars_high:
        ax.text(
            bar.get_x() + bar.get_width() / 2, 
            bar.get_height(), 
            f"{bar.get_height():.2f}", 
            ha='center', va='bottom', fontsize=10
        )
        
    ax.set_xticks(x)
    ax.set_ylim(top=max(high_risk_changes)+2)
    ax.set_xticklabels([f"{p:.1%}" for p in prevalences])
    plt.xlabel("Disease prevalence in reference group")
    plt.ylabel("Change in \ndisease prevalence (%)")

    # Custom legend:
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=colors['low'], alpha = alphas['fill'], label = "Low PGS background"),
                    Patch(facecolor=colors['high'], alpha=alphas['fill'], label = "High PGS background")]
    fig.legend(handles=legend_elements, 
            loc = 'upper center', bbox_transform=fig.transFigure,
            title="Variant Effect", ncols=2,
            bbox_to_anchor = (0.5, 1.1))

    plt.tight_layout()
    plt.show()
    plt.savefig("genetic-backgrounds.png", dpi=300, bbox_inches="tight")

if __name__ == "__main__":
    main()