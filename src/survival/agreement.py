import numpy as np
from scipy.stats import binomtest, norm
import matplotlib.pyplot as plt
import pandas as pd
import pingouin as pg

def paired_permutation_test(x, y, n_permutations=10000, alternative="two-sided", random_state=None):
    """
    Paired permutation test using sign flipping.

    Parameters
    ----------
    x, y : array-like
        Paired samples (same length)
    n_permutations : int
        Number of permutations
    alternative : {"two-sided", "greater", "less"}
    random_state : int or None

    Returns
    -------
    stat : float
        Observed test statistic (mean difference)
    p_value : float
    """
    rng = np.random.default_rng(random_state)

    x = np.asarray(x)
    y = np.asarray(y)
    d = x - y

    # Observed statistic
    stat_obs = np.mean(d)

    # Generate permutation distribution via sign flipping
    signs = rng.choice([-1, 1], size=(n_permutations, len(d)))
    perm_stats = np.mean(signs * d, axis=1)

    if alternative == "two-sided":
        p_value = np.mean(np.abs(perm_stats) >= abs(stat_obs))
    elif alternative == "greater":
        p_value = np.mean(perm_stats >= stat_obs)
    elif alternative == "less":
        p_value = np.mean(perm_stats <= stat_obs)
    else:
        raise ValueError("Invalid alternative")

    return stat_obs, p_value

def sign_test(x, y, alternative="two-sided"):
    """
    Sign test for paired samples.

    Parameters
    ----------
    x, y : array-like
        Paired samples
    alternative : {"two-sided", "greater", "less"}

    Returns
    -------
    n_pos : int
        Number of positive differences
    n : int
        Number of non-zero differences
    p_value : float
    """
    x = np.asarray(x)
    y = np.asarray(y)

    d = x - y
    d = d[d != 0]  # remove ties

    n = len(d)
    n_pos = np.sum(d > 0)

    result = binomtest(n_pos, n, p=0.5, alternative=alternative)

    return n_pos, n, result.pvalue

def paired_permutation_test_median(x, y, n_permutations=20000, 
                                  alternative="two-sided", random_state=None):
    """
    Paired permutation test using the median difference.
    """
    rng = np.random.default_rng(random_state)

    x = np.asarray(x)
    y = np.asarray(y)
    d = x - y

    stat_obs = np.median(d)

    signs = rng.choice([-1, 1], size=(n_permutations, len(d)))
    perm_stats = np.median(signs * d, axis=1)

    if alternative == "two-sided":
        p_value = np.mean(np.abs(perm_stats) >= abs(stat_obs))
    elif alternative == "greater":
        p_value = np.mean(perm_stats >= stat_obs)
    elif alternative == "less":
        p_value = np.mean(perm_stats <= stat_obs)
    else:
        raise ValueError("Invalid alternative")

    return stat_obs, p_value

def paired_bootstrap_ci(x, y, statistic="median", 
                        n_bootstrap=20000, ci=0.95, random_state=None):
    """
    Bootstrap CI for paired differences.
    """
    rng = np.random.default_rng(random_state)

    x = np.asarray(x)
    y = np.asarray(y)
    d = x - y
    n = len(d)

    stats = []

    for _ in range(n_bootstrap):
        idx = rng.integers(0, n, n)
        sample = d[idx]

        if statistic == "median":
            stats.append(np.median(sample))
        elif statistic == "mean":
            stats.append(np.mean(sample))
        else:
            raise ValueError("Invalid statistic")

    alpha = 1 - ci
    lower = np.percentile(stats, 100 * alpha / 2)
    upper = np.percentile(stats, 100 * (1 - alpha / 2))

    return lower, upper

def paired_bootstrap_ci_bca(x, y, statistic="median", 
                           n_bootstrap=20000, ci=0.95, random_state=None):
    rng = np.random.default_rng(random_state)

    x = np.asarray(x)
    y = np.asarray(y)
    d = x - y
    n = len(d)

    # statistic function
    if statistic == "median":
        stat_fn = np.median
    elif statistic == "mean":
        stat_fn = np.mean
    else:
        raise ValueError("Invalid statistic")

    # observed statistic
    theta_hat = stat_fn(d)

    # bootstrap samples
    boot_stats = np.empty(n_bootstrap)
    for i in range(n_bootstrap):
        idx = rng.integers(0, n, n)
        boot_stats[i] = stat_fn(d[idx])

    # bias correction (z0)
    z0 = norm.ppf(np.mean(boot_stats < theta_hat))

    # jackknife for acceleration (a)
    jack_stats = np.empty(n)
    for i in range(n):
        jack_sample = np.delete(d, i)
        jack_stats[i] = stat_fn(jack_sample)

    jack_mean = np.mean(jack_stats)
    num = np.sum((jack_mean - jack_stats)**3)
    den = 6 * (np.sum((jack_mean - jack_stats)**2)**1.5)
    a = num / den if den != 0 else 0.0

    # adjusted quantiles
    alpha = 1 - ci
    z_low = norm.ppf(alpha / 2)
    z_high = norm.ppf(1 - alpha / 2)

    def adjusted_quantile(z):
        return norm.cdf(z0 + (z0 + z) / (1 - a * (z0 + z)))

    q_low = adjusted_quantile(z_low)
    q_high = adjusted_quantile(z_high)

    lower = np.quantile(boot_stats, q_low)
    upper = np.quantile(boot_stats, q_high)

    return lower, upper

def bland_altman_analysis(method1, method2, plots_path, name1="pPCI", name2="sPCI"):
    mean     = (method1 + method2) / 2
    diff     = method1 - method2
    mean_diff = np.mean(diff)
    std_diff  = np.std(diff, ddof=1)
    loa_upper = mean_diff + 1.96 * std_diff
    loa_lower = mean_diff - 1.96 * std_diff

    # --- Stats ---
    print("=== Bland-Altman Results ===")
    print(f"Mean bias ({name1} - {name2}): {mean_diff:.3f}")
    print(f"SD of differences:             {std_diff:.3f}")
    print(f"Upper limit of agreement:      {loa_upper:.3f}")
    print(f"Lower limit of agreement:      {loa_lower:.3f}")

    # --- Plot ---
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(mean, diff, color="steelblue", edgecolors="white", s=70, zorder=3)

    ax.axhline(mean_diff, color="black",  linestyle="-",  linewidth=1.5, label=f"Mean bias: {mean_diff:.2f}")
    ax.axhline(loa_upper, color="crimson", linestyle="--", linewidth=1.2, label=f"+1.96 SD: {loa_upper:.2f}")
    ax.axhline(loa_lower, color="crimson", linestyle="--", linewidth=1.2, label=f"-1.96 SD: {loa_lower:.2f}")
    ax.axhline(0,         color="gray",   linestyle=":",  linewidth=1.0, label="Zero line")

    ax.fill_between(ax.get_xlim(), loa_lower, loa_upper, alpha=0.05, color="crimson")

    ax.set_xlabel(f"Mean of {name1} and {name2}", fontsize=12)
    ax.set_ylabel(f"Difference ({name1} − {name2})", fontsize=12)
    ax.set_title("Bland-Altman Plot", fontsize=14)
    ax.legend(fontsize=10)
    plt.tight_layout()
    plt.savefig(f"{plots_path}/bland_altman.png", dpi=600)

    # --- Scatter plot with line of identity ---
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(method1, method2, color="steelblue", edgecolors="white", s=70, zorder=3)

    lims = (min(np.min(method1), np.min(method2)) - 1,
            max(np.max(method1), np.max(method2)) + 1)
    ax.plot(lims, lims, "k--", linewidth=1.2, label="Line of identity (y = x)")

    ax.set_xlabel(name1, fontsize=12)
    ax.set_ylabel(name2, fontsize=12)
    ax.set_title(f"{name1} vs {name2}", fontsize=14)
    ax.set_xlim(lims); ax.set_ylim(lims)
    ax.legend(fontsize=10)
    ax.set_aspect("equal")
    plt.tight_layout()
    plt.savefig(f"{plots_path}/scatter_identity.png", dpi=150)


    return mean_diff, std_diff, loa_lower, loa_upper

def intraclass_correlation(x, y):
    """
    Calculate the intraclass correlation coefficient (ICC) for two raters.

    Parameters
    ----------
    x, y : array-like
        Paired measurements from two raters

    Returns
    -------
    icc : float
        Intraclass correlation coefficient
    """

    # Structure data in long format
    df = pd.DataFrame({
        "subject": list(range(len(x))) * 2,
        "rater":   ["rater1"] * len(x) + ["rater2"] * len(y),
        "score":   list(x) + list(y)
    })
    
    icc_results = pg.intraclass_corr(
        data=df,
        targets="subject",
        raters="rater",
        ratings="score",
        #ci=0.95
    )

    icc_results = icc_results[["Type", "ICC", "CI95", "pval"]].copy()

    icc21 = icc_results[icc_results["Type"] == "ICC(A,1)"].iloc[0]
    
    return icc21

def icc_bias_corrected(x, y, method="median"):
    x = np.asarray(x)
    y = np.asarray(y)

    d = x - y

    if method == "mean":
        delta = np.mean(d)
    elif method == "median":
        delta = np.median(d)
    else:
        raise ValueError("method must be 'mean' or 'median'")

    # remove bias from x
    x_corr = x - delta

    # correct long-format construction (same as yours)
    df = pd.DataFrame({
        "subject": list(range(len(x))) * 2,
        "rater":   ["rater1"] * len(x) + ["rater2"] * len(y),
        "score":   list(x_corr) + list(y)
    })

    icc = pg.intraclass_corr(
        data=df,
        targets="subject",
        raters="rater",
        ratings="score"
    )

    # ICC2 = ICC(A,1)
    icc_a1 = icc[icc["Type"] == "ICC(A,1)"]

    return delta, icc_a1[["Type", "ICC", "CI95", "pval"]]

def bland_altman_analysis_robust(method1, method2):
    import numpy as np
    from scipy.stats import linregress

    mean = (method1 + method2) / 2
    diff = method1 - method2

    # Classical
    mean_diff = np.mean(diff)
    std_diff = np.std(diff, ddof=1)
    loa_upper = mean_diff + 1.96 * std_diff
    loa_lower = mean_diff - 1.96 * std_diff

    # Robust
    median_diff = np.median(diff)
    loa_lower_q = np.percentile(diff, 2.5)
    loa_upper_q = np.percentile(diff, 97.5)

    # Proportional bias
    slope, intercept, r_value, p_value, _ = linregress(mean, diff)

    print("=== Bland-Altman Results ===")
    print(f"Mean bias:              {mean_diff:.3f}")
    print(f"Median bias:            {median_diff:.3f}")
    print(f"SD of differences:      {std_diff:.3f}")
    print(f"Classical LoA:          [{loa_lower:.3f}, {loa_upper:.3f}]")
    print(f"Robust LoA (quantile):  [{loa_lower_q:.3f}, {loa_upper_q:.3f}]")
    print(f"Slope (bias vs mean):   {slope:.3f} (p={p_value:.3e})")

    return {
        "mean_bias": mean_diff,
        "median_bias": median_diff,
        "loa_classical": (loa_lower, loa_upper),
        "loa_robust": (loa_lower_q, loa_upper_q),
        "slope": slope,
        "p_value": p_value
    }