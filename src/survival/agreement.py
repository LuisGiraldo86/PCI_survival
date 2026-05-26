import numpy as np
from scipy.stats import binomtest, norm, linregress
from scipy import stats
import matplotlib.pyplot as plt
import pandas as pd
import pingouin as pg
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm

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

    # Include observed statistic in the reference distribution (Phipson & Smyth 2010)
    # so that p = (count + 1) / (B + 1), preventing p = 0.
    all_stats = np.concatenate([perm_stats, [stat_obs]])

    if alternative == "two-sided":
        p_value = np.mean(np.abs(all_stats) >= abs(stat_obs))
    elif alternative == "greater":
        p_value = np.mean(all_stats >= stat_obs)
    elif alternative == "less":
        p_value = np.mean(all_stats <= stat_obs)
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

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    d = x - y

    stat_obs = np.median(d)

    signs = rng.choice([-1, 1], size=(n_permutations, len(d)))
    perm_stats = np.median(signs * d, axis=1)

    # Include observed statistic in the reference distribution (Phipson & Smyth 2010)
    # so that p = (count + 1) / (B + 1), preventing p = 0.
    all_stats = np.concatenate([perm_stats, [stat_obs]])

    if alternative == "two-sided":
        p_value = np.mean(np.abs(all_stats) >= np.abs(stat_obs))
    elif alternative == "greater":
        p_value = np.mean(all_stats >= stat_obs)
    elif alternative == "less":
        p_value = np.mean(all_stats <= stat_obs)
    else:
        raise ValueError("Invalid alternative")

    return stat_obs, p_value

def paired_bootstrap_ci(x, y, statistic="median", 
                        n_bootstrap=20000, ci=0.95, random_state=None):
    rng = np.random.default_rng(random_state)
    d = np.asarray(x) - np.asarray(y)
    n = len(d)

    # Vectorized resampling
    resamples = rng.choice(d, size=(n_bootstrap, n), replace=True)
    
    if statistic == "median":
        stats = np.median(resamples, axis=1)
    elif statistic == "mean":
        stats = np.mean(resamples, axis=1)
    else:
        raise ValueError("Invalid statistic")

    # Clean percentile calculation
    alpha = 1 - ci
    lower, upper = np.percentile(stats, [100 * alpha / 2, 100 * (1 - alpha / 2)])

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

    q_low = np.clip(adjusted_quantile(z_low), 0, 1)
    q_high = np.clip(adjusted_quantile(z_high), 0, 1)

    lower = np.quantile(boot_stats, q_low)
    upper = np.quantile(boot_stats, q_high)

    return lower, upper

def bland_altman_analysis(method1, method2, plots_path, name1="sPCI", name2="pPCI"):
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
    plt.savefig(f"{plots_path}/scatter_identity.png", dpi=600)


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

    # long-format construction
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

def bland_altman_analysis_robust(method1, method2, plots_path=None, name1="sPCI", name2="pPCI"):
    method1 = np.asarray(method1, dtype=float)
    method2 = np.asarray(method2, dtype=float)

    mean = (method1 + method2) / 2
    diff = method1 - method2

    n = len(diff)

    # Classical
    mean_diff = np.mean(diff)
    std_diff = np.std(diff, ddof=1)
    loa_upper = mean_diff + 1.96 * std_diff
    loa_lower = mean_diff - 1.96 * std_diff

    # Robust (empirical quantiles are unstable for small n)
    median_diff = np.median(diff)
    loa_lower_q = np.percentile(diff, 2.5)
    loa_upper_q = np.percentile(diff, 97.5)

    # Proportional bias
    slope, intercept, r_value, slope_pvalue, _ = linregress(mean, diff)

    print("=== Bland-Altman Results ===")
    print(f"  Methods:                {name1} - {name2}  (n={n})")
    print(f"  Mean bias:              {mean_diff:.3f}")
    print(f"  Median bias:            {median_diff:.3f}")
    print(f"  SD of differences:      {std_diff:.3f}")
    print(f"  Classical LoA:          [{loa_lower:.3f}, {loa_upper:.3f}]")
    print(f"  Robust LoA (quantile):  [{loa_lower_q:.3f}, {loa_upper_q:.3f}]")
    print(f"  Proportional bias:      slope={slope:.3f}, intercept={intercept:.3f}, r={r_value:.3f} (p={slope_pvalue:.3e})")

    if plots_path is not None:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.scatter(mean, diff, color="steelblue", edgecolors="white", s=70, zorder=3)

        x_line = np.linspace(mean.min(), mean.max(), 200)
        ax.plot(x_line, intercept + slope * x_line, color="darkorange", linestyle="-",
                linewidth=1.2, label=f"Reg. line (slope={slope:.3f}, p={slope_pvalue:.3e})")

        ax.axhline(float(mean_diff),   color="black",     linestyle="-",  linewidth=1.5, label=f"Mean bias: {mean_diff:.2f}")
        ax.axhline(float(loa_upper),   color="crimson",   linestyle="--", linewidth=1.2, label=f"\u00b11.96 SD: [{loa_lower:.2f}, {loa_upper:.2f}]")
        ax.axhline(float(loa_lower),   color="crimson",   linestyle="--", linewidth=1.2, label="_nolegend_")
        ax.axhline(float(loa_upper_q), color="steelblue", linestyle=":",  linewidth=1.2, label=f"2.5/97.5th pct: [{loa_lower_q:.2f}, {loa_upper_q:.2f}]")
        ax.axhline(float(loa_lower_q), color="steelblue", linestyle=":",  linewidth=1.2, label="_nolegend_")
        ax.axhline(0,           color="gray",      linestyle=":",  linewidth=1.0, label="Zero line")

        ax.fill_between(ax.get_xlim(), float(loa_lower), float(loa_upper), alpha=0.05, color="crimson")

        ax.set_xlabel(f"Mean of {name1} and {name2}", fontsize=12)
        ax.set_ylabel(f"Difference ({name1} \u2212 {name2})", fontsize=12)
        ax.set_title("Bland-Altman Plot (Robust)", fontsize=14)
        ax.legend(fontsize=9)
        plt.tight_layout()
        plt.savefig(f"{plots_path}/bland_altman_robust.png", dpi=600)
        plt.close(fig)

    return {
        "mean_bias": mean_diff,
        "median_bias": median_diff,
        "loa_classical": (loa_lower, loa_upper),
        "loa_robust": (loa_lower_q, loa_upper_q),
        "slope": slope,
        "intercept": intercept,
        "r_value": r_value,
        "slope_pvalue": slope_pvalue,
    }

def bland_altman_regression_loa(method1, method2, plots_path=None, name1="sPCI", name2="pPCI"):
    """
    Bland-Altman with regression-based LoA (Bland & Altman, 1999).

    When proportional bias exists (diff correlates with mean), fixed LoA are
    misleading. This fits d ~ mean for the bias and |residuals| ~ mean for the
    SD, producing limits that vary with measurement magnitude.
    """
    method1 = np.asarray(method1, dtype=float)
    method2 = np.asarray(method2, dtype=float)

    mean_vals = (method1 + method2) / 2
    diff = method1 - method2
    n = len(diff)

    # Bias model: diff ~ mean
    slope_b, intercept_b, r_b, p_b, _ = linregress(mean_vals, diff)
    fitted_bias = intercept_b + slope_b * mean_vals
    residuals = diff - fitted_bias

    # Heteroscedasticity: |residuals| ~ mean
    # For normal residuals E[|e|] = sigma * sqrt(2/pi), so sigma = E[|e|] * sqrt(pi/2)
    abs_res = np.abs(residuals)
    slope_s, intercept_s, r_s, p_s, _ = linregress(mean_vals, abs_res)
    heteroscedastic = p_s < 0.05

    x_range = np.linspace(mean_vals.min(), mean_vals.max(), 300)
    bias_line = intercept_b + slope_b * x_range

    if heteroscedastic:
        sd_line = np.maximum(intercept_s + slope_s * x_range, 1e-6) * np.sqrt(np.pi / 2)
    else:
        sd_line = np.full(300, np.std(residuals, ddof=2))

    loa_upper_line = bias_line + 1.96 * sd_line
    loa_lower_line = bias_line - 1.96 * sd_line

    print("=== Bland-Altman: Regression-Based LoA ===")
    print(f"  n = {n}")
    print(f"  Bias:              d = {intercept_b:.3f} + {slope_b:.3f}·mean  (r={r_b:.3f}, p={p_b:.3e})")
    print(f"  Heteroscedastic:   {'Yes' if heteroscedastic else 'No'}  (|residuals| ~ mean: slope={slope_s:.4f}, p={p_s:.3e})")
    print()
    print(f"  {'mean':>6}  {'bias':>7}  {'LoA lower':>10}  {'LoA upper':>10}")
    for m in [5, 10, 15, 20, 25]:
        bias_m = intercept_b + slope_b * m
        if heteroscedastic:
            sd_m = max(intercept_s + slope_s * m, 1e-6) * np.sqrt(np.pi / 2)
        else:
            sd_m = np.std(residuals, ddof=2)
        print(f"  {m:>6}  {bias_m:>7.2f}  {bias_m - 1.96*sd_m:>10.2f}  {bias_m + 1.96*sd_m:>10.2f}")

    if plots_path is not None:
        sd_label = "heteroscedastic" if heteroscedastic else "homoscedastic"

        fig, ax = plt.subplots(figsize=(8, 5))
        ax.scatter(mean_vals, diff, color="steelblue", edgecolors="white", s=70, zorder=3, alpha=0.7)
        ax.plot(x_range, bias_line, color="black", linewidth=1.5,
                label=f"Bias: {intercept_b:.2f} + {slope_b:.3f}·mean")
        ax.plot(x_range, loa_upper_line, color="crimson", linestyle="--", linewidth=1.2,
                label="Regression LoA (±1.96 SD)")
        ax.plot(x_range, loa_lower_line, color="crimson", linestyle="--", linewidth=1.2,
                label="_nolegend_")
        ax.fill_between(x_range, loa_lower_line, loa_upper_line, alpha=0.07, color="crimson")
        ax.axhline(0, color="gray", linestyle=":", linewidth=1.0, label="Zero line")

        ax.set_xlabel(f"Mean of {name1} and {name2}", fontsize=12)
        ax.set_ylabel(f"Difference ({name1} − {name2})", fontsize=12)
        ax.set_title(f"Bland-Altman Plot – Regression LoA ({sd_label} SD)", fontsize=13)
        ax.legend(fontsize=10)
        plt.tight_layout()
        plt.savefig(f"{plots_path}/bland_altman_regression_loa.png", dpi=600)
        plt.close(fig)

    return {
        "slope_bias": slope_b,
        "intercept_bias": intercept_b,
        "r_bias": r_b,
        "p_bias": p_b,
        "heteroscedastic": heteroscedastic,
        "slope_sd": slope_s,
        "intercept_sd": intercept_s,
        "p_sd": p_s,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Level-2 Bland-Altman group analysis
# ─────────────────────────────────────────────────────────────────────────────

def _ba_stats(diff: np.ndarray, alpha: float = 0.05) -> dict:
    """Bland-Altman statistics for a single array of differences, with CIs."""
    n         = len(diff)
    mean_diff = np.mean(diff)
    std_diff  = np.std(diff, ddof=1)
    se_mean   = std_diff / np.sqrt(n)
    se_loa    = std_diff * np.sqrt(3 / n)
    t_crit    = stats.t.ppf(1 - alpha / 2, df=n - 1)
    z_loa     = 1.96

    loa_upper = mean_diff + z_loa * std_diff
    loa_lower = mean_diff - z_loa * std_diff

    return {
        "n"               : n,
        "bias"            : mean_diff,
        "sd"              : std_diff,
        "loa_upper"       : loa_upper,
        "loa_lower"       : loa_lower,
        "ci_bias_lower"   : mean_diff - t_crit * se_mean,
        "ci_bias_upper"   : mean_diff + t_crit * se_mean,
        "ci_loa_upper_lo" : loa_upper - t_crit * se_loa,
        "ci_loa_upper_hi" : loa_upper + t_crit * se_loa,
        "ci_loa_lower_lo" : loa_lower - t_crit * se_loa,
        "ci_loa_lower_hi" : loa_lower + t_crit * se_loa,
    }


def test_bias_by_group(df: pd.DataFrame, diff_col: str, group_col: str) -> dict:
    """
    Test whether systematic bias differs across groups.

    Fits: d_i = mu + beta_g * G_i + epsilon_i  (OLS, Type-II ANOVA).
    A significant beta_g means bias depends on group membership.
    """
    data = df[[diff_col, group_col]].copy()
    data[group_col] = data[group_col].astype("category")

    formula   = f"{diff_col} ~ C({group_col})"
    model     = smf.ols(formula, data=data).fit()
    anova_tbl = anova_lm(model, typ=2)
    overall_p = anova_tbl.loc[f"C({group_col})", "PR(>F)"]

    group_means = []
    for g in data[group_col].cat.categories:
        subset = data.loc[data[group_col] == g, diff_col]
        ba = _ba_stats(subset.values)
        group_means.append({
            "group"         : g,
            "n"             : ba["n"],
            "bias"          : ba["bias"],
            "ci_bias_lower" : ba["ci_bias_lower"],
            "ci_bias_upper" : ba["ci_bias_upper"],
        })

    return {
        "model"       : model,
        "anova_table" : anova_tbl,
        "group_means" : pd.DataFrame(group_means).set_index("group"),
        "overall_pval": overall_p,
    }


def test_variance_by_group(df: pd.DataFrame, diff_col: str, group_col: str) -> dict:
    """
    Test whether the spread of differences is homogeneous across groups.

    Levene's test (robust to non-normality) and Bartlett's test (optimal
    under normality).
    """
    groups     = df[group_col].unique()
    group_data = [df.loc[df[group_col] == g, diff_col].values.astype(float)
                  for g in groups]

    lev_stat, lev_p = stats.levene(*group_data, center="median")
    try:
        bar_stat, bar_p = stats.bartlett(*group_data)
    except ValueError:
        bar_stat, bar_p = float("nan"), float("nan")

    group_sd = pd.DataFrame({
        "group": groups,
        "n"    : [len(g) for g in group_data],
        "sd"   : [np.std(g, ddof=1) for g in group_data],
        "var"  : [np.var(g, ddof=1) for g in group_data],
    }).set_index("group")

    return {
        "levene"  : {"statistic": lev_stat,  "p_value": lev_p},
        "bartlett": {"statistic": bar_stat, "p_value": bar_p},
        "group_sd": group_sd,
    }


def group_loa(
    df: pd.DataFrame,
    diff_col: str,
    group_col: str,
    alpha: float = 0.05,
) -> pd.DataFrame:
    """Compute full Bland-Altman statistics (bias, SD, LoA + CIs) per group."""
    records = []
    for g, sub in df.groupby(group_col):
        ba = _ba_stats(sub[diff_col].values, alpha=alpha)
        ba["group"] = g
        records.append(ba)

    cols = [
        "group", "n", "bias",
        "ci_bias_lower", "ci_bias_upper",
        "sd",
        "loa_upper", "ci_loa_upper_lo", "ci_loa_upper_hi",
        "loa_lower", "ci_loa_lower_lo", "ci_loa_lower_hi",
    ]
    return pd.DataFrame(records)[cols].set_index("group")


def bland_altman_group_analysis(
    df: pd.DataFrame,
    method1_col: str,
    method2_col: str,
    group_col: str,
    alpha: float = 0.05,
) -> dict:
    """
    Level-2 Bland-Altman group analysis (Bland & Altman 1999).

    2a. Linear model test: does bias differ across groups?
    2b. Levene + Bartlett tests: does variance differ across groups?
    Plus group-specific LoA with confidence intervals.
    """
    data          = df.copy()
    data["_diff"] = data[method1_col] - data[method2_col]
    data["_mean"] = (data[method1_col] + data[method2_col]) / 2

    bias_res = test_bias_by_group(data, "_diff", group_col)
    var_res  = test_variance_by_group(data, "_diff", group_col)
    loa_tbl  = group_loa(data, "_diff", group_col, alpha=alpha)

    sep = "─" * 60
    print(sep)
    print(f"BLAND-ALTMAN GROUP ANALYSIS  [{group_col}]")
    print(sep)

    print("\n[2a] BIAS DIFFERENCES ACROSS GROUPS")
    print(f"     d ~ C({group_col}),  overall p = {bias_res['overall_pval']:.4f}", end="  ")
    print("*** significant" if bias_res["overall_pval"] < alpha else "(not significant)")
    print("\n     Group-level bias (95 % CI):")
    print(bias_res["group_means"].round(3).to_string())

    lev = var_res["levene"]
    bar = var_res["bartlett"]
    print(f"\n[2b] VARIANCE DIFFERENCES ACROSS GROUPS")
    print(f"     Levene   W = {lev['statistic']:.4f},  p = {lev['p_value']:.4f}", end="  ")
    print("*** significant" if lev["p_value"] < alpha else "(not significant)")
    print(f"     Bartlett K = {bar['statistic']:.4f},  p = {bar['p_value']:.4f}", end="  ")
    print("*** significant" if bar["p_value"] < alpha else "(not significant)")
    print("\n     Group-level SD:")
    print(var_res["group_sd"].round(3).to_string())

    print(f"\nGROUP-SPECIFIC LoA  (alpha = {alpha})")
    print(loa_tbl.round(2).to_string())
    print(sep)

    return {
        "data"         : data,
        "bias_test"    : bias_res,
        "variance_test": var_res,
        "loa_table"    : loa_tbl,
    }