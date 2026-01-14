import pandas as pd
from scipy.stats import levene
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd

def run_anova_with_posthoc(df, dv, iv):
    """
    Performs a standard statistical pipeline: Levene's Test for homogeneity 
    of variance, Type II ANOVA, and Tukey's HSD post-hoc test.

    Args:
        df (pd.DataFrame): The dataset containing the data.
        dv (str): The column name of the dependent variable (continuous).
        iv (str): The column name of the independent variable (categorical).

    Returns:
        pd.DataFrame: A dataframe containing Tukey HSD pairwise comparisons 
        augmented with ANOVA and Levene test p-values.

    Example:
        >>> data = {
        ...     'diet': ['A', 'A', 'B', 'B', 'C', 'C'],
        ...     'weight_loss': [5, 4, 8, 9, 2, 3]
        ... }
        >>> df = pd.DataFrame(data)
        >>> results = analyze_variance_pipeline(df, dv='weight_loss', iv='diet')
        >>> print(results)
    """
    # 1. Homogeneity of Variance (Levene's Test)
    # Group the DV by the IV and unpack the resulting arrays into levene()
    groups = [group[dv].values for _, group in df.groupby(iv)]
    levene_stat, levene_p = levene(*groups)
    
    # 2. ANOVA (Type II)
    # Using Type II as it is generally more robust for unbalanced designs
    formula = f"{dv} ~ C({iv})"
    model = ols(formula, data=df).fit()
    anova_table = sm.stats.anova_lm(model, typ=2)
    anova_p = anova_table['PR(>F)'].iloc[0]

    # 3. Post-hoc (Tukey HSD)
    tukey = pairwise_tukeyhsd(df[dv], df[iv])
    
    # Convert Tukey summary to DataFrame
    tukey_df = pd.DataFrame(
        data=tukey.summary().data[1:], 
        columns=tukey.summary().data[0]
    )
    
    # 4. Append statistical metadata
    tukey_df["anova_p"] = anova_p
    tukey_df["levene_p"] = levene_p
    tukey_df["is_homoscedastic"] = levene_p > 0.05
    
    return tukey_df

# Example of applying it to a grouped dataframe (like your original cell_group)
# results = df.groupby('cell_type').apply(lambda x: analyze_variance_pipeline(x, 'value', 'age'))
