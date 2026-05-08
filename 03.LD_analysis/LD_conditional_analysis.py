import pandas as pd
import numpy as np
from scipy.stats import norm
from collections import defaultdict
from pathlib import Path



# Statistics helpers
def calculate_p_from_z(z):
    """
    Two-sided P-value from a Z-score.
    Uses asymptotic approximation for |Z| >= 8 to avoid floating-point underflow.
    """
    abs_z = abs(z)
    if abs_z < 8:
        return 2 * (1 - norm.cdf(abs_z))
    return max(1e-300, (2 / (np.sqrt(2 * np.pi) * abs_z)) * np.exp(-0.5 * abs_z**2))


def pseudo_conditional_z(lead_z, ld_pairs, selected_variants):
    """
    Compute the pseudo-conditional Z-score for a test variant given already-
    selected variants it has moderate LD with.

    Params:
    lead_z: absolute Z-score of the test variant
    ld_pairs: list of (variant_id, r2) for selected variants in moderate LD
    selected_variants : dict mapping variant_id -> absolute Z-score

    Returns: Conditional Z-score if > 0, else 0.
    """
    sum_term = sum(selected_variants[snp] * np.sqrt(r2) for snp, r2 in ld_pairs)
    delta_z = lead_z - sum_term
    prod_term = np.prod([np.sqrt(max(1 - r2, 1e-6)) for _, r2 in ld_pairs])
    return max(delta_z, 0) / prod_term

# LD loading & indexing
def load_ld_dictionaries(ld_file, r2_lower_bound, r2_high_threshold):
    """
    Read a pairwise LD file and return three dictionaries keyed by variant ID:
        ld_dict: all pairs with R² >= r2_lower_bound
        high_ld_dict: pairs with R² >  r2_high_threshold
        moderate_ld_dict: pairs with R² [r2_lower_bound, r2_high_threshold)
    """
    print(f"Loading LD file: {ld_file}")
    raw = pd.read_csv(ld_file, delim_whitespace=True)

    # Deduplicate symmetric pairs
    raw["_pair"] = raw.apply(lambda r: tuple(sorted([r["SNP_A"], r["SNP_B"]])), axis=1)
    raw = raw.drop_duplicates(subset=["_pair", "R2"]).drop(columns=["_pair"])

    ld_dict = defaultdict(list)
    for _, row in raw.iterrows():
        a, b, r2 = row["SNP_A"], row["SNP_B"], row["R2"]
        if r2 >= r2_lower_bound:
            ld_dict[a].append((b, r2))
            ld_dict[b].append((a, r2))

    ld_dict = dict(ld_dict)

    high_ld_dict = {
        snp: [(v, r2) for v, r2 in pairs if r2 > r2_high_threshold]
        for snp, pairs in ld_dict.items()
    }
    moderate_ld_dict = {
        snp: [(v, r2) for v, r2 in pairs if r2_lower_bound <= r2 < r2_high_threshold]
        for snp, pairs in ld_dict.items()
    }
    print(f"LD dictionary built for {len(ld_dict):,} variants.")
    return ld_dict, high_ld_dict, moderate_ld_dict

# Summary statistics loading

def load_sumstats(path):
    """
    Load summary statistics, drop variants not found in any study, and add an absolute Z column.
    """
    df = pd.read_csv(path, sep="\t")
    not_found = df[df["study_condition_on"] == "not_found"]
    if len(not_found):
        print(f"NOTE:Dropping {len(not_found)} variants with no summary stats.")
    df = df[~df["ID"].isin(not_found["ID"])].copy()
    df["abs_Z_adj"] = df["Z_adj"].astype(float).abs()
    return df


# Stepwise conditional selection

def run_conditional_selection(sumstats, step0_variants, high_ld_dict, moderate_ld_dict, p_cond_threshold):
    """
    Greedy stepwise pseudo-conditional fine-mapping.

    Returns: variant_step_status : dict mapping variant_id -> (step, status_str, p_cond_or_dot, ld_pairs_or_dot)
    """
    selected_variants = {}   # variant_id -> abs_Z
    skipped_due_to_ld = defaultdict(list)
    variant_step_status = {}

    #  Step 0: seed with index variants ----
    for _, row in sumstats.iterrows():
        if row["ID"] in step0_variants:
            selected_variants[row["ID"]] = row["abs_Z_adj"]

    for v in selected_variants:
        variant_step_status[v] = (0, "selected", ".", ".")
    print(f"  Step 0: {len(selected_variants)} index variant(s) seeded.")

    # interative selection
    candidate_snps = sumstats.copy()
    step = 1

    while not candidate_snps.empty:
        # Remove already-selected variants
        candidate_snps = candidate_snps[~candidate_snps["ID"].isin(selected_variants)]
        candidate_snps = candidate_snps.sort_values("abs_Z_adj", ascending=False)
        print(f"  Step {step}: {len(candidate_snps)} candidate SNPs remaining.")

        # Prune variants in high LD with any selected variant
        high_ld_variants = set()
        for selected_snp in selected_variants:
            for linked_variant, r2 in high_ld_dict.get(selected_snp, []):
                if linked_variant in candidate_snps["ID"].values:
                    skipped_due_to_ld[linked_variant].append((step, selected_snp, r2))
                    variant_step_status[linked_variant] = (step, "skipped_r2>0.6_with_selected", ".", ".")
                    high_ld_variants.add(linked_variant)

        candidate_snps = candidate_snps[~candidate_snps["ID"].isin(high_ld_variants)]
        if candidate_snps.empty:
            break

        lead_row = candidate_snps.iloc[0]
        lead_variant = lead_row["ID"]
        lead_z = lead_row["abs_Z_adj"]
        print(f"  Step {step}: top candidate = {lead_variant}")

        # Pseudo-conditional analysis
        moderate_ld_pairs = [
            (snp, r2)
            for snp, r2 in moderate_ld_dict.get(lead_variant, [])
            if snp in selected_variants
        ]

        if moderate_ld_pairs:
            z_cond = pseudo_conditional_z(lead_z, moderate_ld_pairs, selected_variants)
            p_cond = calculate_p_from_z(z_cond)
        else:
            p_cond = calculate_p_from_z(lead_z)

        if p_cond < p_cond_threshold:
            selected_variants[lead_variant] = lead_z
            variant_step_status[lead_variant] = (step, "selected", p_cond, moderate_ld_pairs)
            print(f"  Step {step}: SELECTED {lead_variant} (P_cond={p_cond:.2e}, conditioned on {len(moderate_ld_pairs)} variant(s))")
        else:
            variant_step_status[lead_variant] = (step, "failed_condition", p_cond, moderate_ld_pairs)
            print(f"  Step {step}: SKIPPED {lead_variant} (P_cond={p_cond:.2e})")
            candidate_snps = candidate_snps[candidate_snps["ID"] != lead_variant]

        step += 1

    return variant_step_status


# Result annotation

def annotate_results(sumstats, variant_step_status, ld_dict, r2_lower_bound):
    """Attach status columns and full LD neighbourhood to the summary stats table."""

    def get_field(variant_id, idx):
        return variant_step_status.get(variant_id, (None, None, None, None))[idx]

    result = sumstats.copy()
    result["tested_step"]          = result["ID"].apply(lambda x: get_field(x, 0))
    result["status"]               = result["ID"].apply(lambda x: get_field(x, 1))
    result["P_cond"]               = result["ID"].apply(lambda x: get_field(x, 2))
    result["condition_on_var(r2)"] = result["ID"].apply(lambda x: get_field(x, 3))

    all_ld = {
        snp: [(v, r2) for v, r2 in pairs if r2 >= r2_lower_bound]
        for snp, pairs in ld_dict.items()
    }
    result[f"all_ld(R2>{r2_lower_bound})"] = result["ID"].apply(lambda x: all_ld.get(x, ""))
    return result



# run analysis per trait 

def run_trait_analysis(trait, step0_variants, ld_dict, high_ld_dict, moderate_ld_dict,
                       r2_lower_bound, p_cond_threshold, r2_high_ld_threshold, output_dir):

    print(f"\n{'-'*30}\nRunning trait: {trait.upper()}\n{'-'*30}")

    sumstats_path = Path("data/summary_stats") / f"combined_sumstats_{trait}_with_Zadj.tsv"
    sumstats = load_sumstats(sumstats_path)
    print(f"  Loaded {len(sumstats):,} variants after filtering.")

    variant_step_status = run_conditional_selection(
        sumstats, step0_variants, high_ld_dict, moderate_ld_dict, p_cond_threshold
    )

    result = annotate_results(sumstats, variant_step_status, ld_dict, r2_lower_bound)

    print(f"\n  Status breakdown for {trait.upper()}:")
    print(result["status"].value_counts().to_string())

    Path(output_dir).mkdir(parents=True, exist_ok=True)
    out_path = Path(output_dir) / f"{trait}_conditional_result_{r2_lower_bound}-{r2_high_ld_threshold}.csv"
    result.to_csv(out_path, index=False)
    print(f"  Results written to: {out_path}")

    return result


def main():
    
    traits              = ["cd", "uc", "ibd"]
    r2_lower_bound      = 0.001   
    r2_high_threshold   = 0.6     
    p_cond_threshold    = 3e-7   
    # paths
    ld_file             = Path("data/ukbb_wgs_5mb.ld")
    index_variants_file = Path("data/Index_variants_for_conditional_analysis.xlsx")
    output_dir          = Path("conditional_result")

    
    ld_dict, high_ld_dict, moderate_ld_dict = load_ld_dictionaries(
        ld_file, r2_lower_bound, r2_high_threshold
    )
    step0_variants = pd.read_excel(index_variants_file)["ID"].tolist()
    print(f"Loaded {len(step0_variants)} index variants for step 0.")

    
    results = {}
    for trait in traits:
        results[trait] = run_trait_analysis(
            trait, step0_variants,
            ld_dict, high_ld_dict, moderate_ld_dict,
            r2_lower_bound, p_cond_threshold, r2_high_threshold, output_dir
        )

    print("\nAll traits completed.")
    return results


if __name__ == "__main__":
    main()