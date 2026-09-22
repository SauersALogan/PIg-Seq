# =============================================================================
# Script:       07a_fraction_rf_classifier.py
# Author:       Katharine Z. Coyte
# Repository:   github.com/SauersALogan/PIg-Seq
# Description:  Binary random forest classifier predicting whether a MAG's
#               relative expression profile came from the IgA-positive or
#               IgA-negative fraction. DNA relative abundance profiles used
#               as control. 100 iterations with ~20% bin-level holdout.
# Dependencies: numpy, pandas, scikit-learn, scipy, openpyxl
# Input:        RNA_gene_count_controlled_subsampled.csv (https://doi.org/10.48420/33951214)
#               DNA_gene_count_controlled_subsampled.csv (https://doi.org/10.48420/33951214)
#               bin_IDs.xlsx (https://doi.org/10.48420/33951214)
#               Taxon_scores_02-10-2025.csv (https://doi.org/10.48420/33951214)
# Output:       accuracy_data.csv
#               roc_curve_dna.csv
#               roc_curve_rna.csv
#               auc_summary.csv
# =============================================================================

import os
import numpy as np
import pandas as pd
import random
from tqdm import tqdm
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, roc_curve, auc
import scipy

# ---- Paths (update for your system) ----
data_dir = "/path/to/data"
out_dir = "/path/to/output"

# =====================================================================
# Functions
# =====================================================================
def create_rf_matrix(file_name, all_scores):
    df = pd.read_csv(file_name, sep="\t")

    for ix in df.index:
        df.loc[ix, "my_bin"] = df.loc[ix, "Chr"].split("_")[0]

    filter = df["product"].str.contains("ibosomal")
    df = df[~filter]

    save_all = pd.DataFrame()

    for my_bin in np.unique(df.my_bin):
        foo = df.loc[df.my_bin == my_bin, :].copy()
        foo = foo.rename(columns={"product": "my_product"})
        foo = foo.loc[foo.my_product != "hypothetical protein", :]
        foo = foo.drop(["Geneid", "Chr", "Length"], axis=1)
        foo = foo.iloc[:, :-1]
        foo = foo.groupby(by="my_product").sum()
        foo.columns = foo.columns + "_" + my_bin
        save_all = pd.concat([save_all, foo], axis=1)

    save_all = save_all.fillna(0)

    filter = save_all.columns.str.contains("_Native_")
    for_random_forest = save_all.loc[:, ~filter].T

    for ix in for_random_forest.index:
        tmp = ix.split("_")
        lookup_index = tmp[0] + "_Native_" + tmp[2]
        for_random_forest.loc[ix, "current_bin"] = tmp[2]
        try:
            for_random_forest.loc[ix, "score_lookup"] = all_scores.loc[
                lookup_index, "scores"]
        except:
            "nothing"
    return for_random_forest


def train_and_evaluate(X, y, num_runs=100):
    accuracies = []
    roc_store = []

    unique_bins = np.unique(X.current_bin)

    for ix in tqdm(range(num_runs)):
        random_state = random.randint(0, 100000)

        test_bins = pd.Series(unique_bins).sample(20,
            random_state=random_state)

        X_test = X.loc[X.current_bin.isin(test_bins), :]
        X_test = X_test.iloc[:, :-2].astype(float)
        y_test = y.loc[X_test.index]

        X_train = X.loc[~X.current_bin.isin(test_bins), :]
        X_train = X_train.iloc[:, :-2].astype(float)
        y_train = y.loc[X_train.index]

        clf = RandomForestClassifier(random_state=random_state,
            n_estimators=100)
        clf.fit(X_train, y_train)

        y_pred = clf.predict(X_test)
        y_proba = clf.predict_proba(X_test)

        accuracy = accuracy_score(y_test, y_pred)
        accuracies.append(accuracy)

        fpr, tpr, thresholds = roc_curve(y_test, y_proba[:, 1])
        roc = auc(fpr, tpr)
        roc_store.append(roc)

    return accuracies, roc_store


def single_rf(X, y, random_state=0):
    unique_bins = np.unique(X.current_bin)
    test_bins = pd.Series(unique_bins).sample(20, random_state=random_state)

    X_test = X.loc[X.current_bin.isin(test_bins), :]
    X_test = X_test.iloc[:, :-2].astype(float)
    y_test = y.loc[X_test.index]

    X_train = X.loc[~X.current_bin.isin(test_bins), :]
    X_train = X_train.iloc[:, :-2].astype(float)
    y_train = y.loc[X_train.index]

    clf = RandomForestClassifier(random_state=random_state, n_estimators=100)
    clf.fit(X_train, y_train)

    y_pred = clf.predict(X_test)
    y_proba = clf.predict_proba(X_test)

    fpr, tpr, thresholds = roc_curve(y_test, y_proba[:, 1])
    roc = auc(fpr, tpr)

    return fpr, tpr, thresholds, roc, y_test, y_pred


# =====================================================================
# Data loading
# =====================================================================
bin_list = pd.read_excel(os.path.join(data_dir, "bin_IDs.xlsx"))
for ix in bin_list.index:
    tmp = bin_list.loc[ix, "Bin Loci ID"].split(".locus")[0]
    bin_list.loc[ix, "new_bin_name"] = tmp
bin_list = bin_list.set_index("Bin name")

scores = pd.read_csv(os.path.join(data_dir, "Taxon_scores_02-10-2025.csv"))
scores = scores.set_index("Unnamed: 0")

renamed_scores = pd.concat([bin_list, scores], axis=1)

all_scores = pd.DataFrame()
for ix in renamed_scores.index:
    tmp = renamed_scores.loc[ix, "S3":"S453"].T
    mag = renamed_scores.loc[ix, "new_bin_name"]
    tmp.index = tmp.index + "_Native_" + mag
    tmp = pd.DataFrame(tmp)
    tmp.columns = ["scores"]
    tmp["bin_id"] = mag
    all_scores = pd.concat([all_scores, tmp], axis=0)

all_scores = all_scores.dropna()

## ---- Bin IgA score classification ----
mean = all_scores.scores.mean()
std = all_scores.scores.std()
all_scores.loc[all_scores.scores > mean + 1.5 * std, "scores"] = 1
all_scores.loc[all_scores.scores < mean - 1.5 * std, "scores"] = 2
all_scores.loc[
    (all_scores.scores != 1) & (all_scores.scores != 2), "scores"] = 0

# =====================================================================
# Build feature matrices
# =====================================================================
for_rf_rna = create_rf_matrix(
    os.path.join(data_dir, "RNA_gene_count_controlled_subsampled.csv"),
    all_scores)
for_rf_dna = create_rf_matrix(
    os.path.join(data_dir, "DNA_gene_count_controlled_subsampled.csv"),
    all_scores)

X_dna = for_rf_dna.dropna()
X_rna = for_rf_rna.dropna()

y_dna = X_dna.index.to_series().apply(lambda x: 1 if '_Pos' in x else 0)
y_rna = X_rna.index.to_series().apply(lambda x: 1 if '_Pos' in x else 0)

# =====================================================================
# Run 100 bootstrap iterations
# =====================================================================
N = 100
DNA_accuracies, DNA_roc = train_and_evaluate(X_dna, y_dna, num_runs=N)
RNA_accuracies, RNA_roc = train_and_evaluate(X_rna, y_rna, num_runs=N)

# =====================================================================
# Single RF for ROC curves
# =====================================================================
fpr_dna, tpr_dna, _, roc_dna, _, _ = single_rf(X_dna, y_dna, random_state=0)
fpr_rna, tpr_rna, _, roc_rna, _, _ = single_rf(X_rna, y_rna, random_state=0)

# =====================================================================
# Export
# =====================================================================
ROC_data = pd.DataFrame({
    'accuracy': DNA_roc + RNA_roc,
    'Type': ['DNA'] * len(DNA_roc) + ['RNA'] * len(RNA_roc)
})
ROC_data.to_csv(os.path.join(out_dir, "accuracy_data.csv"), index=False)

pd.DataFrame({'fpr': fpr_dna, 'tpr': tpr_dna}).to_csv(
    os.path.join(out_dir, "roc_curve_dna.csv"), index=False)
pd.DataFrame({'fpr': fpr_rna, 'tpr': tpr_rna}).to_csv(
    os.path.join(out_dir, "roc_curve_rna.csv"), index=False)

auc_summary = pd.DataFrame({
    'Type': ['DNA', 'RNA'],
    'mean_AUC': [np.mean(DNA_roc), np.mean(RNA_roc)],
    'std_AUC': [np.std(DNA_roc), np.std(RNA_roc)],
    'median_AUC': [np.median(DNA_roc), np.median(RNA_roc)]
})
auc_summary.to_csv(os.path.join(out_dir, "auc_summary.csv"), index=False)

t_stat, p_value = scipy.stats.ranksums(DNA_roc, RNA_roc)
print(f"Wilcoxon rank-sum: t={t_stat:.4f}, p={p_value:.6f}")
print("Export complete.")
