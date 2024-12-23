import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy
import anndata as ad
import sklearn.mixture
import os

import filenames

from IPython.display import display

# Note this writes to disk in the current directory!
def get_targeted_dnam():
    df_testset, ages_testset = preprocess_test_set()

    ad10_fname = "./ad10.h5ad"
    ad100_fname = "./ad100.h5ad"
    if os.path.exists(ad10_fname and ad100_fname):
        ad10_all = ad.read_h5ad(ad10_fname)
        ad100_all = ad.read_h5ad(ad100_fname)
    else:
        ad10_all = preprocess_cgmatrix(filenames.TARGETED_10X_TSV, filenames.TARGETED_FROG_ANNOTATION)
        ad100_all = preprocess_cgmatrix(filenames.TARGETED_100X_TSV, filenames.TARGETED_FROG_ANNOTATION)
        ad10_all.write("./ad10.h5ad")
        ad100_all.write("./ad100.h5ad")
    
    ad10 = ad10_all[~ad10_all.obs['is_genetic'], :]
    ad100 = ad100_all[~ad100_all.obs['is_genetic'], :]
    return ad10_all, ad100_all, ad10, ad100, df_testset, ages_testset

def show_ad(ad10):
    display(pd.concat([ad10.to_df(), pd.DataFrame(ad10.obs[['spearman', 'spearman_p', 'pearson', 'pearson_p']]) ] , axis=1))

def get_ages(ad10):
    return list(ad10.var['Age (years)'])

def preprocess_cgmatrix(data_fname, annotation_fname):
    df_annotation = pd.read_csv(annotation_fname, sep='\t')  # df_annotation is a master metadata spreadsheet I collected
    df_annotation.index = df_annotation["Unique Frog ID"]
    df_annotation = df_annotation[df_annotation["Has TBSeq data?"] == True]  # ignore annotations we have no data for

    df = pd.read_csv(data_fname, index_col="Site", sep='\t')  # TBSeq data matrix ("Site" column is just the index)
    df = df[list(df_annotation["Unique Frog ID"])]  # SORTS df by age since df_annotation is sorted by age
    assert(set(list(df.columns)) == set(list(df_annotation["Unique Frog ID"]))) # Requires column names of df == df_annotation["Unique Frog ID"]
    
    ad_out = ad.AnnData(df, var=df_annotation, dtype=np.float64)  # create anndata object. per-frog metadata in var.
    ages_cpg_sorted = list(ad_out.var["Age (years)"])   # ASSUMES YEAR IS 365.25 DAYS
    assert(ages_cpg_sorted == sorted(ages_cpg_sorted))

    compute_obs_attributes(ad_out)

    # SORT rows (sites) by age correlation
    sorted_indices = np.abs(ad_out.obs['spearman']).sort_values(ascending=False).index
    ad_out = ad_out[sorted_indices, :]

    assert(ages_cpg_sorted == sorted(ages_cpg_sorted))
    assert(all([a == b for a,b in zip(ad_out.var_names, df_annotation["Unique Frog ID"])]))
    assert(all([a == b for a,b in zip(ages_cpg_sorted, df_annotation["Age (years)"])]))

    return ad_out

def compute_obs_attributes(ad_out):
    ages_cpg_sorted = list(ad_out.var["Age (years)"])   # ASSUMES YEAR IS 365.25 DAYS
    assert(ages_cpg_sorted == sorted(ages_cpg_sorted))

    # metadata per-SITE (row) (anndata.obs)
    ad_out.obs['Site'] = list(ad_out.obs_names)
    ad_out.obs['CHR'] = [int(x.split(':')[0].split('Chr')[1]) for x in ad_out.obs_names]  # used in manhattan plot function
    ad_out.obs['BP'] = [int(x.split(':')[1]) for x in ad_out.obs_names]
    ad_out.obs['spearman'] = [scipy.stats.spearmanr(ages_cpg_sorted, cpgrow)[0] for cpgrow in ad_out.X]  # age correlation
    ad_out.obs['spearman_p'] = [scipy.stats.spearmanr(ages_cpg_sorted, cpgrow)[1] for cpgrow in ad_out.X]
    ad_out.obs['pearson'] = [scipy.stats.pearsonr(ages_cpg_sorted, cpgrow)[0] for cpgrow in ad_out.X]
    ad_out.obs['pearson_p'] = [scipy.stats.pearsonr(ages_cpg_sorted, cpgrow)[1] for cpgrow in ad_out.X]

    # determine which sites are genetic (calculation only based on non-tadpoles)
    ad_out_notadpole = ad_out.copy()
    ad_out_notadpole = ad_out_notadpole[:, ~ad_out_notadpole.var['Is tadpole?']]
    has_zero_meth_rows = (ad_out_notadpole.X < 0.006).any(axis=1)  # rows with at least one value less than 0.006
    is_genetic = []
    for cpgrow, has_zero_meth in zip(ad_out_notadpole.X, has_zero_meth_rows):
        # A genetic site is defined as at least one cpg of 0, and has methylation values that cluster well into 3 gaussians.
        is_genetic.append(has_zero_meth and abs(fit_gmm(cpgrow,3) - fit_gmm(cpgrow,1)) > 90)
    ad_out.obs['is_genetic'] = is_genetic

def fit_gmm(series, n_components):
    gmm = sklearn.mixture.GaussianMixture(n_components=n_components, random_state=0)
    gmm.fit(series.reshape(-1, 1))
    bic = gmm.bic(series.reshape(-1, 1))
    return bic
    
def preprocess_test_set():
    df100_test_raw = pd.read_csv(filenames.TARGETED_TEST_100X_TSV, sep='\t')  # 'Site' column + 192 frog DNAm columns

    ages_test_set_sorted = [2.338124572,2.340862423,2.343600274,3.608487337,4.692676249,4.813141684,4.815879535,6.984257358,7.909650924,7.912388775,8.547570157,8.791238877,10.98151951,10.98425736,10.98699521,10.98973306]
    test_desired_col_order = ["134_XT","137_XT","138_XT","10_XT","26_XT","14_XT","16_XT","61_XT","100_XT","108_XT","30_XT","38_XT","40_XT","43_XT","44_XT","59_XT"]

    df100_test_raw.set_index("Unnamed: 0", inplace=True)
    df100_test_raw = df100_test_raw[test_desired_col_order]
    return df100_test_raw, ages_test_set_sorted
