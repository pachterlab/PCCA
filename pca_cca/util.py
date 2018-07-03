"""Utility module for FAA scripts.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import gzip
import functools

import concurrent.futures
import numpy as np
import pandas as pd
from time import time
from sklearn import preprocessing

from IPython import embed

# pylint: disable=invalid-name
def quantile_normalize(df):
  """Quantile normalizes a pandas DataFrame.

  https://stackoverflow.com/questions/37935920/quantile-normalization-on-pandas-dataframe

  Args:
    df: A pandas DataFrame
  Returns:
    The dataframe where all columns have the same distribution.
  """
  rank_mean = df.stack().groupby(
      df.rank(method='first').stack().astype(int)).mean()
  return df.rank(method='min').stack().astype(int).map(rank_mean).unstack()

def standardize_df(df, axis=0):
  """Standardizes a df such that columns have mean 0 and variance 1.

  Args:
    df: A pandas DataFrame.
  """
  df.loc[:] = preprocessing.scale(df.values)

def gtf_to_trans_gene_hash(gtf_file_name):
  """Parses a GTF file to produce a dict mapping transcripts to genes.
  """
  trans_gene_hash = dict()
  if gtf_file_name[-3:] == '.gz':
    open_fn = functools.partial(gzip.open, mode='rt')
  else:
    open_fn = open
  with open_fn(gtf_file_name) as gtf_file:
    for line in gtf_file:
      if line[0:2] == '##':  # skip comments
        continue
      l = line.strip().split()
      if l[2] == 'transcript':
        gene_name = l[9][1:-2]
        trans_name = l[11][1:-2]
        trans_gene_hash[trans_name] = gene_name
  return trans_gene_hash

def trans_to_gene_df(trans_df, trans_gene_hash):
  """Turns a df of transcript levels into gene levels based on trans_gene_bash.

  Args:
    trans_df: dataframe of samples x transcripts.
    trans_gene_hash: a dictionary mapping columns of trans_df (transcripts) to
       their gene names.
  Returns:
    A dataframe of samples by genes, where gene level is sum of its transcripts.
  """
  genes = set([trans_gene_hash[t] for t in trans_df.columns])
  gene_df = pd.DataFrame(np.zeros((len(trans_df.index), len(genes))),
                         index=trans_df.index, columns=genes)
  for transcript in trans_df.columns:
    gene_df[trans_gene_hash[transcript]] += trans_df[transcript]
  return gene_df

def _keep_gene(chrm, start, stop):
  """Determines whether a gene at chr start stop is autosomal non-MHC.
  """
  mhc_s = 29000000
  mhc_e = 33000000
  chrm_not_six = (chrm in ['chr' + str(i) for i in
                           list(range(1, 6)) + list(range(7, 23))])
  six_not_mhc = ((chrm == 'chr6') and
                 ((int(stop) < mhc_s)|(int(start) > mhc_e)))
  return chrm_not_six | six_not_mhc

def get_genes_to_keep(gtf_file_name, annot='gene'):
  """Returns list of genes in gtf file that are autosomal non-MHC.
  """
  keep = []
  if gtf_file_name[-3:] == '.gz':
    open_fn = functools.partial(gzip.open, mode='rt')
  else:
    open_fn = open
  with open_fn(gtf_file_name) as gtf_file:
    for line in gtf_file:
      if line[0:2] == '##':  # skip comments
        continue
      l = line.strip().split()
      chrm, feature, start, stop = l[0], l[2], l[3], l[4]
      if feature == annot:
        if _keep_gene(chrm, start, stop):
          gene_name = l[9][1:-2]
          keep.append(gene_name)
  return keep


def _cv_worker(df, sample, pc_loc, dim_y1, dim_y2, dim_z, method):
  y1_i = df.loc[sample].values
  df_no_i = df.drop(sample)
  u_y2_no_i = pd.read_csv(
      pc_loc, sep=' ', index_col=0, usecols=range(1, dim_y2+2))
  df_no_i = df_no_i.loc[u_y2_no_i.index]  # be careful with index
  u_y1_no_i, _ = pca(df_no_i.values, dim_y1)
  if method == 'CCA':
    l_y1, l_y2, _ = lr_cca(u_y1_no_i, u_y2_no_i.values, dim_z=dim_z)
    df_i_full, tr_error, te_error = project_i_cca(
      df_no_i, u_y1_no_i, l_y1, y1_i)
  elif method == 'regression':
    l_y1, l_y2 = lr_regression(u_y1_no_i, u_y2_no_i.values)
  return df_i_full, tr_error, te_error


# TODO(brielin): Consider making this function less specific.
def build_cv_df(df, pc_locs, dim_y1, dim_y2, dim_z, threads=1, method='cca'):
  """Builds a df where each row is the projection of a held out ind.

  Args:
    df: A pandas DataFrame of samples x features
    pc_locs: A dict mapping sample names to .eigenvec filenames with (genotype)
      PCs from each holdout.
    threads: number of threads to use, default is 1.
  Returns:
    A dataframe where each row is the projection of the sample first into the df
    (expression) sidee of faa with the (genotype) PCs stored on disk, then back
    into faa space.
  """
  args = [(df, sample, pc_loc, dim_y1, dim_y2, dim_z, method)
          for sample, pc_loc in pc_locs.items()]
  with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
    result = list(executor.map(_cv_worker, *zip(*args)))
  projections, tr_error, te_error = zip(*result)
  index = [arg[1] for arg in args]
  cv_df = pd.DataFrame(list(projections), index=index)
  error_df = pd.DataFrame({'train': tr_error, 'test': te_error}, index=index)
  cv_df = cv_df.loc[df.index]
  error_df = error_df.loc[df.index]
  return cv_df, error_df


def _perm_worker(arr_x, u_x, u_y, dim_z, n_perm, true_var_exp, method):
  n_ind, n_genes = arr_x.shape
  perm_mean = np.zeros(n_genes)
  perm_mean_sq = np.zeros(n_genes)
  tail_counts = np.zeros(n_genes)
  seed = int(divmod(time(), 1)[1]*1e8)
  np.random.seed(seed)
  # U, L, V = np.linalg.svd(arr_x, full_matrices=False)
  # nsvs = np.sum(L > 1e-11)
  # U = U[:, 0:nsvs]
  # L = L[0:nsvs]
  # V = V[0:nsvs].T
  # SP = n_ind*V.dot(np.diag(L**(-2))).dot(V.T)
  for i in range(int(n_perm)):
    if i%1000 == 0: print(i)
    u_y_perm = np.random.permutation(u_y)
    if method == 'CCA':
      x_coef_perm, _, _ = _cca_worker(u_x, u_y_perm, dim_z)
    elif method == 'regression':
      beta_pcs_hat, _, _, _ = np.linalg.lstsq(u_y_perm, u_x)
      x_coef_perm = u_y.dot(beta_pcs_hat)
    B_perm = arr_x.T.dot(x_coef_perm)/np.sqrt(n_ind-1)
    # print(B_perm.T.dot(SP).dot(B_perm))
    perm_var_exp = np.sum(B_perm**2, axis=1)
    #print(abs(arr_xT_coef_perm).sum(), (arr_xT_coef_perm**2).sum())
    perm_mean += perm_var_exp/n_perm
    perm_mean_sq += (perm_var_exp**2)/n_perm
    tail_counts += (perm_var_exp >= true_var_exp)
  return perm_mean, perm_mean_sq, tail_counts


def proj_gene_assoc(df_x, u_x, u_y, dim_z, n_perm=int(1e6), threads=10,
                    alpha=0.05, method='CCA'):
  """Computes association statistics for each gene with the learned projection

  Args:
    df_x: Data matrix.
    u_x: Principal components of arr_x
    u_y: Principal components of arr_y.
    dim_z: Hidden dimension of parent latent space.
    n_perm: Number of permutations to do for significance test.
  Returns
    A df with gene names, z-scores, p-values, true and permuted variance
    explained as well as a column indicating significance under a BHY procesure.
  """

  n_ind, n_genes = df_x.shape
  if method == 'CCA':
    x_coef, _, _ = _cca_worker(u_x, u_y, dim_z)
  elif method == 'regression':
    beta_pcs_hat, _, _, _ = np.linalg.lstsq(u_y, u_x)
    x_coef = u_y.dot(beta_pcs_hat)
  B = df_x.values.T.dot(x_coef)/np.sqrt(n_ind-1)
  true_var_exp = np.sum(B**2, axis=1)
  args = [(df_x.values, u_x, u_y, dim_z, n_perm/threads, true_var_exp, method)
          for _ in range(threads)]
  with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
    res = list(executor.map(_perm_worker, *zip(*args)))
  # Transpose and sum each array we're tracking.
  perm_mean, perm_mean_sq, tail_counts = map(sum, list(map(list, zip(*res))))
  perm_mean = perm_mean/threads
  perm_mean_sq = perm_mean_sq/threads
  perm_std_dev = np.sqrt(perm_mean_sq - perm_mean**2)
  p_vals = tail_counts/n_perm
  z_scores = (true_var_exp - perm_mean)/perm_std_dev
  perm_res = pd.DataFrame(
        {'Z score': z_scores, 'p-value': p_vals,  'Var exp': true_var_exp,
         'Perm mean': perm_mean, 'Perm SD': perm_std_dev}, index=df_x.columns)
  perm_res = perm_res.sort_values('p-value')
  n_genes = perm_res.shape[0]
  cm = 0
  for i in range(1, n_genes+1): cm += 1/i  # There has to be a closed form for this...
  for i, p in enumerate(perm_res['p-value']):
    if p > ((i+1)*alpha)/(n_genes*cm): break
  bhy_cutoff = (i*alpha)/(n_genes*cm)
  print(bhy_cutoff, alpha/n_genes)
  perm_res['significant_bhy'] = (perm_res['p-value'] < bhy_cutoff)
  perm_res['significant_bonf'] = (perm_res['p-value'] < alpha/n_genes)
  return perm_res  


def pca(arr_x, dim=None):
  """Computes first dim principles axes of arr_x.

  Note that this implementation is only equivalent to PCA on a centered
  and scaled data matrix, regardless of whether arr_x is normalized.

  Args:
    arr_x: Data matrix.
    dim: Number of components to compute
  Returns:
    u_x: The principal vectors in column order.
    lambda_x: The loadings vector.
  """
  if dim is None:
    dim = arr_x.shape[0]
  xxt = arr_x.dot(arr_x.T)
  u_x, lambda_x, _ = np.linalg.svd(xxt)
  return u_x[:, 0:dim], np.sqrt(lambda_x[0:dim])


def _cca_worker(u_x, u_y, dim):
  """Performs CCA between u_x and u_y with hidden dimension dim.

  Args:
    u_x: matrix of singular vectors from data x.
    u_y: matrix of singular vectors from data y.
    dim: hidden dimension to use.
  Returns:
    l_x: canonical vectors of u_x.
    l_y: canonical vectors of u_y.
    s_c: singular values of canonical decomposition.
  """
  u_c, rho, v_c = np.linalg.svd(u_x.T.dot(u_y), full_matrices=False)
  l_x = u_x.dot(u_c[:, 0:dim])
  l_y = u_y.dot(v_c.T[:, 0:dim])
  return l_x, l_y, rho


def lr_cca(arr_x, arr_y, dim_x=None, dim_y=None, dim_z=None, standardize=False):
  """Computes (optionally low-rank) CCA between arr_x and arr_y.

  Args:
    arr_x: An n_samples x n_features data matrix.
    arr_y: An n_samples x n_features data matrix.
    dim_x: Keep the highest dim_x PCs of x.
    dim_y: Keep the highest dim_y PCs of y.
    dim_z: Compute dim_z canonical correlations.
    standardize: Whether to standardize the data matrices, not relevant
      if setting dim_x >0 and dim_y > 0.
  Returns:
    xcoef: Estimated coeffcients of arr_x.
    ycoef: Estimated coeffcients of arr_y.
    corr: The correlations
  """
  if standardize:
    arr_x = preprocessing.scale(arr_x)
    arr_y = preprocessing.scale(arr_y)
  if dim_x is not None:
    x_data, _ = pca(arr_x, dim_x)
  else:
    x_data = arr_x
  if dim_y is not None:
    y_data, _ = pca(arr_y, dim_y)
  else:
    y_data = arr_y
  if dim_z is None:
    dim_z = min(x_data.shape[1], y_data.shape[1])
  return _cca_worker(x_data, y_data, dim_z)


def lr_regression(arr_x, arr_y, dim_x=None, dim_y=None, standardize=False):
  if standardize:
    arr_x = preprocessing.scale(arr_x)
    arr_y = preprocessing.scale(arr_y)
  if dim_x is not None:
    x_data, _ = pca(arr_x, dim_x)
  else:
    x_data = arr_x
  if dim_y is not None:
    y_data, _ = pca(arr_y, dim_y)
  else:
    y_data = arr_y
  beta_pcs_hat, _, _, _ = np.linalg.lstsq(y_data, x_data)
  exp_coef = y_data.dot(beta_pcs_hat)
  geno_coef = y_data
  return exp_coef, geno_coef


def project_i_cca(arr_x, u_x, l_x, x_i):
  _arr_x = arr_x.values
  _, s, v = np.linalg.svd(_arr_x, full_matrices=False)
  to_keep = (s > 1e-8)
  s = s[to_keep]
  v = v[to_keep]
  a = x_i.dot((v.T/(s**2)))
  b = v.dot(_arr_x.T).dot(l_x)
  x_i_proj = a.dot(b)
  x_i_full = (x_i_proj.dot(l_x.T)).dot(_arr_x)
  X_full = (l_x.dot(l_x.T)).dot(_arr_x)
  tr_error = np.mean((X_full-_arr_x)**2)
  te_error = np.mean((x_i-x_i_full)**2)
  return x_i_full, tr_error, te_error
