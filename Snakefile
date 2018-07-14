# Snakefile for PCA_CCA
#
# Requires a Python3 environment with snakemake installed. We recommend to do this using
# Conda. For example:
# 
# :$ conda create -n PCA_CCA python=3.6.2
# :$ conda install -n PCA_CCA -c bioconda snakemake
# :$ source activate PCA_CCA
#
# It also requires PLINK1.9 to be installed and available on the commad line. See https://www.cog-genomics.org/plink2


# Configuration
configfile: 'config.yaml'

## Convenience variables
fastq_dir = config['data_dir'] + '/FASTQ'
bed_dir = config['data_dir'] + '/PLINK_1KG'
bed_file = bed_dir + '/1kg_phase1_all.tar.gz'
# Processed from https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/E-GEUV-1.sdrf.txt
sample_info = 'sample_info.txt'
ids_to_process = 'ids_to_process.txt'
abundance_dir = config['working_dir'] + '/transcript_abundance'
genotype_dir = config['working_dir'] + '/genotype'
results_dir = config['working_dir'] + '/results'
kallisto_idx = config['working_dir'] + '/transcript_abundance/transcripts.idx'
working_dir = config['working_dir']

rule results:
  input:
    results_dir + '/1a_uncorrected_pca_pop.png',
    results_dir + '/1b_genotype_pcs.png',
    results_dir + '/1c_uncorrected_pca_lab.png',
    results_dir + '/1d_batch_projection.png',
    results_dir + '/1e_corrected_pca_pop.png',
    results_dir + '/1f_corrected_pca_lab.png',
    results_dir + '/2a_model_exp.png',
    results_dir + '/2b_model_exp_cv_projection.png',
    results_dir + '/3a_p_val_hist.png',
    results_dir + '/3_gene_0_plot.png',
    results_dir + '/3_gene_1_plot.png',
    results_dir + '/3_gene_2_plot.png',
    results_dir + '/projection_associated_genes.tsv'

rule projection_figures:
  input:
    results_dir + '/2a_model_exp.png',
    results_dir + '/2b_model_exp_cv_projection.png'

rule preprocess_data:
  input:
    working_dir + '/expression_corrected.tsv',
    genotype_dir + '/geuvadis.eigenvec'

rule get_data:
  input:
    expand(fastq_dir +  '/{sample}.{string}_1.fastq.gz', zip,
           sample=config['samples'], string=config['samples'].values()), 
    bed_file,
    config['gtf_file'],
    config['transcript_fa']


## Rules for pulling data
rule get_transcript_fa:
  output:
    config['transcript_fa']
  shell:
    ' '.join([
      'wget', '-O', config['transcript_fa'],
      'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_27/gencode.v27.pc_transcripts.fa.gz'])

rule get_gtf:
  output:
    config['gtf_file']
  shell:
    ' '.join([
      'wget', '-O', config['gtf_file'],
      'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_27/gencode.v27.annotation.gtf.gz'])

rule get_sample_fastq:
  output:
    fastq_dir + '/{sample}.{string}_1.fastq.gz',
    fastq_dir + '/{sample}.{string}_2.fastq.gz'
  shell:
    ' '.join([
      'wget', '-P', fastq_dir,
      'ftp://ftp.sra.ebi.ac.uk/vol1/ERA169/ERA169774/fastq/{wildcards.sample}.{wildcards.string}_*.fastq.gz'])

rule get_genotypes:
  output:
    bed_file
  shell:
    ' '.join([
      'wget', '-O', bed_file,
      'https://www.dropbox.com/s/k9ptc4kep9hmvz5/1kg_phase1_all.tar.gz'])

## Kallisto rules
rule kallisto_make_index:
  input:
    config['transcript_fa']
  output:
    kallisto_idx
  shell:
    'kallisto index -i {output} {input}'
 
rule kallisto_quant:
  input:
    idx = kallisto_idx,
    r1 = fastq_dir + '/{sample}.{string}_1.fastq.gz',
    r2 = fastq_dir + '/{sample}.{string}_2.fastq.gz'
  output:
    tsv = abundance_dir + '/{sample}.{string}/abundance.tsv'
  threads: 4
  shell:
    ' '.join(['kallisto quant -i', kallisto_idx, '-o',
              abundance_dir + '/{wildcards.sample}.{wildcards.string}',
              '-t', '{threads}', '{input.r1}', '{input.r2}'])

## Processing rules
rule process_expression:
  input:
    gtf_file = config['gtf_file'],
    ab_tsvs = expand(abundance_dir + '/{sample}.{string}/abundance.tsv', zip,
                     sample=config["samples"], string=config['samples'].values()),
    keep = ids_to_process
  output:
    exp_mat = working_dir + '/expression_normalized.tsv'
  run:
    import pandas as pd
    from pca_cca import util
    # Read each of the abundance tsvs and concatenate them.
    print('Reading kallisto abundances and constructing dataframe.')
    keep = pd.read_csv(input.keep, sep='\t')['sample_name'].values
    exp_trans = []
    for ab_tsv in input.ab_tsvs:
      sample_name = ab_tsv.split('/')[-2].split('.')[0]
      if sample_name in keep:
        ab_df = pd.read_csv(ab_tsv, sep='\t', index_col=0)['tpm']
        ab_df.name = sample_name
        exp_trans.append(ab_df)
    exp_trans = pd.concat(exp_trans, axis=1)
    trans_name, gene_name = zip(*[e.split('|')[0:2] for e in exp_trans.index])
    exp_trans.index = trans_name

    # Drop transcripts that have low mean expression.
    print("Dropping low expressed transcripts and quantile normalizing.")
    print("There were " + str(exp_trans.shape[0]) + " transcripts")
    exp_trans = exp_trans.drop(exp_trans.index[exp_trans.mean(1) < config['min_tpm']])
    print("There are now " + str(exp_trans.shape[0]) + " transcripts")

    # Quantile normalize transcript matrix and transpose
    exp_trans = util.quantile_normalize(exp_trans).T
    
    # Convert transcript levels to gene levels
    print("Converting transcript levels to gene levels.")
    trans_gene_hash = dict(zip(trans_name, gene_name))
    exp_gene = util.trans_to_gene_df(exp_trans, trans_gene_hash)
    # Keep only autosomal non-MHC genes
    print("Dropping non-autosome and MHC genes.")
    genes_to_keep = util.get_genes_to_keep(input.gtf_file)
    exp_gene = exp_gene[exp_gene.columns.intersection(genes_to_keep)]
    print("There are " + str(exp_gene.shape[1]) + " genes remaining.")
    print("Standardizing gene matrix")
    util.standardize_df(exp_gene)
    print("Writing processed gene matrix to file.")
    exp_gene.to_csv(working_dir + '/expression_normalized.tsv', sep='\t')


rule correct_expression:
  input:
    exp_mat = working_dir + '/expression_normalized.tsv',
    sample_info = sample_info
  output:
    exp_mat = working_dir + '/expression_corrected.tsv',
    plot = results_dir + '/1d_batch_projection.png'
  run:
    import numpy as np
    import pandas as pd
    from pca_cca import util
    import matplotlib
    matplotlib.use('Agg')  #workaround for x-windows
    import seaborn as sns
    from IPython import embed
    from sklearn import preprocessing
    print("Loading sample metadata and uncorrected gene expression matrix.")
    sample_info = pd.read_csv(input.sample_info, sep='\t', index_col=0)
    exp_gene = pd.read_csv(input.exp_mat, sep='\t', index_col=0)
    common = sample_info.index.intersection(exp_gene.index)
    sample_info, exp_gene = sample_info.loc[common], exp_gene.loc[common]
    conf_mat = pd.DataFrame(index=sample_info.index)
    conf_mat['sex'] = sample_info['sex'] - 1.0
    conf_mat = conf_mat.join(pd.get_dummies(sample_info['lab']))
    n_conf = conf_mat.shape[1]
    print("Correcting gene expression.")
    if config['correction'] == 'CCA':
      exp_coef, _, _ = util.lr_cca(exp_gene, conf_mat, config['dim_exp'], n_conf, n_conf)
      exp_gene = exp_gene - (exp_coef.dot(exp_coef.T)).dot(exp_gene.values)
    elif config['correction'] == 'regression':
      beta_hat, _, _, _ = np.linalg.lstsq(conf_mat.values, exp_gene)
      exp_gene_hat = conf_mat.values.dot(beta_hat)
      exp_coef, _ = util.pca(exp_gene_hat, 2)
      exp_gene = exp_gene - exp_gene_hat
    exp_gene.loc[:] = preprocessing.scale(exp_gene.values)
    print("Writing corrected gene expression matrix to file.")
    exp_gene.to_csv(working_dir + '/expression_corrected.tsv', sep='\t')

    sample_info['Coordinate 1'] = exp_coef[:, 0]
    sample_info['Coordinate 2'] = exp_coef[:, 1]
    plot_pop = sns.lmplot(
        'Coordinate 1', 'Coordinate 2', hue='lab', data=sample_info, fit_reg=False)
    plot_pop.set(xticks=np.arange(-0.1, 0.15, 0.05))
    plot_pop.savefig(output.plot)


rule unzip_bed_gz:
  input:
    bed_file=bed_file
  output:
    bed = genotype_dir + '/1kg_phase1_all.bed',
    bim = genotype_dir + '/1kg_phase1_all.bim',
    fam = genotype_dir + '/1kg_phase1_all.fam'
  shell:
    'tar -xvf {input.bed_file} -C ' + genotype_dir


rule compute_genotype_pcs:
  input:
    bed = genotype_dir + '/1kg_phase1_all.bed',
    bim = genotype_dir + '/1kg_phase1_all.bim',
    fam = genotype_dir + '/1kg_phase1_all.fam',
    keep = ids_to_process
  output:
    genotype_dir + '/geuvadis.eigenvec',
    bed = genotype_dir + '/geuvadis.bed',
    bim = genotype_dir + '/geuvadis.bim',
    fam = genotype_dir + '/geuvadis.fam',
  shell:
    ' '.join(['bed_with_prefix={input.bed}\n',
              config['plink2_bin'],
              '--bfile ${{bed_with_prefix/.bed/}}',
              '--chr 1-22',
              '--snps-only',
              '--keep {input.keep}',
              '--pca 100 header',
              '--maf', str(config['min_maf']),
              '--make-bed',
              '--out', genotype_dir+'/geuvadis'])


rule model:
  input:
    geno_pcs = genotype_dir + '/geuvadis.eigenvec',
    exp_mat = working_dir + '/expression_corrected.tsv',
    sample_info = sample_info
  output:
    plot = results_dir + '/2a_model_exp.png',
    exp_coef = results_dir + '/exp_coef.txt',
    geno_coef = results_dir + '/geno_coef.txt'
  run:
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')  #workaround for x-windows
    import seaborn as sns
    from pca_cca import util
    exp_mat = pd.read_csv(input.exp_mat, sep='\t', index_col=0)
    geno_pcs = pd.read_csv(
      input.geno_pcs, sep=' ', index_col=0, usecols=range(1,config['dim_snp']+2))
    sample_info = pd.read_csv(input.sample_info, sep='\t', index_col=0)

    common_inds = exp_mat.index.intersection(geno_pcs.index)
    exp_mat = exp_mat.loc[common_inds]
    geno_pcs = geno_pcs.loc[common_inds]
    sample_info = sample_info.loc[common_inds]
    if config['model'] == 'CCA':
        exp_coef, geno_coef, rho = util.lr_cca(exp_mat.values, geno_pcs.values,
             dim_x=config['dim_exp'], dim_z=config['dim_z'])
        print(rho)
    elif config['model'] == 'regression':
        exp_coef, geno_coef = util.lr_regression(exp_mat.values, geno_pcs.values,
             dim_x=config['dim_exp'])
    pd.DataFrame(exp_coef[:,0:config['dim_z']], index=exp_mat.index, 
        columns=['COEF_'+str(i) for i in range(config['dim_z'])]).to_csv(
            output.exp_coef, sep='\t')
    pd.DataFrame(geno_coef[:,0:config['dim_z']], index=geno_pcs.index,
        columns=['COEF_'+str(i) for i in range(config['dim_z'])]).to_csv(
            output.geno_coef, sep='\t')
    sample_info['Coordinate 1'] = exp_coef[:, 0]
    sample_info['Coordinate 2'] = exp_coef[:, 1]
    plot_pop = sns.lmplot(
        'Coordinate 1', 'Coordinate 2', hue='pop', data=sample_info, fit_reg=False)
    plot_pop.savefig(output.plot)


rule compute_loo_genotype_pcs:
  input:
    bed = genotype_dir + '/geuvadis.bed',
    bim = genotype_dir + '/geuvadis.bim',
    fam = genotype_dir + '/geuvadis.fam',
  output:
    genotype_dir + '/loo_pcs/no_{sample}.eigenvec'
  shell:
    ' '.join(['echo 0 {wildcards.sample} >',
              genotype_dir+'/loo_pcs/remove_{wildcards.sample}.txt\n',
              'bed_with_prefix={input.bed}\n',
              config['plink2_bin'],
              '--bfile ${{bed_with_prefix/.bed/}}',
              '--remove', genotype_dir+'/loo_pcs/remove_{wildcards.sample}.txt',
              '--pca 100 header',
              '--maf', str(config['min_maf']),
              '--snps-only',
              '--out', genotype_dir+'/loo_pcs/no_{wildcards.sample}'])


rule find_significant_genes:
  input:
    exp_mat = working_dir + '/expression_corrected.tsv',
    geno_pcs = genotype_dir + '/geuvadis.eigenvec'
  output:
    genes = results_dir + '/projection_associated_genes.tsv'
  threads: 10
  run:
    import pandas as pd
    from pca_cca import util
    if config['model'] == 'regression':
        raise ValueError('Gene finding not supported for regression model.')    
    exp_mat = pd.read_csv(input.exp_mat, index_col=0, sep='\t')
    geno_pcs = pd.read_csv(
      input.geno_pcs, sep=' ', index_col=0, usecols=range(1,config['dim_snp']+2))
    u_exp, _ = util.pca(exp_mat)
    perm_res = util.proj_gene_assoc(
        exp_mat, u_exp, geno_pcs, config['dim_z'], n_perm=config['n_perm'], threads=threads)
    perm_res = perm_res.sort_values(by=['p-value', 'Z score'], ascending = [True, False])
    perm_res.to_csv(output.genes, sep='\t')


rule plot_significant_genes:
  input:
    exp_gene = working_dir + '/expression_corrected.tsv',
    genes_list = results_dir + '/projection_associated_genes.tsv',
    sample_info = sample_info
  output:
    results_dir + '/3a_p_val_hist.png',
    results_dir + '/3_gene_0_plot.png',
    results_dir + '/3_gene_1_plot.png',
    results_dir + '/3_gene_2_plot.png'
  run:
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')  #workaround for x-windows
    import seaborn as sns
    from pca_cca import util
    import matplotlib.pyplot as plt
    n_genes_to_plot = 10
    exp_gene = pd.read_csv(input.exp_gene, index_col=0, sep='\t')
    perm_res = pd.read_csv(input.genes_list, index_col=0, sep='\t')
    sample_info = pd.read_csv(input.sample_info, sep='\t', index_col=0)
    common = sample_info.index.intersection(exp_gene.index)
    sample_info, exp_gene = sample_info.loc[common], exp_gene.loc[common]
    exp_gene['pop'] = sample_info['pop']
    p_val_hist = sns.distplot(perm_res['p-value'], kde=False)
    fig = p_val_hist.get_figure()
    fig.savefig(results_dir + '/3a_p_val_hist.png')
    plt.clf()
    perm_res = perm_res.sort_values(by=['p-value', 'Z score'], ascending=[True, False])
    for i in range(n_genes_to_plot):
      gene_name = perm_res.index[i]
      gene_plot = sns.violinplot(x=gene_name, y='pop', data=exp_gene, inner="point")
      fig = gene_plot.get_figure()
      fig.savefig(results_dir + '/3_gene_' + str(i) + '_plot.png')
      plt.clf()


rule plot_uncorrected_expression:
  input:
    exp_mat = working_dir + '/expression_normalized.tsv',
    sample_info = sample_info
  output:
    results_dir + '/1a_uncorrected_pca_pop.png',
    results_dir + '/1c_uncorrected_pca_lab.png' 
  run:
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')  #workaround for x-windows
    import seaborn as sns
    from pca_cca import util
    sample_info = pd.read_csv(input.sample_info, sep='\t', index_col=0)
    exp_gene = pd.read_csv(input.exp_mat, sep='\t', index_col=0)
    common = sample_info.index.intersection(exp_gene.index)
    sample_info, exp_gene = sample_info.loc[common], exp_gene.loc[common]
    u_gene, _ = util.pca(exp_gene, 2)
    sample_info['PC1'] = u_gene[:, 0]
    sample_info['PC2'] = u_gene[:, 1]
    plot_lab = sns.lmplot('PC1', 'PC2', hue='lab', data=sample_info, fit_reg=False)
    plot_lab.savefig(results_dir + '/1c_uncorrected_pca_lab.png')
    plot_pop = sns.lmplot('PC1', 'PC2', hue='pop', data=sample_info, fit_reg=False)
    plot_pop.savefig(results_dir + '/1a_uncorrected_pca_pop.png')


rule plot_corrected_expression:
  input:
    exp_mat = working_dir + '/expression_corrected.tsv',
    sample_info = sample_info
  output:
    results_dir + '/1e_corrected_pca_pop.png',
    results_dir + '/1f_corrected_pca_lab.png' 
  run:
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')  #workaround for x-windows
    import seaborn as sns
    from pca_cca import util
    sample_info = pd.read_csv(input.sample_info, sep='\t', index_col=0)
    exp_gene = pd.read_csv(input.exp_mat, sep='\t', index_col=0)
    common = sample_info.index.intersection(exp_gene.index)
    sample_info, exp_gene = sample_info.loc[common], exp_gene.loc[common]
    u_gene, _ = util.pca(exp_gene, 2)
    sample_info['PC1'] = u_gene[:, 0]
    sample_info['PC2'] = u_gene[:, 1]
    plot_lab = sns.lmplot('PC1', 'PC2', hue='lab', data=sample_info, fit_reg=False)
    plot_lab.set(xticks=np.arange(-0.2, 0.25, 0.05), yticks=np.arange(-0.2, 0.25, 0.05))
    plot_lab.savefig(results_dir + '/1f_corrected_pca_lab.png')
    plot_pop = sns.lmplot('PC1', 'PC2', hue='pop', data=sample_info, fit_reg=False)
    plot_pop.set(xticks=np.arange(-0.2, 0.25, 0.05), yticks=np.arange(-0.2, 0.25, 0.05))
    plot_pop.savefig(results_dir + '/1e_corrected_pca_pop.png')


rule plot_genotype_pcs:
  input:
    geno_pcs = genotype_dir + '/geuvadis.eigenvec',
    sample_info = sample_info
  output:
    results_dir + '/1b_genotype_pcs.png'
  run:
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')  #workaround for x-windows
    import seaborn as sns
    sample_info = pd.read_csv(input.sample_info, sep='\t', index_col=0)
    geno_pcs = pd.read_csv(
        input.geno_pcs, sep=' ', index_col=0, usecols=range(1,config['dim_snp']+2))
    sample_info = sample_info.loc[sample_info.index.intersection(geno_pcs.index)]
    sample_info['PC1'] = geno_pcs['PC1']
    sample_info['PC2'] = geno_pcs['PC2']
    plot_pop = sns.lmplot('PC1', 'PC2', hue='pop', data=sample_info, fit_reg=False)
    plot_pop.savefig(results_dir + '/1b_genotype_pcs.png')


rule make_cv_projection_plot:
  input:
    exp_mat = working_dir + '/expression_corrected.tsv',
    geno_pcs = expand(genotype_dir + '/loo_pcs/no_{sample}.eigenvec', sample=config['samples']),
    sample_info = sample_info
  output:
    results_dir + '/2b_model_exp_cv_projection.png'
  threads: 8
  run:
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')  #workaround for x-windows
    import seaborn as sns
    from pca_cca import util
    from IPython import embed
    if config['model'] == 'regression':
        raise ValueError('CV not supported for regression model.')
    exp_mat = pd.read_csv(input.exp_mat, index_col=0, sep='\t')
    sample_info = pd.read_csv(input.sample_info, sep='\t', index_col=0)
    common_inds = exp_mat.index.intersection(sample_info.index)
    exp_mat = exp_mat.loc[common_inds]
    sample_info = sample_info.loc[common_inds]
    cv_exp_mat, error_df = util.build_cv_df(
      exp_mat, dict(zip(config['samples'], input.geno_pcs)),
        config['dim_exp'], config['dim_snp'], config['dim_z'], threads)
    u_cv_exp, _ = util.pca(cv_exp_mat.values, dim=2)
    sample_info['Coordinate 1'] = u_cv_exp[:, 0]
    sample_info['Coordinate 2'] = u_cv_exp[:, 1]
    plot_pop = sns.lmplot(
        'Coordinate 1', 'Coordinate 2', hue='pop', data=sample_info, fit_reg=False)
    plot_pop.savefig(results_dir + '/2b_model_exp_cv_projection.png')
    error_df.to_csv(results_dir + '/cv_error_result.txt', sep='\t')
