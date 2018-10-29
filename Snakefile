# Snakefile for PCCA
#
# See README.md for more information.

# Configuration
configfile: 'config.yaml'

# Convenience variables
fastq_dir = config['data_dir'] + '/FASTQ'
bed_dir = config['data_dir'] + '/PLINK_1KG'
bed_file = bed_dir + '/1kg_phase1_all.tar.gz'
# Processed from https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/E-GEUV-1.sdrf.txt
sample_info = 'sample_info.txt'
ids_to_process = 'ids_to_process.txt'
abundance_dir = config['working_dir'] + '/transcript_abundance'
genotype_dir = config['working_dir'] + '/genotype'
results_dir = config['working_dir'] + '/results'
eqtl_dir = config['working_dir'] + '/eqtl_analysis'
kallisto_idx = config['working_dir'] + '/transcript_abundance/transcripts.idx'
working_dir = config['working_dir']

eqtl_window_size = 500000

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
    results_dir + '/projection_associated_genes.tsv',
    results_dir + '/venn_geuvadis_results_comparison.png',
    results_dir + '/venn_geuvadis_diff_pcca_genes.txt',
    results_dir + '/venn_pcca_diff_geuvadis_genes.txt',
    results_dir + '/model_no_eqtl.png'


rule eqtl_analysis:
  input:
    eqtl_dir + '/pca_results_all.txt',
    eqtl_dir + '/pca_results_top.txt',
    eqtl_dir + '/pcca_results_all.txt',
    eqtl_dir + '/pcca_results_top.txt',
    results_dir + '/qq_plot.png',
    results_dir + '/covar_genes_comparison.png',
    results_dir + '/venn_eqtl_comparison.png',
    results_dir + '/venn_pca_diff_pcca_genes.txt',
    results_dir + '/venn_pcca_diff_pca_genes.txt'

rule projection_figures:
  input:
    results_dir + '/2a_model_exp.png',
    results_dir + '/2b_model_exp_cv_projection.png'

rule svplots:
  input:
    results_dir + '/S_exp_svs.png',
    results_dir + '/S_geno_svs.png',

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
    config['transcript_fa'],
    config['data_dir'] + 'EUR373.gene.cis.FDR5.best.rs137.txt.gz',
    config['data_dir'] + 'YRI89.gene.cis.FDR5.best.rs137.txt.gz'


## Rules for pulling data
rule get_transcript_fa:
  output:
    config['transcript_fa']
  shell:
    ' '.join([
      'wget', '-O', config['transcript_fa'],
      'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_27/gencode.v27.pc_transcripts.fa.gz'])

rule get_geuvadis_results:
  output:
    eur_results=config['data_dir'] + '/EUR373.gene.cis.FDR5.best.rs137.txt.gz',
    yri_results=config['data_dir'] + '/YRI89.gene.cis.FDR5.best.rs137.txt.gz'
  run:
    command = ' '.join([
      'wget', '-O', output.eur_results,
      'https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/EUR373.gene.cis.FDR5.best.rs137.txt.gz'])
    shell(command)
    command = ' '.join([
      'wget', '-O', output.yri_results,
      'https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/YRI89.gene.cis.FDR5.best.rs137.txt.gz'])
    shell(command)

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
    from pcca import util
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
    exp_plink = working_dir + '/expression_plink.tsv',
    plot = results_dir + '/1d_batch_projection.png'
  run:
    import numpy as np
    import pandas as pd
    from pcca import util
    import matplotlib
    matplotlib.use('Agg')  #workaround for x-windows
    import seaborn as sns
    from IPython import embed
    from sklearn import preprocessing
    print("Loading sample metadata and uncorrected gene expression matrix.")
    print(config['correction'])
    sample_info = pd.read_csv(input.sample_info, sep='\t', index_col=0)
    exp_gene = pd.read_csv(input.exp_mat, sep='\t', index_col=0)
    common = sample_info.index.intersection(exp_gene.index)
    sample_info, exp_gene = sample_info.loc[common], exp_gene.loc[common]
    conf_mat = pd.DataFrame(index=sample_info.index)
    conf_mat['sex'] = sample_info['sex'] - 1.0
    conf_mat = conf_mat.join(pd.get_dummies(sample_info['lab']))
    n_conf = conf_mat.shape[1]
    print("Correcting gene expression.")
    if config['correction'][0:4] == 'peer':
      gene_cols = exp_gene.columns
      gene_index = exp_gene.index
      exp_gene.to_csv(
          working_dir + '/expression_normalized_raw.csv', index=False, header=False)
      command = ' '.join([
          config['peertool_bin'],
          '-n', config['correction'][4:],
          '-f', working_dir + '/expression_normalized_raw.csv'])
      shell(command)
      exp_gene = pd.read_csv(
          working_dir + '/peer_out/residuals.csv', index_col=False, header=None).T
      exp_gene.index = gene_index
      exp_gene.columns = gene_cols
      exp_coef = pd.read_csv(
          working_dir + '/peer_out/X.csv', index_col=False, header=None).T.values
    elif config['correction'] == 'CCA':
      exp_coef, _, _ = util.lr_cca(exp_gene, conf_mat, config['dim_exp'], n_conf, n_conf)
      exp_gene = exp_gene - (exp_coef.dot(exp_coef.T)).dot(exp_gene.values)
    elif config['correction'] == 'regression':
      beta_hat, _, _, _ = np.linalg.lstsq(conf_mat.values, exp_gene)
      exp_gene_hat = conf_mat.values.dot(beta_hat)
      exp_coef, _ = util.pca(exp_gene_hat, 2)
      exp_gene = exp_gene - exp_gene_hat
    elif config['correction'] == 'None':
       exp_coef = np.zeros((exp_gene.shape[0], 2))
    else:
      raise ValueError('Correction must be CCA, regression, peerN or None.')
    exp_gene.loc[:] = preprocessing.scale(exp_gene.values)
    print("Writing corrected gene expression matrix to file.")
    exp_gene.to_csv(working_dir + '/expression_corrected.tsv', sep='\t')
    exp_gene.insert(0, 'FID', 0)
    exp_gene.insert(1, 'IID', exp_gene.index)
    exp_gene.to_csv(working_dir + '/expression_plink.tsv', sep='\t', index=False)

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
    genotype_dir + '/geuvadis.eigenval',
    bed = genotype_dir + '/geuvadis.bed',
    bim = genotype_dir + '/geuvadis.bim',
    fam = genotype_dir + '/geuvadis.fam',
  run:
    command = ' '.join(
      ['bed_with_prefix={input.bed}\n',
       config['plink2_bin'],
       '--bfile ${{bed_with_prefix/.bed/}}',
       '--chr 1-22',
       '--snps-only',
       '--keep {input.keep}',
       '--pca 999 header',  # TODO(brielin): infer number of inds
       '--maf', str(config['min_maf']),
       '--make-bed',
       '--out', genotype_dir+'/geuvadis'])
    shell(command)
    if config['subsample_geno'] < 1.0:
      command = ' '.join([
        'bed_with_prefix={output.bed}\n',
        config['plink2_bin'],
        '--bfile ${{bed_with_prefix/.bed/}}',
        '--thin', str(config['subsample_geno']),
        '--pca 999 header',  # TODO(brielin): infer number of inds
        '--make-bed',
        '--out', genotype_dir+'/geuvadis'])
      shell(command)

rule plot_svs:
  input:
    geno_evs = genotype_dir + '/geuvadis.eigenval',
    exp_mat = working_dir + '/expression_corrected.tsv'
  output:
    plot_exp = results_dir + '/S_pca_ve_exp.png',
    plot_geno = results_dir + '/S_pca_ve_geno.png'
  run:
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')  #workaround for x-windows
    import matplotlib.pyplot as plt
    import seaborn as sns

    exp_mat = pd.read_csv(input.exp_mat, sep='\t', index_col=0)
    svs = pd.read_csv(input.geno_evs, sep=' ', header=None, names=['geno_val'])
    # sum_geno_vals = [np.sqrt(svs.loc[0:i, 'geno_val']).sum() for i in range(len(svs))]
    sum_geno_vals = [svs.loc[0:i, 'geno_val'].sum() for i in range(len(svs))]
    svs['sum_geno_vals'] = sum_geno_vals
    svs['genotype'] = svs['sum_geno_vals']/np.sqrt(svs['geno_val']).sum()

    u, s, v = np.linalg.svd(exp_mat)
    svs['exp_val'] = s**2  # TODO(brielin): SV vs EV?
    sum_exp_vals = [svs.loc[0:i, 'exp_val'].sum() for i in range(len(svs))]
    svs['sum_exp_vals'] = sum_exp_vals
    svs['expression'] = svs['sum_exp_vals']/svs['exp_val'].sum()
    svs['component'] = svs.index + 1
    # svs_long = pd.melt(
    #     svs[['genotype', 'expression', 'component']],
    #     id_vars=['component'],
    #     value_vars=['genotype', 'expression'],
    #     value_name='Percent',
    #     var_name='Data type')
    plot = sns.lineplot(y='genotype', x='component', data=svs.iloc[0:15])
    plt.title('Genotype')
    plt.xlabel('Number of components')
    plt.ylabel('Percent variance explained')
    plot.figure.savefig(output.plot_geno)
    plt.clf()
    plot = sns.lineplot(y='expression', x='component', data=svs.iloc[0:100])
    plt.title('Gene expression')
    plt.xlabel('Number of components')
    plt.ylabel('Percent variance explained')
    plot.figure.savefig(output.plot_exp)


rule model:
  input:
    geno_pcs = genotype_dir + '/geuvadis.eigenvec',
    exp_mat = working_dir + '/expression_corrected.tsv',
    sample_info = sample_info
  output:
    plot = results_dir + '/2a_model_exp.png',
    exp_coef = results_dir + '/exp_coef.txt',
    geno_coef = results_dir + '/geno_coef.txt',
    covar_pcca = working_dir + '/covar_pcca.txt',
    covar_pca = working_dir + '/covar_pca.txt'
  run:
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')  #workaround for x-windows
    import seaborn as sns
    from pcca import util
    exp_mat = pd.read_csv(input.exp_mat, sep='\t', index_col=0)
    keep_gene = np.random.binomial(n=1, p=config['subsample_exp'], size=len(exp_mat.columns))
    exp_mat = exp_mat.loc[:, exp_mat.columns[keep_gene > 0]]
    geno_pcs = pd.read_csv(
      input.geno_pcs, sep=' ', index_col=0, usecols=range(1, config['dim_snp']+2))
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
    exp_coef = pd.DataFrame(exp_coef[:,0:config['dim_z']], index=exp_mat.index, 
        columns=['exp_coef_'+str(i) for i in range(config['dim_z'])])
    exp_coef.to_csv(output.exp_coef, sep='\t')
    geno_coef = pd.DataFrame(geno_coef[:,0:config['dim_z']], index=geno_pcs.index,
        columns=['geno_coef_'+str(i) for i in range(config['dim_z'])])
    geno_coef.to_csv(output.geno_coef, sep='\t')
    covar = sample_info['sex'].to_frame()
    covar['FID'] = 0
    covar = covar.join(geno_pcs)
    covar = covar.join(exp_coef)
    covar = covar.join(geno_coef)
    covar.reset_index(inplace=True)
    covar[['FID', 'IID', 'sex'] + list(geno_pcs.columns)[0:5]].to_csv(
        output.covar_pca, sep='\t', index=False)
    covar[['FID', 'IID', 'sex'] + list(geno_coef.columns) + list(exp_coef.columns)].to_csv(
        output.covar_pcca, sep='\t', index=False)
    sample_info['Coordinate 1'] = exp_coef.values[:, 0]
    sample_info['Coordinate 2'] = exp_coef.values[:, 1]
    plot_pop = sns.lmplot(
        'Coordinate 1', 'Coordinate 2', hue='pop', data=sample_info, fit_reg=False)
    plot_pop.savefig(output.plot)

rule find_signficant_pc_genes:
  input:
    exp_mat = working_dir + '/expression_corrected.tsv',
  output:
    genes = results_dir + '/pca_associated_genes.tsv'
  threads: 10
  run:
    import pandas as pd
    from pcca import util
    exp_mat = pd.read_csv(input.exp_mat, index_col=0, sep='\t')
    perm_res = util.pca_gene_assoc(exp_mat, 2,
                                   n_perm=config['n_perm'], threads=threads)
    perm_res.to_csv(output.genes, sep='\t')

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
              '--pca 999 header',
              '--maf', str(config['min_maf']),
              '--snps-only',
              '--out', genotype_dir+'/loo_pcs/no_{wildcards.sample}'])


rule find_significant_genes:
  input:
    exp_mat = working_dir + '/expression_corrected.tsv',
    geno_pcs = genotype_dir + '/geuvadis.eigenvec',
    bed = genotype_dir + '/geuvadis.bed'
  output:
    genes = results_dir + '/projection_associated_genes.tsv'
  threads: 10
  run:
    import pandas as pd
    from pcca import util
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


rule run_eqtl_analysis:
  input:
    exp_mat = working_dir + '/expression_plink.tsv',
    covar_pcca = working_dir + '/covar_pcca.txt',
    covar_pca = working_dir + '/covar_pca.txt',
    gtf_file = config['gtf_file']
  output:
    eqtl_dir + '/pca_results_all.txt',
    eqtl_dir + '/pca_results_top.txt',
    eqtl_dir + '/pcca_results_all.txt',
    eqtl_dir + '/pcca_results_top.txt',
  threads: 30
  run:
    from pcca import util
    from glob import glob
    import pandas as pd
    import concurrent.futures
    import matplotlib
    matplotlib.use('Agg')  #workaround for x-windows
    import matplotlib.pyplot as plt
    import seaborn as sns

    genes_to_use = util.get_col_names(input.exp_mat)
    gene_location_dict = util.build_gene_location_dict(input.gtf_file, genes_to_use)
    command_list_pca = []
    command_list_pcca = []
    for gene_name, (chrom, tss) in gene_location_dict.items():
        command = ' '.join([
              config['plink2_bin'],
              '--pheno expression_plink.tsv',
              '--pheno-name', gene_name,
              '--bfile', genotype_dir + '/geuvadis',
              '--maf', str(config['min_maf']),
              '--chr', chrom,
              '--from-bp', str(max(int(tss) - eqtl_window_size, 0)),
              '--to-bp', str(int(tss) + eqtl_window_size),
              '--linear'])
        command_pca = ' '.join([
              command, '--covar', input.covar_pca,
              '--out', eqtl_dir + '/' + gene_name + '_pca'])
        command_pcca = ' '.join([
              command, '--covar', input.covar_pcca,
              '--out', eqtl_dir + '/' + gene_name + '_pcca'])
        command_list_pca.append(command_pca)
        command_list_pcca.append(command_pcca)
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
          executor.map(shell, command_list_pca + command_list_pcca)
    pca_results_all = []
    pca_results_top = []
    for pca_result_file in glob(eqtl_dir + '/*_pca.assoc.linear'):
      gene = pca_result_file.split('/')[-1].split('_')[0]
      result = pd.read_table(pca_result_file, sep='\s+')
      result = result[result['TEST'] == 'ADD'].sort_values('P')
      result['gene'] = gene
      pca_results_all.append(result)
      pca_results_top.append(result.iloc[0,:].to_frame().T)
    pcca_results_all = []
    pcca_results_top = []
    for pcca_result_file in glob(eqtl_dir + '/*_pcca.assoc.linear'):
      gene = pcca_result_file.split('/')[-1].split('_')[0]
      result = pd.read_table(pcca_result_file, sep='\s+')
      result = result[result['TEST'] == 'ADD'].sort_values('P')
      result['gene'] = gene
      pcca_results_all.append(result)
      pcca_results_top.append(result.iloc[0,:].to_frame().T)
    pcca_results_all = pd.concat(pcca_results_all, ignore_index=True)
    pcca_results_top = pd.concat(pcca_results_top, ignore_index=True).sort_values('P')
    pcca_results_all.dropna(inplace=True)
    pca_results_all = pd.concat(pca_results_all, ignore_index=True)
    pca_results_top = pd.concat(pca_results_top, ignore_index=True).sort_values('P')
    pca_results_all.dropna(inplace=True)
    pca_results_all.to_csv(eqtl_dir + '/pca_results_all.txt', sep='\t')
    pca_results_top.to_csv(eqtl_dir + '/pca_results_top.txt', sep='\t')
    pcca_results_all.to_csv(eqtl_dir + '/pcca_results_all.txt', sep='\t')
    pcca_results_top.to_csv(eqtl_dir + '/pcca_results_top.txt', sep='\t')


rule make_eqtl_plots:
  input:
    eqtl_dir + '/pca_results_all.txt',
    eqtl_dir + '/pca_results_top.txt',
    eqtl_dir + '/pcca_results_all.txt',
    eqtl_dir + '/pcca_results_top.txt',
  output:
    results_dir + '/qq_plot.png',
    results_dir + '/covar_genes_comparison.png',
    results_dir + '/venn_eqtl_comparison.png',
    results_dir + '/venn_pca_diff_pcca_genes.txt',
    results_dir + '/venn_pcca_diff_pca_genes.txt'
  run:
    import pandas as pd
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')  #workaround for x-windows
    import seaborn as sns
    import matplotlib.pyplot as plt
    from matplotlib_venn import venn2
    pca_results_all = pd.read_csv(eqtl_dir + '/pca_results_all.txt', sep='\t')
    pca_results_top = pd.read_csv(eqtl_dir + '/pca_results_top.txt', sep='\t')
    pcca_results_all = pd.read_csv(eqtl_dir + '/pcca_results_all.txt', sep='\t')
    pcca_results_top = pd.read_csv(eqtl_dir + '/pcca_results_top.txt', sep='\t')

    qe_pca, qo_pca = util.make_qq_plot_data(pca_results_all['P'].values, points=1000)
    qe_pcca, qo_pcca = util.make_qq_plot_data(
        pcca_results_all['P'].values, points=1000)
    plot_data = pd.DataFrame({
        'expected': np.concatenate([qe_pca, qe_pcca]),
        'observed': np.concatenate([qo_pca, qo_pcca]),
        'method': ['pca']*len(qe_pca) + ['pcca']*len(qe_pcca)})

    ax = sns.lmplot(x='expected', y='observed', data=plot_data, hue='method', fit_reg=False)
    plt.plot(qe_pca, qe_pca, color="k", ls="--")
    plt.savefig(results_dir + '/qq_plot.png')
    plt.clf()

    pv_cutoffs = [1/10**i for i in range(4,10)]
    genes_by_cutoff = pd.DataFrame({
        'n_genes': [(pca_results_top['P'] < cutoff).sum() for cutoff in pv_cutoffs] +
               [(pcca_results_top['P'] < cutoff).sum() for cutoff in pv_cutoffs],
        'method': ['pca']*len(pv_cutoffs) + ['pcca']*len(pv_cutoffs),
        'cutoff': pv_cutoffs*2})
    plot = sns.lineplot('cutoff', 'n_genes', data=genes_by_cutoff, hue='method')
    plt.xscale('log')
    plot.figure.savefig(results_dir + '/covar_genes_comparison.png')
    plt.clf()

    cutoff = float(config['eqtl_threshold'])
    pca_genes = set(pca_results_top.loc[pca_results_top['P'] < cutoff]['gene'].values)
    pcca_genes = set(pcca_results_top.loc[pcca_results_top['P'] < cutoff]['gene'].values)
    venn2([pca_genes, pcca_genes], set_labels=['PCA eQTL genes', 'PCA CCA eQTL genes'])
    plt.savefig(results_dir + '/venn_eqtl_comparison.png')
    plt.clf()
    pd.Series(list(pca_genes.difference(pcca_genes))).to_csv(results_dir + '/venn_pca_diff_pcca_genes.txt')
    pd.Series(list(pcca_genes.difference(pca_genes))).to_csv(results_dir + '/venn_pcca_diff_pca_genes.txt')
    pd.Series(list(pca_genes.intersection(pcca_genes))).to_csv(results_dir + '/venn_intersection_genes.txt')


rule geuvadis_comparison:
  input:
    config['data_dir'] + '/YRI89.gene.cis.FDR5.best.rs137.txt.gz',
    config['data_dir'] + '/EUR373.gene.cis.FDR5.best.rs137.txt.gz',
    results_dir + '/projection_associated_genes.tsv'
  output:
    results_dir + '/venn_geuvadis_results_comparison.png',
    results_dir + '/venn_geuvadis_diff_pcca_genes.txt',
    results_dir + '/venn_pcca_diff_geuvadis_genes.txt'
  run:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import pandas as pd
    from matplotlib_venn import venn2
    geu_res_eur = pd.read_csv(
        config['data_dir'] + '/EUR373.gene.cis.FDR5.best.rs137.txt.gz', sep='\t', header=None)
    geu_res_yri = pd.read_csv(
        config['data_dir'] + '/YRI89.gene.cis.FDR5.best.rs137.txt.gz', sep='\t', header=None)
    pcca_gene_table = pd.read_csv(results_dir + '/projection_associated_genes.tsv', sep='\t', index_col=0)
    eur_genes = set([gene.split('.')[0] for gene in geu_res_eur.iloc[:, 2].values])
    yri_genes = set([gene.split('.')[0] for gene in geu_res_yri.iloc[:, 2].values])
    geuvadis_genes = eur_genes.union(yri_genes)
    pcca_genes = pcca_gene_table.loc[pcca_gene_table['significant_bhy'] == True].index.values
    pcca_genes = set(gene.split('.')[0] for gene in pcca_genes)
    venn2([eur_genes.union(yri_genes), pcca_genes], set_labels=['geuvadis eQTL genes', 'PCA CCA projection genes'])
    plt.savefig(results_dir + '/venn_geuvadis_results_comparison.png')
    plt.clf()
    pd.Series(list(geuvadis_genes.difference(pcca_genes))).to_csv(results_dir + '/venn_geuvadis_diff_pcca_genes.txt')
    pd.Series(list(pcca_genes.difference(geuvadis_genes))).to_csv(results_dir + '/venn_pcca_diff_geuvadis_genes.txt')



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
    from pcca import util
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
    from pcca import util
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
    from pcca import util
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
    from pcca import util
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


rule compute_af_by_pop:
  input:
    bed = genotype_dir + '/geuvadis.bed'
  output:
    genotype_dir + '/gbr_af.frq',
    genotype_dir + '/fin_af.frq',
    genotype_dir + '/tsi_af.frq',
    genotype_dir + '/yri_af.frq',
    genotype_dir + '/all_af.frq'
  shell:
    """
    awk 'BEGIN{{print "fid\tsample_name"}}; $4=="GBR" {{print "0\t"$1}}' \
      sample_info.txt > gbr_ids.txt
    awk 'BEGIN{{print "fid\tsample_name"}}; $4=="FIN" {{print "0\t"$1}}' \
      sample_info.txt > fin_ids.txt
    awk 'BEGIN{{print "fid\tsample_name"}}; $4=="TSI" {{print "0\t"$1}}' \
      sample_info.txt > tsi_ids.txt
    awk 'BEGIN{{print "fid\tsample_name"}}; $4=="YRI" {{print "0\t"$1}}' \
      sample_info.txt > yri_ids.txt
    bed_with_prefix={input.bed}
    {config[plink2_bin]} --bfile ${{bed_with_prefix/.bed/}} \
      --freq --out {genotype_dir}/all_af
    {config[plink2_bin]} --bfile ${{bed_with_prefix/.bed/}} --keep gbr_ids.txt \
      --freq --out {genotype_dir}/gbr_af
    {config[plink2_bin]} --bfile ${{bed_with_prefix/.bed/}} --keep fin_ids.txt \
      --freq --out {genotype_dir}/fin_af
    {config[plink2_bin]} --bfile ${{bed_with_prefix/.bed/}} --keep tsi_ids.txt \
      --freq --out {genotype_dir}/tsi_af
    {config[plink2_bin]} --bfile ${{bed_with_prefix/.bed/}} --keep yri_ids.txt \
      --freq --out {genotype_dir}/yri_af
    """


rule model_no_eqtl:
  input:
    eqtl_eur = config['data_dir'] + '/EUR373.gene.cis.FDR5.best.rs137.txt.gz',
    eqtl_yri = config['data_dir'] + '/YRI89.gene.cis.FDR5.best.rs137.txt.gz',
    af_gbr = genotype_dir + '/gbr_af.frq',
    af_fin = genotype_dir + '/fin_af.frq',
    af_tsi = genotype_dir + '/tsi_af.frq',
    af_yri = genotype_dir + '/yri_af.frq',
    af_all = genotype_dir + '/all_af.frq',
    okg_bim = genotype_dir + '/1kg_phase1_all.bim',
    exp_mat = working_dir + '/expression_corrected.tsv',
    sample_info = working_dir + '/sample_info.txt',
    geno_pcs = genotype_dir + '/geuvadis.eigenvec'
  output:
    plot = results_dir + '/model_no_eqtl.png'
  run:
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')  #workaround for x-windows
    import seaborn as sns
    from pcca import util
    from sklearn import preprocessing
    exp_mat = pd.read_csv(input.exp_mat, sep='\t', index_col=0)
    sample_info = pd.read_csv('sample_info.txt', sep='\t', index_col=0)
    # A2 is ref in original file
    ref_alleles = pd.read_csv('genotype/1kg_phase1_all.bim', usecols=[1, 4, 5], sep='\s+',
                              names=['rsid', 'alt', 'ref'], index_col=0)
    # A2 is MAJOR after processing
    af_all = pd.read_csv(input.af_all, sep='\s+', index_col=1)
    af_gbr = pd.read_csv(input.af_gbr, sep='\s+', index_col=1)
    af_fin = pd.read_csv(input.af_fin, sep='\s+', index_col=1)
    af_tsi = pd.read_csv(input.af_tsi, sep='\s+', index_col=1)
    af_yri = pd.read_csv(input.af_yri, sep='\s+', index_col=1)

    eqtl_eur = pd.read_csv(input.eqtl_eur, sep='\t', header=None, usecols=[0,2,9],
                           names=['rsid', 'gene', 'corr'])
    eqtl_yri = pd.read_csv(input.eqtl_yri, sep='\t', header=None, usecols=[0,2,9],
                           names=['rsid', 'gene', 'corr'])
    geno_pcs = pd.read_csv(
      input.geno_pcs, sep=' ', index_col=0, usecols=range(1, config['dim_snp']+2))
    common_inds = exp_mat.index.intersection(geno_pcs.index)
    exp_mat = exp_mat.loc[common_inds]
    geno_pcs = geno_pcs.loc[common_inds]
    sample_info = sample_info.loc[common_inds]

    eqtl_eur.index = [gene.split('.')[0] for gene in eqtl_eur['gene']]
    eqtl_yri.index = [gene.split('.')[0] for gene in eqtl_yri['gene']]
    
    eqtls = pd.concat([eqtl_eur, eqtl_yri.loc[eqtl_yri.index.difference(eqtl_eur.index)]])
    for gene in eqtl_yri.index.intersection(eqtl_eur.index):
      corr_yri = eqtl_yri.loc[gene, 'corr']
      corr_eur = eqtl_eur.loc[gene, 'corr']
      if abs(corr_yri) > abs(corr_eur):
        eqtls.loc[gene, ['rsid', 'corr']] = eqtl_yri.loc[gene, ['rsid', 'corr']]

    index = af_all.index.intersection(eqtls.rsid)
    ref_alleles = ref_alleles.loc[index]
    af_all = af_all.loc[index]
    af_gbr = af_gbr.loc[index]
    af_fin = af_fin.loc[index]
    af_tsi = af_tsi.loc[index]
    af_yri = af_yri.loc[index]

    ref_alleles['af_alt_all'] = [af if a1 == alt else 1-af for af, a1, alt in zip(af_all['MAF'], af_all['A1'], ref_alleles['alt'])]
    ref_alleles['af_alt_gbr'] = [af if a1 == alt else 1-af for af, a1, alt in zip(af_gbr['MAF'], af_gbr['A1'], ref_alleles['alt'])]
    ref_alleles['af_alt_fin'] = [af if a1 == alt else 1-af for af, a1, alt in zip(af_fin['MAF'], af_fin['A1'], ref_alleles['alt'])]
    ref_alleles['af_alt_tsi'] = [af if a1 == alt else 1-af for af, a1, alt in zip(af_tsi['MAF'], af_tsi['A1'], ref_alleles['alt'])]
    ref_alleles['af_alt_yri'] = [af if a1 == alt else 1-af for af, a1, alt in zip(af_yri['MAF'], af_yri['A1'], ref_alleles['alt'])]
 
    effects = eqtls.join(ref_alleles, on='rsid', how='inner')

    snp_sd = np.sqrt(2*effects['af_alt_all']*(1-effects['af_alt_all']))
    effects['gbr'] = (2*effects['corr']/snp_sd)*effects['af_alt_gbr']
    effects['fin'] = (2*effects['corr']/snp_sd)*effects['af_alt_fin']
    effects['tsi'] = (2*effects['corr']/snp_sd)*effects['af_alt_tsi']
    effects['yri'] = (2*effects['corr']/snp_sd)*effects['af_alt_yri']

    gbr_ids = sample_info.loc[sample_info['pop']=='GBR'].index
    fin_ids = sample_info.loc[sample_info['pop']=='FIN'].index
    tsi_ids = sample_info.loc[sample_info['pop']=='TSI'].index
    yri_ids = sample_info.loc[sample_info['pop']=='YRI'].index

    exp_mat.columns = [gene.split('.')[0] for gene in exp_mat.columns]

    effects = effects.loc[effects.index.intersection(exp_mat.columns)]
    for gene, eff_gbr, eff_fin, eff_tsi, eff_yri in zip(effects.index, effects.gbr, effects.fin, effects.tsi, effects.yri):
      exp_mat.loc[gbr_ids, gene] = exp_mat.loc[gbr_ids, gene] - eff_gbr
      exp_mat.loc[fin_ids, gene] = exp_mat.loc[fin_ids, gene] - eff_fin
      exp_mat.loc[tsi_ids, gene] = exp_mat.loc[tsi_ids, gene] - eff_tsi
      exp_mat.loc[yri_ids, gene] = exp_mat.loc[yri_ids, gene] - eff_yri
    exp_mat.loc[:] = preprocessing.scale(exp_mat.values)
    exp_coef, geno_coef, rho = util.lr_cca(exp_mat.values, geno_pcs.values,
        dim_x=config['dim_exp'], dim_z=config['dim_z'])
    print(rho)
    sample_info['Coordinate 1'] = exp_coef[:, 0]
    sample_info['Coordinate 2'] = exp_coef[:, 1]
    plot_pop = sns.lmplot(
        'Coordinate 1', 'Coordinate 2', hue='pop', data=sample_info, fit_reg=False)
    plot_pop.savefig(output.plot)
