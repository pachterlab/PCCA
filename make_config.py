"""Script to make config.yaml and plink keep file.

Usage:
  python make_config.py --keep_ceu
"""

import argparse
import os
import pandas as pd

KEEP_IMP=False # Not a paramter since we're using phase1 genotypes

def main():
  """Creates config.yaml based on arguments.
  """
  parser = argparse.ArgumentParser()
  parser.add_argument(
      '--working_dir',
      help='Directory to store the output of the script. Default is local dir.',
      default=None)
  parser.add_argument(
      '--data_dir',
      help=('Directory to store reads, genotypes, references and annotations.'
            ' Default is raw_data in the working_dir.'),
      default=None)
  parser.add_argument(
      '--transcript_fa',
      help=('Path to a gzipped transcript fasta file to use. Default is to'
            ' download the gencode protien coding transcript fasta in'
            ' data_dir.'),
      default=None)
  parser.add_argument(
      '--gtf_file',
      help=('Path toa gzipped gtf file to use. Default us to download the'
            ' the gencode v27 gtf in data_dir.'),
      default=None)
  parser.add_argument(
      '--keep_ceu',
      help='Keep CEU individuals for analysis. Default is drop them.',
      action='store_true',
      default=False)
  parser.add_argument(
      '--drop_yri',
      help='Drop YRI individuals for analysis. Default is keep them.',
      action='store_false',
      default=True)
  parser.add_argument(
      '--sample_file',
      help='Path of sample info file, default is local directory.',
      default='sample_info.txt')
  parser.add_argument(
      '--config_file',
      help='Path to store the config file, default is local directory.',
      default='config.yaml')
  parser.add_argument(
      '--geu_id_file',
      help='Path to store id file for plink, default is local directory.',
      default='ids_to_process.txt')
  parser.add_argument(
      '--min_tpm',
      help='Minimum mean TPM to include transcript. Default is 0.0001',
      default=0.0001)
  parser.add_argument(
      '--min_maf',
      help='Minimum MAF to include SNP. Default is 0.01',
      default=0.01)
  parser.add_argument(
      '--dim_exp',
      help='Number of expression hidden dimensions to model. Default is 12',
      default=12)
  parser.add_argument(
      '--dim_snp',
      help='Number of SNP hidden dimensions to model. Default is 8',
      default=8)
  parser.add_argument(
      '--dim_z',
      help='Number of parent latent dimensions to model. Default is 2',
      default=2)
  parser.add_argument(
      '--n_perm',
      help=('Number of permutations to perform to test gene significance.'
            ' Default is 1M'),
      default=int(1e6))
  parser.add_argument(
      '--plink2_bin',
      help=('Location of the plink2 binary or name of it in your PATH.'
            ' Default is plink2 in PATH.'),
      default='plink2')
  parser.add_argument(
      '--correction',
      help='Correction method to use for gene expression. CCA or regression.',
      default='regression')
  parser.add_argument(
      '--model',
      help='Analysis model to use. CCA or regression.',
      default='CCA')
  args = parser.parse_args()

  sample_info = pd.read_csv(args.sample_file, sep='\t')
  if not args.keep_ceu:
    drop = sample_info[sample_info['pop'] == 'CEU'].index
    sample_info = sample_info.drop(drop)
  if not args.drop_yri:
    drop = sample_info[sample_info['pop'] == 'YRI'].index
    sample_info = sample_info.drop(drop)
  if not KEEP_IMP:
    drop = sample_info[sample_info['imputed'] == 1].index
    sample_info = sample_info.drop(drop)
  sample_info['fid'] = 0
  #  Don't rewrite the geu_id_file unless it has changed so we dont trigger
  #  Snakemake to re-run some rules.
  try:
    id_file_old = pd.read_csv(
        args.geu_id_file, sep=' ', header=None, names=['fid', 'sample_name'])
    if set(id_file_old['sample_name']).symmetric_difference(sample_info['sample_name']):
        sample_info[['fid', 'sample_name']].to_csv(
          args.geu_id_file, sep='\t', index=False)
  except IOError:
    sample_info[['fid', 'sample_name']].to_csv(
      args.geu_id_file, sep='\t', index=False)    

  #  TODO(brielin): Rewrite this whole thing to use a yaml or json package.
  with open(args.config_file, 'w') as config_file:
    config_file.write('min_maf: ' + str(args.min_maf) + '\n')
    config_file.write('min_tpm: ' + str(args.min_tpm) + '\n')
    config_file.write('dim_exp: ' + str(args.dim_exp) + '\n')
    config_file.write('dim_snp: ' + str(args.dim_snp) + '\n')
    config_file.write('dim_z: ' + str(args.dim_z) + '\n')
    config_file.write('n_perm: ' + str(args.n_perm) + '\n')
    config_file.write('plink2_bin: ' + args.plink2_bin + '\n')
    config_file.write('correction: ' + args.correction + '\n')
    config_file.write('model: ' + args.model + '\n')
    if args.working_dir is not None:
      config_file.write('working_dir: ' + args.working_dir + '\n')
    else:
      config_file.write('working_dir: ' + os.getcwd() + '\n')
    if args.data_dir is not None:
      config_file.write('data_dir: ' + args.data_dir + '\n')
      data_dir = args.data_dir
    else:
      config_file.write('data_dir: ' + os.getcwd() + '/raw_data\n')
      data_dir = os.getcwd() + '/raw_data'
    if args.transcript_fa is not None:
      config_file.write('transcript_fa: ' + args.transcript_fa + '\n')
    else:
      config_file.write(''.join(
          ['transcript_fa: ', data_dir,
           '/gencode.v27.pc_transcripts.fa.gz\n']))
    if args.gtf_file is not None:
      config_file.write('gtf_file: ' + args.gtf_file + '\n')
    else:
      config_file.write(''.join(
          ['gtf_file: ', data_dir,
           '/gencode.v27.annotation.gtf.gz\n']))
    config_file.write('samples:\n')
    for sample, string in sample_info[['sample_name', 'sample_string']].values:
      config_file.write('    ' + sample + ': ' + string + '\n')

if __name__ == '__main__':
  main()
