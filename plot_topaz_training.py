#! /usr/bin/env python
# Plots progress of Topaz training when run from RELION
# Author Huw Jenkins 235021
from __future__ import print_function
import os
import sys
import argparse
import pandas as pd
import matplotlib.pyplot as plt

def make_plot(star_file, output_file):
  if output_file == 'topaz_training.pdf': # default
    output_file = os.path.join(os.path.split(star_file)[0], output_file)
  with open(star_file) as f:
    for line in f:
      if line.strip() != '' and line[0] != '#':
        if line.startswith('_rlnJobTypeLabel'):
          job = line.split()[1].strip()
        elif line.split()[0].strip() == 'topaz_nr_particles':
          n = int(line.split()[1].strip())
    if n < 0 or job != 'relion.autopick.topaz.train':
      sys.exit("Sorry {} doesn't appear to be from a Topaz training job".format(star_file))
  results_file = os.path.join(os.path.split(star_file)[0],'model_training.txt')
  table = pd.read_csv(results_file, sep='\t')
  table = table.loc[table['split'] == 'test'] # only keep the validation results
  table['auprc'] = table['auprc'].astype(float)
  print('Plotting area under the precision-recall curve for {} epochs of training...'.format(len(table)))
  fig, ax = plt.subplots()
  ax.plot(table['epoch'], table['auprc'], '-o',  label=str(n))
  ax.set_xlabel('Epoch')
  ax.set_ylabel('AUPRC')
  ax.legend(loc='best')
  plt.savefig(output_file, format='pdf')
  print('...written plot to {}'.format(output_file))
  plt.close()

if __name__=='__main__':
  parser = argparse.ArgumentParser(description='Plot progress of Topaz training run through RELION')
  parser.add_argument('star_file', metavar='Autopick/jobNNN/job.star', type=str,
                      help='path from RELION job directory to star file from Topaz training')
  parser.add_argument('--output', required=False, default='topaz_training.pdf', metavar='topaz_training.pdf', type=str,
                      help='output file_name')
  args = parser.parse_args()
  if not os.path.split(args.star_file)[0].startswith('AutoPick'):
    sys.exit('Please run this script from the RELION job directory and supply the path to the star file as Autopick/jobNNN/job.star')
  make_plot(star_file=args.star_file, output_file=args.output)
