#! /usr/bin/env python
# Plot class distributions and log-likelihood from RELION model.star files. Author: Huw Jenkins 19.05.20
from __future__ import print_function
import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt

def make_plot(star_files, output_file):
  n_itns = len(star_files)
  itn = []
  ll = []
  dist = []
  for star_file in sorted(star_files, key=lambda x: int(x[x.find('_it') + 3:x.find('_model.star')])): 
    itn.append(int(star_file[star_file.find('_it') + 3:star_file.find('_model.star')]))
    data = False
    with open(star_file) as f:
      data = False
      c = []
      for line in f:
        if line[0] == '_' and line[line.find('_') + 1:line.find(' ')] == 'rlnLogLikelihood':
          ll.append(float(line.split()[1]))
        elif line[0] == '_' and line[line.find('_') + 1:line.find(' ')] == 'rlnNrClasses':
          n_classes = int(line.split()[1])
        elif 'data_model_classes' in line:
          data = True
          labels = []
        elif data and line[0] == '_':
          labels.append(line[line.find('_') + 1:line.find('#') - 1])
        elif data and '_'.join(os.path.split(star_file)[1].split('_')[:2]) in line:
          c.append(float(line.split()[labels.index('rlnClassDistribution')]))
        elif 'data_model_class_1' in line:
          assert len(c) == n_classes
          dist.append(c)
          break

  i = np.array(itn)
  l = np.array(ll)
  d = np.array(dist).T
  colors = ['#e69f00','#0072b2','#009e73','#cc79a7','#f0e442','#56b4e9','#d55e00','#999999']
  fig, ax = plt.subplots()
  for c in range(n_classes):
    if c < 8:
      ax.plot(i, d[c], '-', linewidth=2, color=colors[c], label='Class {}'.format(c+1))
    else:
      ax.plot(i, d[c], '--', linewidth=2, color=colors[c-8], label='Class {}'.format(c+1))
  ax.set_xlabel('Iteration')
  ax.set_ylabel('Class Distribution')
  ax.legend(loc='upper left', fontsize=10)
  ax.tick_params(axis='x', direction='out')
  ax2 = ax.twinx()
  ax2.plot(itn, l, '-', linewidth=2, color='#000000', label='LogLikelihood')
  ax2.set_ylim(0)
  ax2.set_ylabel('LogLikelihood')
  ax2.legend(loc='upper center', fontsize=10)
  print('Writing results to {}'.format(output_file))
  fig.savefig(output_file, format='pdf')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Plot progress of classification from run_itNNN_model.star files')
  parser.add_argument('star_files', metavar='run_it0*_model.star', type=str, nargs='+',
                      help='list of star files (use * or ?? to match multiple files')
  parser.add_argument('--output', required=False, default='iterations.pdf', metavar='defocus.pdf', type=str,
                      help='output file_name')
  args = parser.parse_args()
  try:
    args.star_files.remove('run_it000_data.star') # classes > nclass
  except ValueError:
    pass
  if len([f for f in args.star_files if 'model' in f]) != len(args.star_files):
    sys.exit('Error: You need to give a list of run_itNNN_model.star files')
  make_plot(star_files=args.star_files, output_file=args.output)
