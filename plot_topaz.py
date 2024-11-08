#! /usr/bin/env python
# Plots progress of Topaz training or Auto-picking when run from RELION
# replaces plot_topaz_training.py and plot_topaz_results.py
# Author Huw Jenkins 061223
# 081124 add Table of No. particles at various FOM thresholds
# 081124 Allow plotting multiple training runs.

from __future__ import print_function
import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

def read_headers(star_file):
  data = False
  with open(star_file) as f:
    for line in f:
      if line[0:5] == 'data_' and line.strip()[5:] in ['', 'coordinate_files']:
        data = True
        labels = []
      elif data and line[0] == '_':
        labels.append(line[line.find('_') + 1:line.find('#') - 1])
      elif data and len(labels) > 0:
        return labels
      else:
         continue

def get_star_files(star_file):
  star_files = []
  labels=read_headers(star_file)
  sf = labels.index('rlnMicrographCoordinates')
  with open(star_file) as f:
    data = False
    for line in f:
      if data and line.strip() != '':
        items = line.split()
        star_files.append(items[sf])
      elif line.startswith('_' + labels[-1]):
        data = True
      else:
        continue
  return star_files

def get_job_type(star_file):
  with open(star_file) as f:
    for line in f:
      if line.strip() != '' and line[0] != '#':
        if line.startswith('_rlnJobTypeLabel'):
          job = line.split()[1].strip()
        elif line.split()[0].strip() == 'topaz_nr_particles':
          n = int(line.split()[1].strip())
  if job not in ['relion.autopick.topaz.train', 'relion.autopick.topaz.pick']:
    sys.exit("Sorry {} doesn't appear to be from a Topaz training or auto-picking job".format(star_file))
  return job, n

def make_FOM_plot(star_file, output_file, min, max, bins):
  foms = []
  star_files = get_star_files(star_file)
  print('Reading FOMs from {} star files...'.format(len(star_files)))
  for sf in star_files:
    labels=read_headers(sf)
    fom = labels.index('rlnAutopickFigureOfMerit')
    with open(sf) as f:
      data = False
      for line in f:
        if data and line.strip() != '':
          items = line.split()
          foms.append(float(items[fom]))
        elif line.startswith('_' + labels[-1]):
          data = True
        else:
          continue
  print('Plotting histogram of FOM from {} picks from {} micrographs...'.format(len(foms), len(star_files)))
  a = np.array(foms)
  print(' FOM  No. ptcls')
  for t in [0.0, -1.0, -1.5,  -2.0, -2.5, -3.0, -3.5, -4.0, -4.5, -5, -6]:
      print(f"{'{:4.1f}'.format(t):>3} {np.sum(a>t):7d}")
  fig, ax1 = plt.subplots()
  ax1.hist(a, bins=bins, range=(min, max))
  ax1.set_xlabel('Predicted score (predicted log-likelihood ratio)')
  ax1.set_ylabel('Number of particles')
  ax1.xaxis.set_minor_locator(AutoMinorLocator())
  plt.grid(True)
  plt.savefig(output_file, format='pdf')
  print('...written plot to {}'.format(output_file))
  plt.close()

def make_training_plot(star_files, n, output_file):
  fig, ax = plt.subplots()
  for star_file, n in zip(star_files,n):
    if n == -1:
      print('WARNING {} has number of expected particles set to -1.'.format(star_file))
      print('WARNING RELION sets the default number of expected particles to 200 but *you* should have set this!')
      print('WARNING Topaz training quality is highly dependent on the value for the number of expected particles')
      print('WARNING You should optimise the choice of value for the number of expected particles!')
    results_file = os.path.join(os.path.split(star_file)[0],'model_training.txt')
    table = pd.read_csv(results_file, sep='\t')
    table = table.loc[table['split'] == 'test'] # only keep the validation results
    table['auprc'] = table['auprc'].astype(float)
    print('Plotting area under the precision-recall curve for {} epochs of training...'.format(len(table)))
    ax.plot(table['epoch'], table['auprc'], '-o',  label=str(n))
  ax.set_xlabel('Epoch')
  ax.set_ylabel('AUPRC')
  ax.legend(loc='best')
  plt.savefig(output_file, format='pdf')
  print('...written plot to {}'.format(output_file))
  plt.close()

def make_plot(star_files, output_file, min, max, bins):
  if len(star_files) == 1:
    star_file=star_files[0]
    job, n = get_job_type(star_file)
    if output_file == 'topaz.pdf': # default
      if job == 'relion.autopick.topaz.pick':
        output_file = os.path.join(os.path.split(star_file)[0], 'topaz_FOM.pdf')
      elif job == 'relion.autopick.topaz.train':
        output_file = os.path.join(os.path.split(star_file)[0], 'topaz_training.pdf')
    if job == 'relion.autopick.topaz.pick':
      star_file = os.path.join(os.path.split(star_file)[0],'autopick.star')
      make_FOM_plot(star_file, output_file, min, max, bins)
    elif job == 'relion.autopick.topaz.train':
      make_training_plot(star_files, [n], output_file)
  else:
    jobs = [get_job_type(sf)[0] for sf in star_files]
    n = [get_job_type(sf)[1] for sf in star_files]
    if 'relion.autopick.topaz.pick' in jobs:
      sys.exit('Sorry only multiple training runs can be plotted together')
    else:
      output_file = 'topaz_training.pdf'
      make_training_plot(star_files, n, output_file)

if __name__=='__main__':
  parser = argparse.ArgumentParser(description='Plot results from Topaz run through RELION')
  parser.add_argument('star_files', metavar='Autopick/jobNNN/job.star', type=str, nargs='+',
                      help='path(s) from RELION job directory to star file(s) from Autopick')
  parser.add_argument('--output', required=False, default='topaz.pdf', metavar='topaz.pdf', type=str,
                      help='output file_name')
  parser.add_argument('--min', required=False, default=-6, metavar='-6', type=float,
                      help='minimum score for FOM plot')
  parser.add_argument('--max', required=False, default=5, metavar='5', type=float,
                      help='maximum score for FOM plot')
  parser.add_argument('--bins', required=False, default=50, metavar='50', type=int,
                      help='number of bins in FOM histogram')
  args = parser.parse_args()
  for star_file in args.star_files:
    if not os.path.split(star_file)[0].startswith('AutoPick'):
      sys.exit('Please run this script from the RELION job directory and supply the path to the job.star file as Autopick/jobNNN/job.star')
    if not os.path.split(star_file)[-1] == 'job.star':
      sys.exit('Please run this script from the RELION job directory and supply the path to the job.star file as Autopick/jobNNN/job.star')
    if not os.path.isfile(star_file):
      sys.exit('Could not find {}'.format(star_file))
  make_plot(star_files=args.star_files, output_file=args.output, min=args.min, max=args.max, bins=args.bin
