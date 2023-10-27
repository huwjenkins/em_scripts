#! /usr/bin/env python
# Plots histogram of Topaz picking scores from RELION autopick job
# Author Huw Jenkins 205021
from __future__ import print_function
import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

def read_headers(star_file):
  got_labels = False
  data = False
  in_loop = False
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

def make_plot(star_file, output_file, min, max, bins):
  if output_file == 'topaz_FOM.pdf': # default
    output_file = os.path.join(os.path.split(star_file)[0], output_file)
  foms = []
  star_files = get_star_files(star_file)
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
  fig, ax1 = plt.subplots()
  ax1.hist(a, bins=bins, range=(min, max))
  ax1.set_xlabel('Predicted score (predicted log-likelihood ratio)')
  ax1.set_ylabel('Number of particles')
  ax1.xaxis.set_minor_locator(AutoMinorLocator())
  plt.grid(True)
  plt.savefig(output_file, format='pdf')
  print('...written plot to {}'.format(output_file))
  plt.close()

if __name__=='__main__':
  parser = argparse.ArgumentParser(description='Plot predicted scores from Topaz run through RELION')
  parser.add_argument('star_file', metavar='Autopick/jobNNN/autopick.star', type=str,
                      help='path from RELION job directory to star file from Autopick')
  parser.add_argument('--output', required=False, default='topaz_FOM.pdf', metavar='topaz_FOM.pdf', type=str,
                      help='output file_name')
  parser.add_argument('--min', required=False, default=-6, metavar='-6', type=float,
                      help='minimum score for plot')
  parser.add_argument('--max', required=False, default=5, metavar='5', type=float,
                      help='maximum score for plot')
  parser.add_argument('--bins', required=False, default=50, metavar='50', type=int,
                      help='number of bins in histogram')
  args = parser.parse_args()
  if not os.path.split(args.star_file)[0].startswith('AutoPick'):
    sys.exit('Please run this script from the RELION job directory and supply the path to the AutoPick star file as Autopick/jobNNN/autopick.star')
  if not os.path.split(args.star_file)[-1] == 'autopick.star':
    sys.exit('Please run this script from the RELION job directory and supply the path to the AutoPick star file as Autopick/jobNNN/autopick.star')
  if not os.path.isfile(args.star_file):
    sys.exit('Could not find {}'.format(args.star_file))
  make_plot(star_file=args.star_file, output_file=args.output, min=args.min, max=args.max, bins=args.bins)
