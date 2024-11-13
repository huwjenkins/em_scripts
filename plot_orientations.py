#!/usr/bin/env python
# Orientation plotter. Author: Huw Jenkins 12.11.24
from __future__ import print_function
import argparse 
import numpy as np
import matplotlib.pyplot as plt

def read_headers(star_file):
  data = False
  with open(star_file) as f:
    for line in f:
      if line[0:5] == 'data_' and line.strip()[5:] in ['', 'particles', 'micrographs']:
        data = True
        labels = []
      elif data and line[0] == '_':
        labels.append(line[line.find('_') + 1:line.find('#') - 1])
      elif data and len(labels) > 0:
        return labels
      else:
         continue 

def make_plots(star_file, output_file, bins):
  labels = read_headers(star_file)
  r, t, p = labels.index('rlnAngleRot'), labels.index('rlnAngleTilt'), labels.index('rlnAnglePsi')
  with open(star_file) as f:
    rot = []
    tilt = []
    psi = []
    data = False
    for line in f:
      if data and line.strip() != '':
        items = line.split()
        rot.append(float(items[r]))
        tilt.append(float(items[t]))
        psi.append(float(items[p]))
      elif line.startswith('_' + labels[-1]):
        data = True
      else:
        continue

  r = np.array(rot)
  t = np.array(tilt)
  p = np.array(psi)
  kwargs = dict(color='#0072b2', histtype='stepfilled', edgecolor='none', alpha=0.75, bins=bins)
  plt.subplot2grid((15,1), (0,0), rowspan=3)
  plt.hist(r, range=(-180,180), **kwargs)
  plt.xticks(ticks=[-180, -135, -90, -45, 0, 45, 90, 135, 180])
  plt.xlabel('rlnAngleRot')
  plt.ylabel('No. particles')
  plt.subplot2grid((15,1), (5,0), rowspan=3)
  plt.hist(t, range=(0,180), **kwargs)
  plt.xticks(ticks=[0, 45, 90, 135, 180])
  plt.xlabel('rlnAngleTilt')
  plt.ylabel('No. particles')
  plt.subplot2grid((15,1), (10,0), rowspan=3)
  plt.hist(p, range=(-180,180), **kwargs)
  plt.xticks(ticks=[-180, -135, -90, -45, 0, 45, 90, 135, 180])
  plt.xlabel('rlnAnglePsi')
  plt.ylabel('No. particles')
  d = {'Rot':r, 'Tilt':t, 'Psi':p}
  for ang in d:
    print(f'Range of rlnAngle{ang}: {np.min(d[ang]):0.2f} - {np.max(d[ang]):0.2f}')
  print('Writing orientation results to {}'.format(output_file))
  plt.savefig(output_file, format='pdf')
  plt.close()

if __name__=='__main__':
  parser = argparse.ArgumentParser(description='Plot range of orientations in star file')
  parser.add_argument('star_file', metavar='[run_data.star]', type=str,
                      help='star file from Refine3D or Class3D')
  parser.add_argument('--output', required=False, default='orientations.pdf', metavar='orientations.pdf', type=str,
                      help='output file_name')
  parser.add_argument('--bins', required=False, default=180, metavar='180', type=int,
                      help='number of bins in histogram')
  args = parser.parse_args()
  make_plots(star_file=args.star_file, output_file=args.output, bins=args.bins)
