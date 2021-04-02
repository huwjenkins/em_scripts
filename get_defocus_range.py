#! /usr/bin/env python
# Get median defocus for particles per micrograph in RELION star file. Author Huw Jenkins 2020
# Better header reading 02.04.21
from __future__ import print_function
import os
import sys
import argparse
import numpy as np

def read_headers(star_file):
  got_labels = False
  data = False
  in_loop = False
  with open(star_file) as f:
    for line in f:
      if line[0:5] == 'data_' and line.strip()[5:] in ['', 'particles']:
        data = True
        labels = []
      elif data and line[0] == '_':
        labels.append(line[line.find('_') + 1:line.find('#') - 1])
      elif data and len(labels) > 0:
        return labels
      else:
         continue 

def print_defocus_range(star_file, cutoff):
  results = {}
  labels = read_headers(star_file)
  m, u, v = labels.index('rlnMicrographName'), labels.index('rlnDefocusU'), labels.index('rlnDefocusV')
  with open(star_file) as f:
    data = False
    for line in f:
      if data and line.strip() != '':
        items = line.split()
        try:
          results[items[m].split('/')[-1]]['defocusU_results'].append(float(items[u]))
          results[items[m].split('/')[-1]]['defocusV_results'].append(float(items[v]))
        except KeyError:
          results[items[m].split('/')[-1]] = {'defocusU_results':[float(items[u])]}
          results[items[m].split('/')[-1]].update({'defocusV_results':[float(items[v])]})
      elif line.startswith('_' + labels[-1]):
        data = True
      else:
        continue

  if cutoff is not None:
    print('Micrograph                                                         median   mean     max      num > cutoff')
  else:
    print('Micrograph                                                         median   mean     max      no. ptcls')
  summary = {}
  for mic in results:
    d = (np.array(results[mic]['defocusU_results']) + np.array(results[mic]['defocusV_results']))/2.0
    five_num = np.percentile(d, [0, 25, 50, 75, 100])
    mean = np.mean(d)
    
    if cutoff is not None:
      summary.update({mic:[five_num[2],mean,five_num[4],(d>cutoff).sum(), d.size]})
    else:
      summary.update({mic:[five_num[2],mean,five_num[4],d.size]})
  for mic in sorted(summary.items(), key=lambda x: x[1][0]):
    if cutoff is not None:
      print(mic[0], '{:8.1f} {:8.1f} {:8.1f} {:4d}/{:4d}'.format(mic[1][0],mic[1][1],mic[1][2],mic[1][3],mic[1][4]))
    else:
      print(mic[0], '{:8.1f} {:8.1f} {:8.1f} {:4d}'.format(mic[1][0],mic[1][1],mic[1][2],mic[1][3]))

if __name__=='__main__':
  parser = argparse.ArgumentParser(description='Print per micrograph defocus spread')
  parser.add_argument('star_file', metavar='[run_data.star, particles_ctf_refine.star]', type=str,
                      help='star file with refined CTF parameters')
  parser.add_argument('--cutoff', required=False, default=None, metavar='15000', type=int,
                      help='write number of particles with defocus below this cutoff on each micrograph')
  args = parser.parse_args()
  print_defocus_range(star_file=args.star_file, cutoff=args.cutoff)
