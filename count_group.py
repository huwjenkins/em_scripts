#! /usr/bin/env python
# Count particles in groups in RELION star file. Author Huw Jenkins 2019
# Better header reading 02.04.21
from __future__ import print_function
import os
import sys
import argparse

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

def count_group(star_file, output_file, cutoff):
  groups = {}
  mics = {}
  regrouped = False
  labels = read_headers(star_file)
  try:
    mic, n = labels.index('rlnMicrographName'), labels.index('rlnGroupNumber')
  except ValueError:
    mic, n = labels.index('rlnMicrographName'), labels.index('rlnGroupName')
    regrouped = True
  total = 0
  with open(star_file) as f:
    data = False
    for line in f:
      if data and line.strip() != '':
        items = line.split()
        try:
          if not regrouped:
            groups[int(items[n])] +=1
          else:
            groups[items[n]] +=1
          total += 1
        except KeyError:
          if not regrouped:
            groups[int(items[n])] = 1
            mics[int(items[n])] = items[mic].split('/')[-1]
          else:
            groups[items[n]] = 1
          total += 1
      elif line.startswith('_' + labels[-1]):
        data = True
      else:
        continue

  running_total = 0
  print('Group   #ptcls    total  Micrograph')
  reject = []
  for grp in (sorted(groups.items(), key=lambda x: x[1],reverse=True)):
    if not regrouped:
      print('{:<5d} {:8d} {:8d}  {}'.format(grp[0], groups[grp[0]], total - running_total, mics[grp[0]]))
      if cutoff is not None and groups[grp[0]] < cutoff:
          reject.append(mics[grp[0]])
    else:
      print('{} {:8d} {:8d}'.format(grp[0], groups[grp[0]], total - running_total))
    running_total += groups[grp[0]]
  if cutoff is not None and len(reject) > 0:
    print ('Writing micrographs with fewer than {} particles to {}'.format(cutoff, output_file))
    with open(output_file, 'w') as f:
      for mic in reject:
        f.write(mic+'\n')
if __name__=='__main__':
  parser = argparse.ArgumentParser(description='Count number of particles in each group (micrograph)')
  parser.add_argument('star_file', metavar='[run_data.star, shiny.star, particles_ctf_refine.star]', type=str,
                      help='star file')
  parser.add_argument('--cutoff', required=False, default=None, metavar='50', type=int,
                      help='write list of micrographs with fewer than this many particles')
  parser.add_argument('--output', required=False, default='reject.txt', metavar='reject.txt', type=str,
                      help='output file_name')
  args = parser.parse_args()
  count_group(star_file=args.star_file, output_file=args.output, cutoff=args.cutoff)
