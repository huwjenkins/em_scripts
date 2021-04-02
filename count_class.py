#! /usr/bin/env python
# Count particles in classes in RELION star file. Author Huw Jenkins 2019
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

def count_particles(star_files, sort_reso):
  for star_file in sorted(star_files, key=lambda x: int(x[x.find('_it') + 3:x.find('_data.star')])): 
    iteration = int(star_file[star_file.find('_it') + 3:star_file.find('_data.star')])
    classes = {}
    labels = read_headers(star_file)
    n = labels.index('rlnClassNumber')
    with open(star_file) as f:
      data = False
      for line in f:
        if data and line.strip() != '':
          items = line.split()
          try:
            classes[int(items[n])] +=1
          except KeyError:
            classes[int(items[n])] = 1
        elif line.startswith('_' + labels[-1]):
          data = True
        else:
          continue

    total = 0
    with open(star_file.replace('data','model'), 'r') as f:
      cls = 0
      data = False
      for line in f:
        if 'data_model_classes' in line:
          data = True
          labels = []
        elif data and line[0] == '_':
          labels.append(line[line.find('_') + 1:line.find('#') - 1])
        elif data and '_'.join(os.path.split(star_file)[1].split('_')[:2]) in line:
          cls += 1
          try:
            classes[cls] = (classes[cls], float(line.split()[labels.index('rlnEstimatedResolution')].replace('inf', '9999.99')))
          except KeyError:
            classes[cls] = (0, 9999.99) 
        elif 'data_model_class_1' in line:
          break

    print('Itn {:3d} Class #ptcls  Resn'.format(iteration))
    unclassified = 0
    for cls in sorted(classes):
      try:
        print('        {:3d}   {:7d} {:8.5f}'.format(cls, classes[cls][0], classes[cls][1]))
      except TypeError: #Fast subsets
        unclassified += classes[cls]
    if unclassified > 0: # fast subsets
      print('          -   {:7d}'.format(unclassified))
  if len(star_files) == 1 and unclassified == 0: 
    # only print sorted list for single data.star file
    print('\nItn {:3d} Class #ptcls  Resn'.format(iteration))
    for cls in (sorted(classes.items(), key=lambda x: x[1][sort_reso], reverse=not(sort_reso))):
      print('        {:3d}   {:7d} {:8.5f}'.format(cls[0], classes[cls[0]][0], classes[cls[0]][1]))
if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Count particles in each class in run_itNNN_data.star files')
  parser.add_argument('star_files', metavar='run_it0*_data.star', type=str, nargs='+',
                      help='list of star files (use * or ?? to match multiple files')
  parser.add_argument('--reso', required=False, default=False, action='store_true',
                      help='sort classes by resolution (instead of by No. of particles)')
  args = parser.parse_args()
  try:
    args.star_files.remove('run_it000_data.star') # classes > nclass
  except ValueError:
    pass
  if len([f for f in args.star_files if '*' in f]) > 0:
    sys.exit('Error: incorrect wildcard specified')
  if len([f for f in args.star_files if 'data' in f]) != len(args.star_files):
    sys.exit('Error: You need to give a list of run_itNNN_data files')
  count_particles(star_files=args.star_files, sort_reso=args.reso)
