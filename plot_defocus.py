#!/usr/bin/env python
# Defocus plotter. Author: Huw Jenkins 19.06.19
# Updated to plot refined defocus from particles star file 19.02.20
# Added class option 03.11.20
# Fixed width bins 08.12.20
# Better header reading 02.04.21
from __future__ import print_function
import sys
import argparse 
import numpy as np
import matplotlib.pyplot as plt

def read_headers(star_file):
  got_labels = False
  data = False
  in_loop = False
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

def make_plots(star_file, output_file, cutoff, bins, select):
  defocusU_results = []
  defocusV_results = []
  astigmatism_results = []
  classes = {}
  if cutoff is None:
    cutoff = 999999.99
  labels = read_headers(star_file)
  u = v = n = None
  try:
    n, u, v = labels.index('rlnClassNumber'), labels.index('rlnDefocusU'), labels.index('rlnDefocusV')
    data_particles = True
  except ValueError:
    u, v = labels.index('rlnDefocusU'), labels.index('rlnDefocusV')
    data_particles = False

  with open(star_file) as f:
    data = False
    for line in f:
      if data and line.strip() != '':
        items = line.split()
        if u is not None and v is not None:
          a = abs(float(items[u]) - float(items[v]))
          if a < cutoff and n is None:
            defocusU_results.append(float(items[u]))
            defocusV_results.append(float(items[v]))
            astigmatism_results.append(a)
          elif n is not None:
            try:
              classes[int(items[n])]['defocusU_results'].append(float(items[u]))
              classes[int(items[n])]['defocusV_results'].append(float(items[v]))
            except KeyError:
              classes[int(items[n])] = {'defocusU_results':[float(items[u])]}
              classes[int(items[n])].update({'defocusV_results':[float(items[v])]})
          else:
            sys.exit('Sorry could not find _rlnDefocusU or _rlnDefocusV in {}'.format(star_file))
      elif line.startswith('_' + labels[-1]):
        data = True
      else:
        continue

  assert len(defocusU_results) == len(defocusV_results)
  if data_particles or any(n in star_file for n in ['data', 'particles', 'shiny']): 
    colors = ['#e69f00','#0072b2','#009e73','#cc79a7','#f0e442','#56b4e9','#d55e00','#999999']
    if n is None:
      classes = {'1':{'defocusU_results':defocusU_results, 'defocusV_results':defocusV_results}}
    if select is not None:
      classes = {select:classes[select]}
      if output_file == 'defocus.pdf':
        output_file = 'defocus_class{}.pdf'.format(select)
    if len(classes) == 1:
      colors = ['#0072b2']
    dmin = 99999.0
    dmax = 0.0 
    data = []
    labels = []
    for cls in sorted(classes, key=lambda c: len(classes[c]['defocusU_results']), reverse=True):
      u = np.array(classes[cls]['defocusU_results'])
      v = np.array(classes[cls]['defocusV_results'])
      d = (u+v)/2.0
      if len(classes) == 1:
        pct = [0.1,5,10,25,50,75,90,95,99.9]
        pc = np.percentile(d,pct)
        for j, p in enumerate(pct):
          print('{:4.1f}% particles have defocus < {:.0f} A'.format(p, pc[j]))
      if np.min(d) < dmin: dmin = np.min(d)
      if np.max(d) > dmax: dmax = np.max(d)
      data.append(d)
      labels.append('class {}'.format(cls))
    kwargs = dict(histtype='stepfilled', edgecolor='none', alpha=0.75, bins=bins, range=(dmin, dmax))
    for i, d in enumerate(data):
      if i < 9:
        plt.hist(d, color=colors[i], label=labels[i], **kwargs)
      else: 
        plt.hist(d, **kwargs)
    plt.xlabel('Defocus ($\AA$)')
    plt.ylabel('Number of particles')
    if len(classes) > 1:
      plt.legend(loc='best', fontsize=10)
  else:
    u = np.array(defocusU_results)
    v = np.array(defocusV_results)
    a = np.array(astigmatism_results)
    plt.subplot2grid((2,2), (0,0), colspan=2)
    plt.scatter(u,v)
    plt.title(star_file)
    plt.xlabel('Defocus U ($\AA$)')
    plt.ylabel('Defocus V ($\AA$)')
    plt.subplot2grid((2,2), (1,0), colspan=2)
    plt.hist(a, bins=bins)
    plt.xlabel('Astigmatism ($\AA$)')
    plt.ylabel('Number of micrographs')
  if cutoff != 999999.99:
    print('Writing defocus results with astigmatism lower than than {:0.2f} to {}'.format(cutoff, output_file))
  else:
    print('Writing defocus results to {}'.format(output_file))
  plt.savefig(output_file, format='pdf')
  plt.close()

if __name__=='__main__':
  parser = argparse.ArgumentParser(description='Plot per micrograph or per particle defocus')
  parser.add_argument('star_file', metavar='[micrographs_ctf.star, particles_ctf_refine.star]', type=str,
                      help='star file from CtfFind or with refined CTF parameters')
  parser.add_argument('--output', required=False, default='defocus.pdf', metavar='defocus.pdf', type=str,
                      help='output file_name')
  parser.add_argument('--cutoff', required=False, default=None, metavar='999.9', type=float,
                      help='maximum astigmatism to include in plots (helpful to reject very poor micrographs)')
  parser.add_argument('--bins', required=False, default=60, metavar='60', type=int,
                      help='number of bins in histogram')
  parser.add_argument('--select_class', required=False, default=None, metavar='1', type=int,
                      help='just plot this class')
  args = parser.parse_args()
  make_plots(star_file=args.star_file, output_file=args.output, bins=args.bins, cutoff=args.cutoff, select=args.select_class)
