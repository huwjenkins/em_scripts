#!/usr/bin/env python
# Defocus plotter. Author: Huw Jenkins 19.06.19
# Updated to plot refined defocus from particles star file 19.02.20
from __future__ import print_function
import sys
import argparse 
import numpy as np
import matplotlib.pyplot as plt

def make_plots(star_file, output_file, cutoff, bins, density):
  defocusU_results = []
  defocusV_results = []
  astigmatism_results = []
  classes = {}
  u = None
  v = None
  n = None
  if cutoff is None:
    cutoff = 999999.99
  relion31 = False
  data_particles= False
  with open(star_file,'r') as f:
    for line in f:
      if line[0] == '#' and line.split()[1] == 'version' and line.split()[2][0] == '3':
        relion31 = True
        data_micrographs = False
        data_particles= False
      if relion31 and not (data_micrographs or data_particles):
        if 'data_micrographs' in line:
          data_micrographs = True
        elif 'data_particles' in line:
          data_particles = True
      else:
        items = line.split()
        if len(items) == 2 and items[0] == '_rlnClassNumber':
            n = int(items[1].replace('#',''))
        elif len(items) == 2 and items[0] == '_rlnDefocusU':
          u = int(items[1].replace('#',''))
        elif len(items) == 2 and items[0] == '_rlnDefocusV':
          v = int(items[1].replace('#',''))
        elif len(items) > 2 and items[0] != '#':
          if u is not None and v is not None:
            a = abs(float(items[u - 1]) - float(items[v - 1]))
            if a < cutoff and n is None:
              defocusU_results.append(float(items[u - 1]))
              defocusV_results.append(float(items[v - 1]))
              astigmatism_results.append(a)
            elif n is not None:
              try:
                classes[int(items[n - 1])]['defocusU_results'].append(float(items[u - 1]))
                classes[int(items[n - 1])]['defocusV_results'].append(float(items[v - 1]))
              except KeyError:
                classes[int(items[n - 1])] = {'defocusU_results':[float(items[u - 1])]}
                classes[int(items[n - 1])].update({'defocusV_results':[float(items[v - 1])]})
          else:
            sys.exit('Sorry could not find _rlnDefocusU or _rlnDefocusV in {}'.format(star_file))
  assert len(defocusU_results) == len(defocusV_results)
  if data_particles or any(n in star_file for n in ['data', 'particles', 'shiny']): 
    colors = ['#e69f00','#0072b2','#009e73','#cc79a7','#f0e442','#56b4e9','#d55e00','#999999']
    if n is None:
      classes = {'1':{'defocusU_results':defocusU_results, 'defocusV_results':defocusV_results}}
    if len(classes) == 1:
      colors = ['#0072b2']
    for i, cls in enumerate(sorted(classes, key=lambda c: len(classes[c]['defocusU_results']), reverse=True)):
      kwargs = dict(histtype='stepfilled', edgecolor='none', alpha=0.75, bins=bins, label='class {}'.format(cls), normed=density)
      u = np.array(classes[cls]['defocusU_results'])
      v = np.array(classes[cls]['defocusV_results'])
      d = (u+v)/2.0
      if i < 9:
        plt.hist(d, color=colors[i], **kwargs)
      else: 
        plt.hist(d, **kwargs)
    plt.xlabel('Defocus ($\AA$)')
    plt.ylabel('Number of particles')
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
  parser.add_argument('--density', required=False, default=False, action='store_true',
                      help='normalise (i.e. form a probability density')
  args = parser.parse_args()
  make_plots(star_file=args.star_file, output_file=args.output, bins=args.bins, cutoff=args.cutoff, density=args.density)
