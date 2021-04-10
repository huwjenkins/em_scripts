#! /usr/bin/env python
# Plot FSC curves from RELION postprocess.star files. Author Huw Jenkins 03.07.20
# Better header reading 02.04.21
from __future__ import print_function
import os
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
      if line[0:5] == 'data_' and line.strip()[5:] in ['', 'fsc']:
        data = True
        labels = []
      elif data and line[0] == '_':
        labels.append(line[line.find('_') + 1:line.find('#') - 1])
      elif data and len(labels) > 0:
        return labels
      else:
         continue 

def make_plot(star_files, output_file, show_legend, legend, colors):
  n_itns = len(star_files)
  curves = []
  invresols = []
  fscs = []
  for star_file in star_files: 
    labels = read_headers(star_file)
    r, f = labels.index('rlnResolution'), labels.index('rlnFourierShellCorrelationCorrected')
    with open(star_file) as sf:
      data = False
      invres = []
      fsc = []
      for line in sf: 
        if data and line.strip() != '' and line[0] != '#':
          if line[0:12] != 'data_guinier':
            items = line.split()
            invres.append(float(items[r]))
            fsc.append(float(items[f]))
          elif 'data_guinier' in line: 
            break
        elif line[0] == '_' and line[line.find('_') + 1:line.find(' ')] == 'rlnUnfilteredMapHalf1':
          curves.append('_'.join(line.split()[1].split('/')[0:2]))
        elif line.startswith('_' + labels[-1]):
          data = True
        else:
          continue

    invresols.append(invres)
    fscs.append(fsc)
    invresols = list({tuple(r) for r in invresols})
  i = np.vstack(invresols)
  f = np.array(fscs)
  if i.shape[0] != 1:
    sys.exit('Error: plotting FSCs with different resolution bins is not supported!')
  if i.shape[1] != f.shape[1]:
    sys.exit('Error: Mis-match between no. of resolution shells and FSC values!')
  if colors is None:
    colors = ['#0072b2','#e69f00','#009e73','#cc79a7','#f0e442','#56b4e9','#d55e00','#999999']
  if legend is not None:
    curves = legend
  fig, ax = plt.subplots()
  for c in range(len(curves)):
    if c < 9:
      ax.plot(i[0], f[c], '-', linewidth=2, color=colors[c], label=curves[c])
    else:
      ax.plot(i[0], f[c], '-', linewidth=2, label=curves[c])
  ax.axhline(y=0.143, ls='--', color='#000000')
  ax.set_xlabel('Resolution (1/$\AA$)')
  ax.set_ylabel('FSC')
  if show_legend:
    ax.legend(loc='center left', fontsize=10)
  ax.tick_params(axis='x', direction='out')
  print('Writing results to {}'.format(output_file))
  fig.savefig(output_file, format='pdf')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Plot FSC curves from postprocess.star files')
  parser.add_argument('star_files', metavar='jobNNN/postprocess.star', type=str, nargs='+',
                      help='list of star files')
  parser.add_argument('--output', required=False, default='FSC.pdf', metavar='FSC.pdf', type=str,
                      help='output file_name')
  parser.add_argument('--no_legend', required=False, action='store_true', default=False,
                      help='hide legend')
  parser.add_argument('--colors', required=False, default=None, metavar='"#969696,#0072b2,#e69f00"', type=str,
                      help='comma separated list of colours')
  parser.add_argument('--legend', required=False, default=None, metavar='"curve 1,curve 2,curve 3"', type=str,
                      help='comma separated list for curves in legend')
  args = parser.parse_args()
  if len([f for f in args.star_files if 'postprocess.star' in f]) != len(args.star_files):
    sys.exit('Error: You need to give a list of postprocess.star files')
  legend = args.legend.split(',') if args.legend is not None else None
  colors = args.colors.split(',') if args.colors is not None else None
  if legend is not None and len(legend) != len(args.star_files):
      sys.exit('Error: Mismatch between number of labels and number of star files')
  if colors is not None and len(colors) != len(args.star_files):
      sys.exit('Error: Mismatch between number of colours and number of star files')
  make_plot(star_files=args.star_files, output_file=args.output, show_legend=not(args.no_legend), legend=legend, colors=colors)
