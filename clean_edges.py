#!/usr/bin/env python
# Remove particles that will lie close to edge or outside micrograph after recentring. Author: Huw Jenkins 27.11.24

import argparse
import numpy as np
from scipy.spatial.transform import Rotation as R

def euler_angles2matrix_scipy(alpha, beta, gamma):
  # This reproduces result of Euler_angles2matrix() from RELION src/euler.cpp
  alpha = np.radians(alpha)
  beta  = np.radians(beta)
  gamma = np.radians(gamma)
  A = R.from_euler('ZYZ', (alpha, beta, gamma))
  return A.as_matrix().T

def print_info(particle_angpix, orig_angpix, recenter_x, recenter_y, recenter_z, mic_x, mic_y, distance):
  print(f"Micrographs have dimensions: {mic_x} x {mic_y} px and pixel size of {orig_angpix} A")
  if mic_x != mic_y:
    print(f"WARNING! Micrographs are not square - check orientation is correct!!")
  recenter = [recenter_x, recenter_y, recenter_z]
  recenter_ang = [c * particle_angpix for c in recenter]
  print(f"After applying recentring of {recenter} px ({recenter_ang} A assuming {particle_angpix} A/px in reference) ")
  print(f"Particles closer than {distance} px ({distance * orig_angpix} A) to edge of micrographs will be removed)")
  print(f"Remaining particles will have center in range {distance} - {mic_x - distance - 1} in X and {distance} - {mic_y - distance - 1} in Y")

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

def filter_particles(star_file, output_file, particle_angpix, orig_angpix, recenter_x, recenter_y, recenter_z, mic_x, mic_y, distance, verbose, debug):
  print(f"Reading particles from {star_file}....")
  rescale = particle_angpix / orig_angpix
  labels = read_headers(star_file)
  xo, yo = labels.index('rlnOriginXAngst'), labels.index('rlnOriginYAngst')
  x, y = labels.index('rlnCoordinateX'), labels.index('rlnCoordinateY')
  r, t, p = labels.index('rlnAngleRot'), labels.index('rlnAngleTilt'), labels.index('rlnAnglePsi')
  print_info(particle_angpix, orig_angpix, recenter_x, recenter_y, recenter_z, mic_x, mic_y, distance)
  n_rejected = 0
  n_retained = 0
  n_particles = 0
  with open(output_file, 'w') as fout:
    with open(star_file) as f:
      data = False
      for line in f:
        if data and line.strip() != '':
          n_particles+=1
          items = line.split()
          xoff, yoff = float(items[xo]), float(items [yo])
          rot, tilt, psi = float(items[r]), float(items[t]), float(items[p])
          xcoord, ycoord = float(items[x]), float(items[y])
          # This code reproduces result of getCoordinateMetaDataTable() from RELION src/preprocessing.cpp
          xoff /= particle_angpix # now in px
          yoff /= particle_angpix # now in px
          center = np.array([recenter_x, recenter_y, recenter_z])
          transform = euler_angles2matrix_scipy(rot, tilt, psi)
          projected_center = np.dot(transform, center)
          xoff -= projected_center[0]
          yoff -= projected_center[1]
          xoff *= rescale # now in micrograph px
          yoff *= rescale # now in micrograph px
          xcoord -= round(xoff, 0)
          ycoord -= round(yoff, 0)
          if debug:
            print(f"center: {center}")
            print("transform:")
            print(transform)
            print(f"projected_center: {projected_center}")
            print(f"coordinates: [{xcoord:4.0f}, {ycoord:4.0f}, 0]")
          if xcoord < distance or xcoord >= mic_x - distance or ycoord < distance or ycoord >= mic_y - distance:
            if verbose:
              print(f"Particle with centre: {xcoord:4.0f} {ycoord:4.0f} removed")
            n_rejected+=1
          else:
            n_retained+=1
            fout.write(line)
        elif line.startswith('_' + labels[-1]):
          data = True
          fout.write(line)
        else:
          fout.write(line)
          continue
  print(f"{n_rejected} of {n_particles} particles removed.")
  print(f"...{n_retained} particles written to {output_file}")

if __name__=='__main__':
  parser = argparse.ArgumentParser(description='Remove particles that will be close to edge (or off edge) of micrograph after recentring')
  parser.add_argument('star_file', metavar='run_data.star', type=str,
                      help='star file from Refine3D or Class3D')
  parser.add_argument('--output_file', required=False, default='filtered.star', metavar='filtered.star', type=str,
                      help='star file from Refine3D or Class3D')
  parser.add_argument('--particle_angpix', required=True, default=None, metavar='1.0', type=float,
                      help='Reference A/pix')
  parser.add_argument('--orig_angpix', required=True, default=None, metavar='1.0', type=float,
                      help='A/pix in micrographs')
  parser.add_argument('--mic_x', required=True, default=None, metavar='4096', type=int,
                      help='Micrograph size in X in pixels (Falcon 4: 4096, K2: 3838, K3: 5760)')
  parser.add_argument('--mic_y', required=True, default=None, metavar='4096', type=int,
                      help='Micrograph size in Y in pixels (Falcon 4: 4096, K2: 3710, K3: 4092)')
  parser.add_argument('--distance', required=False, default=0, metavar='0', type=int,
                      help='Remove particles with center closer than this distance to any edge of the micrograph')
  parser.add_argument('--recenter_x', required=False, default=0, metavar='0', type=int,
                      help='X-coordinate (in px of the reference) to recenter on')
  parser.add_argument('--recenter_y', required=False, default=0, metavar='0', type=int,
                      help='Y-coordinate (in px of the reference) to recenter on')
  parser.add_argument('--recenter_z', required=False, default=0, metavar='0', type=int,
                      help='Z-coordinate (in px of the reference) to recenter on')
  parser.add_argument('--verbose', required=False, default=False, action='store_true',
                      help='list particles that are removed')
  parser.add_argument('--debug', required=False, default=False, action='store_true',
                      help='print debugging information')
  args = parser.parse_args()
  filter_particles(star_file=args.star_file, 
                   output_file=args.output_file,
                   particle_angpix=args.particle_angpix,
                   orig_angpix=args.orig_angpix,
                   mic_x=args.mic_x,
                   mic_y=args.mic_y,
                   distance=args.distance,
                   recenter_x=args.recenter_x,
                   recenter_y=args.recenter_y,
                   recenter_z=args.recenter_z,
                   verbose=args.verbose,
                   debug=args.debug
                  )