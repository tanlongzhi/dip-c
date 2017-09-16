# Python functions for NucDynamics

import numpy as np
import dyn_util

def load_ncc_file(file_path):
  """Load chromosome and contact data from NCC format file, as output from NucProcess"""
  
  if file_path.endswith('.gz'):
    import gzip
    file_obj = gzip.open(file_path)
  
  else:
    file_obj = open(file_path) 
  
  # Observations are treated individually in single-cell Hi-C,
  # i.e. no binning, so num_obs always 1 for each contact
  num_obs = 1  
    
  contact_dict = {}
  chromosomes = set()
    
  for line in file_obj:
    chr_a, f_start_a, f_end_a, start_a, end_a, strand_a, chr_b, f_start_b, f_end_b, start_b, end_b, strand_b, ambig_group, pair_id, swap_pair = line.split()
    
    if strand_a == '+':
      pos_a = int(f_start_a)
    else:
      pos_a = int(f_end_a)
    
    if strand_b == '+':
      pos_b = int(f_start_b)       
    else:
      pos_b = int(f_end_b)
 
    if chr_a > chr_b:
      chr_a, chr_b = chr_b, chr_a
      pos_a, pos_b = pos_b, pos_a
    
    if chr_a not in contact_dict:
      contact_dict[chr_a] = {}
      chromosomes.add(chr_a)
      
    if chr_b not in contact_dict[chr_a]:
      contact_dict[chr_a][chr_b] = [] 
      chromosomes.add(chr_b)
        
    contact_dict[chr_a][chr_b].append((pos_a, pos_b, num_obs, int(ambig_group)))
   
  file_obj.close()
  
  chromo_limits = {}
    
  for chr_a in contact_dict:
    for chr_b in contact_dict[chr_a]:
      contacts = np.array(contact_dict[chr_a][chr_b]).T
      contact_dict[chr_a][chr_b] = contacts
      
      seq_pos_a = contacts[1]
      seq_pos_b = contacts[2]
      
      min_a = min(seq_pos_a)
      max_a = max(seq_pos_a)
      min_b = min(seq_pos_b)
      max_b = max(seq_pos_b)
        
      if chr_a in chromo_limits:
        prev_min, prev_max = chromo_limits[chr_a]
        chromo_limits[chr_a] = [min(prev_min, min_a), max(prev_max, max_a)]
      else:
        chromo_limits[chr_a] = [min_a, max_a]
      
      if chr_b in chromo_limits:
        prev_min, prev_max = chromo_limits[chr_b]
        chromo_limits[chr_b] = [min(prev_min, min_b), max(prev_max, max_b)]
      else:
        chromo_limits[chr_b] = [min_b, max_b]
         
  chromosomes = sorted(chromosomes)      
        
  return chromosomes, chromo_limits, contact_dict
  
digits_upper = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
digits_lower = digits_upper.lower()
digits_upper_values = dict([pair for pair in zip(digits_upper, range(36))])
digits_lower_values = dict([pair for pair in zip(digits_lower, range(36))])

def encode_pure(digits, value):
  "encodes value using the given digits"
  assert value >= 0
  if (value == 0): return digits[0]
  n = len(digits)
  result = []
  while (value != 0):
    rest = value // n
    result.append(digits[value - rest * n])
    value = rest
  result.reverse()
  return "".join(result)

def decode_pure(digits_values, s):
  "decodes the string s using the digit, value associations for each character"
  result = 0
  n = len(digits_values)
  for c in s:
    result *= n
    result += digits_values[c]
  return result

def hy36encode(width, value):
  "encodes value as base-10/upper-case base-36/lower-case base-36 hybrid"
  i = value
  if (i >= 1-10**(width-1)):
    if (i < 10**width):
      return ("%%%dd" % width) % i
    i -= 10**width
    if (i < 26*36**(width-1)):
      i += 10*36**(width-1)
      return encode_pure(digits_upper, i)
    i -= 26*36**(width-1)
    if (i < 26*36**(width-1)):
      i += 10*36**(width-1)
      return encode_pure(digits_lower, i)
  raise ValueError("value out of range.")

def hy36decode(width, s):
  "decodes base-10/upper-case base-36/lower-case base-36 hybrid"
  if (len(s) == width):
    f = s[0]
    if (f == "-" or f == " " or f.isdigit()):
      try: return int(s)
      except ValueError: pass
      if (s == " "*width): return 0
    elif (f in digits_upper_values):
      try: return decode_pure(
        digits_values=digits_upper_values, s=s) - 10*36**(width-1) + 10**width
      except KeyError: pass
    elif (f in digits_lower_values):
      try: return decode_pure(
        digits_values=digits_lower_values, s=s) + 16*36**(width-1) + 10**width
      except KeyError: pass
  raise ValueError("invalid number literal.")

def export_pdb_coords(file_path, coords_dict, seq_pos_dict, particle_size, scale=1.0, extended=True):
  """
  Write chromosome particle coordinates as a PDB format file
  """

  alc = ' '
  ins = ' '
  prefix = 'HETATM'
  line_format = '%-80.80s\n'
  
  if extended:
    pdb_format = '%-6.6s%s %4.4s%s%3.3s %s%4.1d%s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2.2s  %10d\n'
    ter_format = '%-6.6s%s      %s %s%4.1d%s                                                     %10d\n'
  else:
    pdb_format = '%-6.6s%5.1d %4.4s%s%3.3s %s%4.1d%s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2.2s  \n'
    ter_format = '%-6.6s%5.1d      %s %s%4.1d%s                                                     \n'

  file_obj = open(file_path, 'w')
  write = file_obj.write
  
  chromosomes = list(seq_pos_dict.keys())
  num_models = len(coords_dict[chromosomes[0]])
  title = 'NucDynamics genome structure export'
      
  pos_chromo = {}
  for m in range(num_models):
    
    c = 0
    j = 1
    seqPrev = None
    for chromo in seq_pos_dict:
      
      if chromo.isdigit():
        idx = int(chromo) - 1
        chain_code = chr(ord('A')+idx)
      
      elif len(chromo) == 1:
        chain_code = chromo.upper()
        
      else:
        idx = chromosomes.index(chromo)
        chain_code = chr(ord('a')+idx)
      
      tlc = chromo
      while len(tlc) < 2:
        tlc = '_' + tlc
      
      if len(tlc) == 2:
        tlc = 'C' + tlc
      
      if len(tlc) > 3:
        tlc = tlc[:3]
      
      chromo_model_coords = coords_dict[chromo][m]
      
      if not len(chromo_model_coords):
        continue
      
      pos = seq_pos_dict[chromo]
      
      for i, seqPos in enumerate(pos):
        c += 1
 
        seqMb = int(seqPos//1e6) + 1
        
        if seqMb == seqPrev:
          j += 1
        else:
          j = 1
        
        el = 'C'
        a = 'C%d' % j
          
        aName = '%-3s' % a
        x, y, z = chromo_model_coords[i] #XYZ coordinates
         
        seqPrev = seqMb
        pos_chromo[c] = chromo
                 
        line = str(chromo)+'\t'+str(seqPos)+'\t'+str(x)+'\t'+str(y)+'\t'+str(z)+'\n'
        write(line)

  file_obj.close()


def remove_isolated_contacts(contact_dict, threshold=int(2e6)):
  """
  Select only contacts which are within a given sequence separation of another
  contact, for the same chromosome pair
  """
    
  for chromoA in contact_dict:
    for chromoB in contact_dict[chromoA]:
      contacts = contact_dict[chromoA][chromoB]
      positions = np.array(contacts[:2], np.int32).T
      
      if len(positions): # Sometimes empty e.g. for MT, Y chromos 
        active_idx = dyn_util.getSupportedPairs(positions, np.int32(threshold))
        contact_dict[chromoA][chromoB] = contacts[:,active_idx]
     
  return contact_dict
  

def remove_violated_contacts(contact_dict, coords_dict, particle_seq_pos, particle_size, threshold=5.0):  
  """
  Remove contacts whith structure distances that exceed a given threshold
  """
  
  for chr_a in contact_dict:
    for chr_b in contact_dict[chr_a]:
      contacts = contact_dict[chr_a][chr_b]
      
      contact_pos_a = contacts[0].astype(np.int32)
      contact_pos_b = contacts[1].astype(np.int32)
      
      coords_a = coords_dict[chr_a]
      coords_b = coords_dict[chr_b]

      struc_dists = []
 
      for m in range(len(coords_a)):
        coord_data_a = dyn_util.getInterpolatedCoords([chr_a], {chr_a:contact_pos_a}, particle_seq_pos, coords_a[m])
        coord_data_b = dyn_util.getInterpolatedCoords([chr_b], {chr_b:contact_pos_b}, particle_seq_pos, coords_b[m])
 
        deltas = coord_data_a - coord_data_b
        dists = np.sqrt((deltas*deltas).sum(axis=1))
        struc_dists.append(dists)
      
      # Average over all conformational models
      struc_dists = np.array(struc_dists).T.mean(axis=1)
      
      # Select contacts with distances below distance threshold
      indices = (struc_dists < threshold).nonzero()[0]
      contact_dict[chr_a][chr_b] = contacts[:,indices]
        
  return contact_dict
  
  
def get_random_coords(pos_dict, chromosomes, num_models, radius=10.0):
  """
  Get random, uniformly sampled coorinate positions, restricted to
  a sphere of given radius
  """
    
  from numpy.random import uniform
    
  num_particles = sum([len(pos_dict[chromo]) for chromo in chromosomes])
  coords = np.empty((num_models, num_particles, 3))
  r2 = radius*radius
    
  for m in range(num_models):
    
    for i in range(num_particles):
      x = y = z = radius

      while x*x + y*y + z*z >= r2:
        x = radius * (2*uniform(0,1) - 1)
        y = radius * (2*uniform(0,1) - 1)
        z = radius * (2*uniform(0,1) - 1)
      
      coords[m,i] = [x,y,z]

  return coords
  
  
def pack_chromo_coords(coords_dict, chromosomes):
  """
  Place chromosome 3D coordinates stored in a dictionary keyed by
  chromosome name into a single, ordered array. The chromosomes argument
  is required to set the correct array storage order.
  """

  chromo_num_particles = [len(coords_dict[chromo][0]) for chromo in chromosomes]
  n_particles = sum(chromo_num_particles)
  n_models = len(coords_dict[chromosomes[0]])  
  coords = np.empty((n_models, n_particles, 3), float)
  
  j = 0
  for i, chromo in enumerate(chromosomes):
    span = chromo_num_particles[i]
    coords[:,j:j+span] = coords_dict[chromo]
    j += span
      
  return coords
  
 
def unpack_chromo_coords(coords, chromosomes, seq_pos_dict):
  """
  Exctract coords for multiple chromosomes stored in a single array into
  a dictionary, keyed by chromosome name. The chromosomes argument is required
  to get the correct array storage order.
  """

  chromo_num_particles = [len(seq_pos_dict[chromo]) for chromo in chromosomes]
  n_seq_pos = sum(chromo_num_particles)
  n_models, n_particles, dims = coords.shape

  if n_seq_pos != n_particles:
    msg = 'Model coordinates must be an array of num models x %d' % (n_seq_pos,)
    raise(Exception(msg))  
  
  coords_dict = {}
        
  j = 0
  for i, chromo in enumerate(chromosomes):
    span = chromo_num_particles[i]
    coords_dict[chromo] = coords[:,j:j+span] # all models, slice
    j += span
 
  return coords_dict
  
  
def anneal_genome(chromosomes, contact_dict, num_models, particle_size,
                  general_calc_params, anneal_params,
                  prev_seq_pos_dict=None, start_coords=None):
    """
    Use chromosome contact data to generate distance restraints and then
    apply a simulated annealing protocul to generate/refine coordinates.
    Starting coordinates may be random of from a previous (e.g. lower
    resolution) stage.
    """
    
    from numpy import random
    from math import log, exp, atan, pi
    import gc
    
    random.seed(general_calc_params['random_seed'])
    particle_size = np.int32(particle_size)
    
    # Calculate distance restrains from contact data   
    restraint_dict, seq_pos_dict = dyn_util.calc_restraints(chromosomes, contact_dict, particle_size,
                                                   scale=1.0, exponent=general_calc_params['dist_power_law'],
                                                   lower=general_calc_params['contact_dist_lower'],
                                                   upper=general_calc_params['contact_dist_upper'])
    
    # Concatenate chromosomal data into a single array of particle restraints
    # for structure calculation. Add backbone restraints between seq. adjasent particles.
    restraint_indices, restraint_dists = dyn_util.concatenate_restraints(restraint_dict, seq_pos_dict, particle_size,
                                                                general_calc_params['backbone_dist_lower'],
                                                                general_calc_params['backbone_dist_upper'])
 
    # Setup starting structure
    if (start_coords is None) or (prev_seq_pos_dict is None):
      coords = get_random_coords(seq_pos_dict, chromosomes, num_models,
                                 general_calc_params['random_radius'])
      
      num_coords = coords.shape[1]
            
    else:
      # Convert starting coord dict into single array
      coords = pack_chromo_coords(start_coords, chromosomes)
      num_coords = sum([len(seq_pos_dict[c]) for c in chromosomes])
        
      if coords.shape[1] != num_coords: # Change of particle_size
        interp_coords = np.empty((num_models, num_coords, 3))
        
        for m in range(num_models): # Starting coords interpolated from previous particle positions
          interp_coords[m] = dyn_util.getInterpolatedCoords(chromosomes, seq_pos_dict, prev_seq_pos_dict, coords[m])
        
        coords = interp_coords
        
    # Equal unit masses and radii for all particles
    masses = np.ones(num_coords,  float)
    radii = np.ones(num_coords,  float)
    
    # Ambiguiity strides not used here, so set to 1
    num_restraints = len(restraint_indices)
    ambiguity = np.ones(num_restraints,  np.int32)
        
    # Below will be set to restrict memory allocation in the repusion list
    # (otherwise all vs all can be huge)
    n_rep_max = np.int32(0)
    
    # Annealing parameters
    temp_start = anneal_params['temp_start']
    temp_end = anneal_params['temp_end']
    temp_steps = anneal_params['temp_steps']
    
    # Setup annealig schedule: setup temps and repulsive terms
    adj = 1.0 / atan(10.0)
    decay = log(temp_start/temp_end)    
    anneal_schedule = []
    
    for step in range(temp_steps):
      frac = step/float(temp_steps)
    
      # exponential temp decay
      temp = temp_start * exp(-decay*frac)
    
      # sigmoidal repusion scheme
      repulse = 0.5 + adj * atan(frac*20.0-10) / pi 
      
      anneal_schedule.append((temp, repulse))  
        
    # Paricle dynamics parameters
    # (these need not be fixed for all stages, but are for simplicity)    
    dyn_steps = anneal_params['dynamics_steps']
    time_step = anneal_params['time_step']
       
    # Update coordinates in the annealing schedule
    time_taken = 0.0
    
    for m in range(num_models): # For each repeat calculation
      model_coords = coords[m]
        
      for temp, repulse in anneal_schedule:
        gc.collect() # Try to free some memory
        
        # Update coordinates for this temp
        dt, n_rep_max = dyn_util.runDynamics(model_coords, masses, radii, restraint_indices, restraint_dists,
                                    ambiguity, temp, time_step, dyn_steps, repulse, nRepMax=n_rep_max)
      
        n_rep_max = np.int32(1.05 * n_rep_max) # Base on num in prev cycle, plus an overhead
        time_taken += dt
  
      # Center
      model_coords -= model_coords.mean(axis=0)
      
      # Update
      coords[m] = model_coords
    
    # Convert from single coord array to dict keyed by chromosome
    coords_dict = unpack_chromo_coords(coords, chromosomes, seq_pos_dict)
    
    return coords_dict, seq_pos_dict


def calc_genome_structure(ncc_file_path, pdb_file_path, general_calc_params, anneal_params,
                          particle_sizes, num_models=5, isolation_threshold=2e6):

  from time import time

  # Load single-cell Hi-C data from NCC contact file, as output from NucProcess
  chromosomes, chromo_limits, contact_dict = load_ncc_file(ncc_file_path)

  # Only use contacts which are supported by others nearby in sequence, in the initial instance
  remove_isolated_contacts(contact_dict, threshold=isolation_threshold)

  # Initial coords will be random
  start_coords = None

  # Record partile positions from previous stages
  # so that coordinates can be interpolated to higher resolution
  prev_seq_pos = None

  for stage, particle_size in enumerate(particle_sizes):
 
      print("Running structure caculation stage %d (%d kb)" % (stage+1, (particle_size/1e3)))
 
      # Can remove large violations (noise contacts inconsistent with structure)
      # once we have a resonable resolution structure
      
      if stage > 0:
        if particle_size < 0.5e6:
            remove_violated_contacts(contact_dict, coords_dict, particle_seq_pos,
                                     particle_size, threshold=6.0)
        elif particle_size < 0.25e6:
            remove_violated_contacts(contact_dict, coords_dict, particle_seq_pos,
                                     particle_size, threshold=5.0)
 
      coords_dict, particle_seq_pos = anneal_genome(chromosomes, contact_dict, num_models, particle_size,
                                                    general_calc_params, anneal_params,
                                                    prev_seq_pos, start_coords)
 
      # Next stage based on previous stage's 3D coords
      # and thier respective seq. positions
      start_coords = coords_dict
      prev_seq_pos = particle_seq_pos

  # Save final coords as PDB format file
  export_pdb_coords(pdb_file_path, coords_dict, particle_seq_pos, particle_size)
  print('Saved structure file to: %s' % pdb_file_path)


PROG_NAME = 'nuc_dynamics'
DESCRIPTION = 'Single-cell Hi-C genome and chromosome structure calculation module for Nuc3D and NucTools'

def warn(msg, prefix='WARNING'):

  print('%8s : %s' % (prefix, msg))


def test_imports(gui=False):
  import sys
  from distutils.core import run_setup
  
  critical = False
  
  try:
    import numpy
  except ImportError as err:
    critical = True
    warn('Critical Python module "numpy" is not installed or accessible')

  try:
    import cython
  except ImportError as err:
    critical = True
    warn('Critical Python module "cython" is not installed or accessible')
  
  try:
    import dyn_util
  except ImportError as err:
    warn('Utility C/Cython code not compiled. Attempting to compile now...')    
    run_setup('setup_cython.py', ['build_ext', '--inplace'])
    
    try:
      import dyn_util
    except ImportError as err:
      critical = True
      warn('Utility C/Cython code compilation/import failed')   
    
  if critical:
    warn('NucDynamics cannot proceed because critical Python modules are not available', 'ABORT')
    sys.exit(0)
    
    
def demo_calc_genome_structure():
  """
  Example of settings for a typical genome structure calculation from input single-cell
  Hi-C contacts stored in an NCC format file (as output from NucProcess)
  """
  
  from nuc_dynamics import calc_genome_structure
  
  ncc_file_path = 'example_chromo_data/Cell_1_contacts.ncc'
  save_path = 'example_chromo_data/Cell_1_structure.pdb'
  
  # Number of alternative conformations to generate from repeat calculations
  # with different random starting coordinates
  num_models = 2

  # Parameters to setup restraints and starting coords
  general_calc_params = {'dist_power_law':-0.33,
                         'contact_dist_lower':0.8, 'contact_dist_upper':1.2,
                         'backbone_dist_lower':0.1, 'backbone_dist_upper':1.1,
                         'random_seed':int(time()), 'random_radius':10.0}

  # Annealing & dyamics parameters: the same for all stages
  # (this is cautious, but not an absolute requirement)
  anneal_params = {'temp_start':5000.0, 'temp_end':10.0, 'temp_steps':500,
                   'dynamics_steps':100, 'time_step':0.001}

  # Hierarchical scale protocol
  particle_sizes = [8e6, 4e6, 2e6, 4e5, 2e5, 1e5]
  
  # Contacts must be clustered with another within this separation threshold
  # (at both ends) to be considered supported, i.e. not isolated
  isolation_threshold=2e6
  
  calc_genome_structure(ncc_file_path, pdb_file_path, general_calc_params, anneal_params,
                        particle_sizes, num_models, isolation_threshold)

    
test_imports()

if __name__ == '__main__':
  
  import os, sys
  from time import time
  from argparse import ArgumentParser
  from nuc_dynamics import calc_genome_structure
  
  epilog = 'For further help on running this program please email tjs23@cam.ac.uk'
  
  arg_parse = ArgumentParser(prog=PROG_NAME, description=DESCRIPTION,
                             epilog=epilog, prefix_chars='-', add_help=True)

  arg_parse.add_argument('ncc_path', nargs=1, metavar='NCC_FILE',
                         help='Input NCC format file containing single-cell Hi-C contact data, e.g. use the demo data at example_chromo_data/Cell_1_contacts.ncc')

  arg_parse.add_argument('-o', metavar='PDB_FILE',
                         help='Optional name of PDB format output file for 3D coordinates. If not set this will be auto-generated from the input file name')

  arg_parse.add_argument('-m', default=1, metavar='NUM_MODELS',
                         type=int, help='Number of alternative conformations to generate from repeat calculations with different random starting coordinates: Default: 1')

  arg_parse.add_argument('-s', nargs='+', default=[8.0,4.0,2.0,0.4,0.2,0.1], metavar='Mb_SIZE', type=float,
                         help='One or more sizes (Mb) for the hierarchical structure calculation protocol (will be used in descending order). Default: 8.0 4.0 2.0 0.4 0.2 0.1')

  arg_parse.add_argument('-iso', default=2.0, metavar='Mb_SIZE', type=float,
                         help='Contacts must be near another, within this (Mb) separation threshold (at both ends) to be considered supported: Default 2.0')

  arg_parse.add_argument('-pow', default=-0.33, metavar='FLOAT',
                         type=float, help='Distance power law for combining multiple Hi-C contacts between the same particles. Default: -0.33')

  arg_parse.add_argument('-lower', default=0.8, metavar='DISTANCE',
                         type=float, help='Lower limit for a contact-derived distance restraint, as a fraction of the ideal distance. Default: 0.8')

  arg_parse.add_argument('-upper', default=1.2, metavar='DISTANCE',
                         type=float, help='Upper limit for a contact-derived distance restraint, as a fraction of the ideal distance. Default: 1.2')

  arg_parse.add_argument('-bb_lower', default=0.1, metavar='DISTANCE',
                         type=float, help='Lower limit for sequential particle backbone restraints, as a fraction of the ideal distance. Default: 0.1')

  arg_parse.add_argument('-bb_upper', default=1.1, metavar='DISTANCE',
                         type=float, help='Upper limit for sequential particle backbone restraints, as a fraction of the ideal distance. Default: 1.1')

  arg_parse.add_argument('-ran', metavar='INT',
                         type=int, help='Seed for psuedo-random number generator')

  arg_parse.add_argument('-rad', default=10.0, metavar='DISTANCE',
                         type=float, help='Radius of sphere for random starting coordinates. Default: 10.0')

  arg_parse.add_argument('-hot', default=5000.0, metavar='TEMP_KELVIN',
                         type=float, help='Start annealing temperature in pseudo-Kelvin units. Default: 5000')

  arg_parse.add_argument('-cold', default=10.0, metavar='TEMP_KELVIN',
                         type=float, help='End annealing temperature in pseudo-Kelvin units. Default: 10')

  arg_parse.add_argument('-temps', default=500, metavar='NUM_STEPS',
                         type=int, help='Number of temperature steps in annealing protocol between start and end temperature. Default: 500')

  arg_parse.add_argument('-dyns', default=100, metavar='NUM_STEPS',
                         type=int, help='Number of particle dynamics steps to apply at each temperature in the annealing protocol. Default: 100')

  arg_parse.add_argument('-time_step', default=0.001, metavar='TIME_DELTA',
                         type=float, help='Simulation time step between re-calculation of particle velocities. Default: 0.001')
  
  args = vars(arg_parse.parse_args())
  
  ncc_file_path = args['ncc_path'][0]
  
  save_path = args['o']
  if save_path is None:
    save_path = os.path.splitext(ncc_file_path)[0] + '.pdb'
    
  particle_sizes = args['s']
  particle_sizes = sorted([x * 1e6 for x in particle_sizes if x > 0], reverse=True)
  if not particle_sizes:
    warn('No positive particle sizes (Mb) specified', 'ABORT')
    sys.exit(0)
  
  num_models = args['m']
  dist_power_law = args['pow']
  contact_dist_lower = args['lower']
  contact_dist_upper = args['upper']
  backbone_dist_lower = args['bb_lower']
  backbone_dist_upper = args['bb_upper']
  random_radius = args['rad']
  random_seed = args['ran']
  temp_start = args['hot']
  temp_end = args['cold']
  temp_steps = args['temps']
  dynamics_steps = args['dyns']
  time_step = args['time_step']
  isolation_threshold = args['iso']
  
  for val, name, sign in ((num_models,         'Number of conformational models', '+'),
                          (dist_power_law,     'Distance power law', '-0'),
                          (contact_dist_lower, 'Contact distance lower bound', '+'),
                          (contact_dist_upper, 'Contact distance upper bound', '+'),
                          (backbone_dist_lower,'Backbone distance lower bound', '+'),
                          (backbone_dist_upper,'Backbone distance upper bound', '+'),
                          (random_radius,      'Random-start radius', '+'),
                          (temp_start,         'Annealing start temperature', '+'),
                          (temp_end,           'Annealing end temperature', '+0'),
                          (temp_steps,         'Number of annealing temperature steps', '+'),
                          (dynamics_steps,     'Number of particle dynamics steps', '+'),
                          (time_step,          'Particle dynamics time steps', '+'),
                          (isolation_threshold,'Contact isolation threshold', '+0'),
                         ):
                          
    if '+' in sign:
      if '0' in sign:
        if val < 0.0:
          warn('%s must be non-negative' % name, 'ABORT')
          sys.exit(0)
      
      else:
        if val <= 0.0:
          warn('%s must be positive' % name, 'ABORT')
          sys.exit(0)
      
    elif '-' in sign:  
      if '0' in sign:
        if val > 0.0:
          warn('%s must be non-positive' % name, 'ABORT')
          sys.exit(0)
      
      else:
        if val >= 0.0:
          warn('%s must be negative' % name, 'ABORT')
          sys.exit(0)
     

  contact_dist_lower, contact_dist_upper = sorted([contact_dist_lower, contact_dist_upper])
  backbone_dist_lower, backbone_dist_upper = sorted([backbone_dist_lower, backbone_dist_upper])
  temp_end, temp_start = sorted([temp_end, temp_start])
  
  if not random_seed:
    random_seed = int(time())
  
  general_calc_params = {'dist_power_law':dist_power_law,
                         'contact_dist_lower':contact_dist_lower,
                         'contact_dist_upper':contact_dist_upper,
                         'backbone_dist_lower':backbone_dist_lower,
                         'backbone_dist_upper':backbone_dist_upper,
                         'random_seed':random_seed,
                         'random_radius':random_radius}

  anneal_params = {'temp_start':temp_start, 'temp_end':temp_end, 'temp_steps':temp_steps,
                   'dynamics_steps':dynamics_steps, 'time_step':time_step}

  
  isolation_threshold *= 1e6
  
  calc_genome_structure(ncc_file_path, save_path, general_calc_params, anneal_params,
                        particle_sizes, num_models, isolation_threshold)

# TO DO
# -----
# Allow chromosomes to be specified
# Add an official NucDynamics struture output format
# Allow starting structures to be input


