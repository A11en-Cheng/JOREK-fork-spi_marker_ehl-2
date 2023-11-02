# ------------------------------------------------ #
# Python script to be used for comparing the       #
# histograms of two JOREK particle populations.    #
# The JOREK particles should be stored in two      #
# different HDF5 files.                            #
# ------------------------------------------------ #

# Convert JOREK particle phase space data into table
def convert_jorek_particle_velocities_charges_to_table(group):
  from numpy import array,transpose,zeros,float64
  velocity = array([]);  charge = array([]);
  if('q' in group):
    charge = transpose(array(group['q']))
  if('v' in group):
    velocity = transpose(array(group['v']))
  elif('mu' in group):
    if('Vpar' in group):
      dummy_array = array(group['Vpar'])
    elif('E' in group):
      dummy_array = array(group['E'])
    mu       = array(group['mu'])
    velocity = zeros((2,mu.size),dtype=float64)
    velocity[0,:] = dummy_array; velocity[1,:] = mu;   
  return velocity,charge

# Read a JOREK particle restart file. Zero ended
# bytes structure ('S'-type) are identified and
# transformed in strings. The JOREK particle data 
# are homogenized in a N-dimensional table.
def read_jorek_particle_restart_file(filename,filepath,separator):
  from h5py import File
  from numpy import array,float64,transpose,zeros,floor_divide,pi
  particle_tables = []; mass = []; weights = [];
  fhandler = File("".join([filepath,separator,filename]))
  for group in fhandler['groups'].values():
    positions          = transpose(array(group['x']))
    velocities,charges = convert_jorek_particle_velocities_charges_to_table(group)
    ndim = positions.shape[0]
    if(velocities.shape[0]==positions.shape[1]):
      ndim = ndim+1
    else:
      ndim = ndim + velocities.shape[0]
    if(charges.size>0):
      ndim = ndim + 1
    particles = zeros((ndim,positions.shape[1]),dtype=float64)
    particles[0:positions.shape[0],:] = positions
    if(velocities.size>0):
      particles[positions.shape[0]:positions.shape[0]+velocities.shape[0],:] = velocities
    if(charges.size>0):
      particles[-1,:] = charges
    particle_tables.append(particles)
    mass.append(float64(group['mass']))
    weights.append(array(group['weight']))
  sim_time = (float64(fhandler['time']))
  fhandler.close()
  return particle_tables,weights,array(mass),sim_time

# compute the range of coordinates for each dimensions of each particle table
def compute_particle_coordinate_ranges(particle_table):
  from numpy import amin,amax,array
  min_coord = amin(particle_table,axis=1)
  max_coord = amax(particle_table,axis=1)
  coordinate_ranges = array([array(min_coord),array(max_coord)])
  return coordinate_ranges

# assemble common coordinate ranges for both test and reference histograms
# requires the test and reference particle simulations to have the same
# number of particle tables
def assemble_common_coordinate_ranges(test_tables,ref_tables):
  from numpy import amax,amin,array,transpose
  common_coordintate_ranges = []
  for test_table_id,test_table in enumerate(test_tables):
    test_coordinate_ranges = compute_particle_coordinate_ranges(test_table)
    ref_coordinate_ranges = compute_particle_coordinate_ranges(ref_tables[test_table_id])
    lower_bound = amin(array([test_coordinate_ranges[0,:],ref_coordinate_ranges[0,:]]),axis=0)
    upper_bound = amax(array([test_coordinate_ranges[1,:],ref_coordinate_ranges[1,:]]),axis=0)
    common_coordintate_ranges.append(array([lower_bound,upper_bound]))
  common_coordintate_ranges = transpose(array(common_coordintate_ranges),axes=[0,2,1])
  print('histogram range: ',common_coordintate_ranges)
  return common_coordintate_ranges

# compute the multivariate histograms for all particle tables    
def compute_particle_histograms(particle_tables,weights,coordinate_ranges,nbins):
  from numpy import histogramdd,transpose
  histograms = []; edges = [];
  for particle_id,particles in enumerate(particle_tables):
    histogram,edge = histogramdd(transpose(particles),bins=nbins[0:particles.shape[0]],\
    range=coordinate_ranges[particle_id],density=None,weights=weights[particle_id])
    histograms.append(histogram); edges.append(edge);
  return histograms,edges

# compute L error normalised to the number of bins of a given order
def compute_L_error(diff_histogram,order=2,axis=None,bins=[]):
  from numpy import multiply
  from numpy import sum as npsum
  from numpy import abs as npabs
  from numpy import sqrt as npsqrt
  if(order==1):
    L_error = npsum(npabs(diff_histogram),axis=axis,dtype=diff_histogram.dtype,\
    initial=0e0)
  elif(order==2):
    L_error = npsqrt(npsum(multiply(diff_histogram,diff_histogram),axis=axis,\
    dtype=diff_histogram.dtype,initial=0e0))
  else:
    raise ValueError("unrecognised L norm order for computing the error")
  if(axis==None):
    L_error = L_error/diff_histogram.size
  elif(isinstance(axis,int)):
    L_error = L_error/bins[axis]
  else:
    L_error = L_error/prod(bins[axis])
  return L_error 

# compute L-order error between histograms. The number of particle tables between
# test and reference simulations must be the same
def compute_histogram_L_error(test_histograms,ref_histograms,orders=[1,2],axis=None,bins=[]):
  from numpy import array
  L_errors_glob = []
  for histogram_id,histogram in enumerate(test_histograms):
    diff_histogram = histogram-ref_histograms[histogram_id]; L_errors = [];
    for order in orders:
      L_errors.append(compute_L_error(diff_histogram,order=order,axis=axis,bins=bins))
    L_errors_glob.append(array(L_errors))
  return array(L_errors_glob)

# log overall histogram L errors
def log_histogram_L_errors(histogram_L_errors,orders):
  for L_errors_id,L_errors in enumerate(histogram_L_errors):
    for L_error_id,L_error in enumerate(L_errors):
      print("".join(['Particle table N#: ',str(L_errors_id),' histogram L-',\
      str(orders[L_error_id]),' error: ',str(L_error)]))

# log the histogram L errors
def show_histogram_errors(histogram_L_errors,histogram_edges,orders):
  if(histogram_L_errors.ndim==2):
    log_histogram_L_errors(histogram_L_errors,orders)

# main function
def compute_error_between_histograms(test_filename,ref_filename,filepath,separator,\
nbins,L_error_orders):
  # read particle tables
  test_tables,test_weights,test_mass,test_sim_time = read_jorek_particle_restart_file(\
  test_filename,filepath,separator)
  ref_tables,ref_weights,ref_mass,ref_sim_time = read_jorek_particle_restart_file(\
  ref_filename,filepath,separator)
  # compute coordinate ranges
  coordinate_ranges = assemble_common_coordinate_ranges(test_tables,ref_tables)
  # compute histograms
  test_histograms,test_edges = compute_particle_histograms(test_tables,test_weights,\
  coordinate_ranges,nbins)
  ref_histograms,ref_edges   = compute_particle_histograms(ref_tables,ref_weights,\
  coordinate_ranges,nbins)
  # clean up the particle tables for saving memory
  del test_tables,test_weights,ref_tables,ref_weights
  # compute histograms overall error 
  L_errors_all = compute_histogram_L_error(test_histograms,ref_histograms,\
  orders=L_error_orders,axis=None,bins=[])
  # log L-norm errors
  show_histogram_errors(L_errors_all,test_edges,L_error_orders)

# argument parser
def generate_argument_parser():
  from argparse import ArgumentParser
  parser = ArgumentParser(description='read and compute metrics for statistical coherency of particle histograms')
  parser.add_argument('--filename','-f',type=str,action='store',required=False,\
  dest='filename',default='part_restart.h5',\
  help='name of the particle restart file, default: part_restart.h5')
  parser.add_argument('--ref-filename','-rf',type=str,action='store',required=False,\
  dest='ref_filename',default='part_restart_original.h5',\
  help='name of the reference particle restart file, default: part_restart_original.h5') 
  parser.add_argument('--filepath','-fpath',type=str,action='store',required=False,\
  dest='filepath',default='.',help='path of the file to be read, default: .')
  parser.add_argument('--separator','-sep',type=str,action='store',required=False,\
  dest='separator',default='/',help='file separator, default: /')
  parser.add_argument('--nbins','-nb',type=int,nargs='*',action='store',required=False,\
  default=[10,10,10,10,10,10,10],dest='nbins',help='number of bins for each particle coordinate, default: 10')
  parser.add_argument('--Lorder','-lo',type=int,nargs='*',action='store',required=False,\
  default=[1,2],dest='L_error_orders',help='L norm orders for computing histogram errors, default: [1,2]')
  return parser.parse_args()

# Run main --------------------------------------- #
if __name__ == "__main__":
  args = generate_argument_parser()
  compute_error_between_histograms(args.filename,args.ref_filename,args.filepath,\
  args.separator,args.nbins,args.L_error_orders)

# ------------------------------------------------ #
