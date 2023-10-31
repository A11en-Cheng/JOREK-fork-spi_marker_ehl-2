# ----------------------------------------------------- #
# The function plots the particle properties read by
# a JOREK restart file
# ----------------------------------------------------- #
# Transform cylindrical in cartesian coordinates
def cylindrical_to_cartesian(RZphi):
  from numpy import cos,sin
  return np.array([RZphi[0]*cos(-RZphi[2]),RZphip[0]*sin(-RZphi[2]),RZphi[1]])

# compute the radial volumes for cylindrical coordinates
def compute_cylindrical_radial_volumes(r):
  from numpy import power
  return 5e-1*(power(r[2:],2)-power(r[1:-2],2))
# compute cartesian volumes
def computes_cartesian_volumes(x):
  return x[2,:]-x[1:-2]
# compute the radial volumes for spherical coordinates
def compute_spherical_radial_volumes(r):
  from numpy import power
  return (power(r[2:],3)-power(r[1:-2],3))/3e0
# compute spherical azimuthal volume
def compute_spherical_azimuthal_volumes(theta):
  from numpy import cos
  return cos(theta[1:-2])-cos(theta[2:])

# Return the correct x labels and titles in phase space
def define_histogram_labels_titles(key,ptype):
  if(key=='x'):
    titles = ['Major radius','Vertical coordinate','Toroidal angle']
    labels = ['R [m]','Z [m]','phi [r]']
  elif(key=='v'):
    if(ptype=='particle_kinetic_relativistic'):
      titles = ['x-cartesian momentum','y-cartesian momentum','z-cartesian momentum']
      labels = ['px [AMUm/s]','py [AMUm/s]','pz [AMUm/s]']
    elif(ptype=='particle_gc_relativistic'):
      titles = ['Parallel momentum','Magnetic moment']
      labels = ['p_par [AMUm/s]','mu [AMUm^2/Cs]']
    elif(ptype=='particle_gc'):
      titles = ['Energy','Magnetic moment']
      labels = ['E [eV]','mu [m^2/s]']
    else:
      titles = ['velocity 1','velocity 2','velocity 3']
      labels = ['v1','v2','v3']
  return titles,labels

# Return the correct x,y labels, titles and aspect ratio in phase space
def define_histogram_labels_titles_aspectratio_2d(key,ptype):
  if(key=='x'):
    titles  = ['Major radius - vertical coordinate','Major radius - toroidal angle',\
    'Vertical coordinate - toroidal angle']
    xlabels = ['R [m]','R [m]','Z [m]']
    ylabels = ['Z [m]','phi [r]','phi [r]']
    aspectratio = [True,False,False]
  elif(key=='v'):
    if(ptype=='particle_kinetic_relativistic'):
      titles  = ['Momenta: px-py','Momenta: px-pz','Momenta: py-pz']
      xlabels = ['px [AMUm/s]','px [AMUm/s]','py [AMUm/s]']
      ylabels = ['py [AMUm/s]','pz [AMUm/s]','pz [AMUm/s]']
      aspectratio = [True,True,True]
    elif(ptype=='particle_gc_relativistic'):
      titles  = ['Parallel momentum - magnetic moment']
      xlabels = ['p_par [AMUm/s]']
      ylabels = ['mu [AMUm^2/Cs]']
      aspectratio = [False]
    elif(ptype=='particle_gc'):
      titles   = ['Energy - magnetic moment']
      xlabels = ['E [eV]']
      ylabels = ['mu [m^2/s]']
      aspectratio = [False]
    else:
      titles  = ['velocities v1-v2','velocities v1-v3','velocities v2-v3']
      xlabels = ['v1','v1','v2']
      ylabels = ['v2','v3','v3']
      aspectratio = [False,False,False]
  return titles,xlabels,ylabels,aspectratio

# return same velocity structure as a function of the particle type
def homogenize_particle_structure(groups):
  from numpy import zeros
  for group_id,group in enumerate(groups):
    if(group['type']=='particle_gc'):
      v = zeros((2,group['E'].size),dtype=group['E'].dtype)
      v[0,:] = group['E']; v[1,:] = group['mu']
      group['v'] = v; del group['E'], group['mu']
      groups[group_id] = group
  return groups

# Read a jorek particle restart file. Zero ended bytes structre 
# ('S'-type) are identified and transformed in strings
def read_jorek_particle_restart_file(filename,filepath,separator):
  import h5py
  from numpy import array,float64,transpose,floor_divide,pi
  groups = []
  fhandler = h5py.File("".join([filepath,separator,filename]))
  for group in fhandler['groups'].values():
    particles = {}
    for part_data_name,part_data in group.items():
      particles[part_data_name] = array(part_data)
      if((part_data_name=='x') or (part_data_name=='v')):
        particles[part_data_name] = transpose(particles[part_data_name])
      if(part_data_name=='x'):
        particles[part_data_name][2,:] = particles[part_data_name][2,:]-\
        2e0*pi*floor_divide(particles[part_data_name][2,:],2e0*pi)
      if('S' in str(particles[part_data_name].dtype)):
        particles[part_data_name] = str(particles[part_data_name])
        particles[part_data_name] = particles[part_data_name][2:-1]
    groups.append(particles)
  sim_time = float64(fhandler['time'])
  fhandler.close()
  return groups,sim_time

# compute the amount of uninitialised particles w.r.t. the total number of particles
def compute_uninitialised_particles(groups,deathflags,deathvalues):
  for group_id,group in enumerate(groups):
    for deathflag_id,deathflag in enumerate(deathflags):
      values = group[deathflag]
      n_dead_particles = len(values[values<=deathvalues[deathflag_id]]) 
      print("group id: ",group_id," number of particles with ",\
      deathflag,"<=",deathvalues[deathflag_id],1e2*float(n_dead_particles)/len(values),"%, ",n_dead_particles)

# generate 1d histograms for a set of positions (physical or velocity space)
# results are stored in a list having elements of the form: 
# [histogram,histogram_edges]
def create_phase_space_1d_histograms(data_array,p_weights,flags,deathvalue,bins=[]):
  from numpy import histogram
  hists = []
  for ids,data in enumerate(data_array):
    histo,edges = histogram(data[flags>deathvalue],bins=bins[ids],weights=p_weights[flags>deathvalue])
    hists.append([histo,edges])
  return hists

# generate 2d histograms for a set of positions (physical or velocity space)
# results are stored in a list of the form 
def create_phase_space_2d_histograms(data_array,p_weights,flags,deathvalue,bins2d=[]):
  from numpy import histogram2d,amin,amax
  n_histograms = 0
  histos = []
  for id1,data1 in enumerate(data_array):
    local_histos = []
    for id2,data2 in enumerate(data_array[id1+1:]):
      histo,xedges,yedges = histogram2d(data1[flags>deathvalue],data2[flags>deathvalue],\
      bins=[bins2d[id1],bins2d[id1+id2+1]],weights=p_weights[flags>deathvalue])
      local_histos.append([histo,xedges,yedges])
    n_histograms = n_histograms + len(local_histos)
    histos.append(local_histos)
  return histos,n_histograms

# generate 2d histograms for a set of positions slicing at a specific
# interval of one data id
def create_phase_space_2d_sliced_histograms(data_array,p_weights,slice_data_ids,\
slice_data_intervals,flags,deathvalues,bins2d=[]):
  from numpy import histogram2d,amin,amax
  n_histograms = 0
  histos = []
  for slice_id,slice_data_id in enumerate(slice_data_ids):
    slice_mask = [(data>deathvalues and data >= slice_data_intervals[slice_id,0] and \
    data < slice_data_intervals[slice_id,1]) for data in data_array[slice_data_id]]
    compute_ids = [ids for ids,data in enumerate(data_array)]
    compute_ids.remove(slice_data_id)
    for id1,data1_id in  enumerate(compute_ids):
      local_histos = []
      for data2_id in compute_ids[id1+1:]:
        histo,xedges,yedges = histogram2d(data_array[data1_id,slice_mask],\
        data_array[data2_id,slice_mask],bins=[bins2d[data1_id],bins2d[data2_id]],\
        weights=p_weights[slice_mask])
        local_histos.append([histo,xedges,yedges])
      n_histograms = n_histograms + len(local_histos)
      histos.append(local_histos)
  return histos,n_histograms 

# plot 1d histograms
def plot_1d_histograms(hists,titles,xlabels,ylabels,fontsize=18):
  from matplotlib.pyplot import figure,stairs
  figs = []; axs = [];
  for histo_id,histo in enumerate(hists):
    figs.append(figure(histo_id+1,facecolor='white',edgecolor='white'))
    axs.append([figs[histo_id].add_subplot(111)])
    axs[histo_id][0].stairs(histo[0],edges=histo[1],fill=True)
    if(len(titles)>histo_id):
      axs[histo_id][0].set_title(titles[histo_id],fontsize=fontsize)
    if(len(xlabels)>histo_id):
      axs[histo_id][0].set_xlabel(xlabels[histo_id],fontsize=fontsize)
    if(len(ylabels)>histo_id):
      axs[histo_id][0].set_ylabel(ylabels[histo_id],fontsize=fontsize)
    axs[histo_id][0].tick_params(axis='x',labelsize=fontsize)
    axs[histo_id][0].tick_params(axis='y',labelsize=fontsize)
    axs[histo_id][0].grid()
  return figs,axs

# plot 2d histograms
def plot_2d_histograms(hists2d,n_histos,titles,xlabels,ylabels,aspectequal,\
fontsize=18,colormap='inferno',colorbarmaxval=[1.,1.,1.,1.]):
  from numpy import amax,transpose
  from matplotlib.pyplot import figure,pcolormesh
  count = 0; figs = []; axs = [];
  for histos_id,histos in enumerate(hists2d):
    for histo_id,histo in enumerate(histos):
      scaling = 1e0;
      if(count<len(colorbarmaxval)):
        scaling = colorbarmaxval[count]
      figs.append(figure(count+1,facecolor='white',edgecolor='white'))
      axs.append([figs[count].add_subplot(111)])
      im = axs[count][0].pcolormesh(histo[1],histo[2],transpose(histo[0]),\
      cmap=colormap,vmin=0.,vmax=scaling*amax(histo[0]),edgecolors='none',shading='flat')
      if(len(titles)>count):
        axs[count][0].set_title(titles[count],fontsize=fontsize)
      if(len(xlabels)>count):
        axs[count][0].set_xlabel(xlabels[count],fontsize=fontsize)
      if(len(ylabels)>count):
        axs[count][0].set_ylabel(ylabels[count],fontsize=fontsize)
      axs[count][0].tick_params(axis='x',labelsize=fontsize)
      axs[count][0].tick_params(axis='y',labelsize=fontsize)
      if(len(aspectequal)>count):
        if(aspectequal[count]):
          axs[count][0].set_aspect('equal',adjustable='datalim')
      figs[count].colorbar(im,ax=axs[count][0])
      count = count + 1
  return figs,axs
      
# generate and plot 1d histograms from the jorek particle distribution
def computes_1d_histograms_jorek_particles(groups,key,deathvalue,bins=[100,100,100],\
ylabels=['','',''],fontsize=18):
  from matplotlib.pyplot import show
  for group in groups: 
    # extract titles and x lables
    titles,xlabels = define_histogram_labels_titles(key,group['type'])
    # compute histograms physical space and plot it
    hists = create_phase_space_1d_histograms(group[key],group['weight'],group['i_elm'],\
    deathvalue,bins=bins)
    figs,axs = plot_1d_histograms(hists,titles,xlabels,ylabels,fontsize=fontsize)
    show()

# generate and plot 2d histograms from the jorek particle distribution
def computes_2d_histograms_jorek_particles(groups,key,deathvalue=0,bins=[100,100,100],\
fontsize=18,colormap='inferno',colorbarmaxval=[1.,1.,1.,1.]):
  from matplotlib.pyplot import show
  for group in groups:
    # extract titles, xlabels, ylabels and plot aspect ratio
    titles,xlabels,ylabels,aspectratio = define_histogram_labels_titles_aspectratio_2d(key,group['type'])
    # compute histograms physical space and plot it
    hists2d,n_histos = create_phase_space_2d_histograms(group[key],group['weight'],group['i_elm'],\
    deathvalue,bins2d=bins)
    fig,axs = plot_2d_histograms(hists2d,n_histos,titles,xlabels,ylabels,aspectratio,\
    fontsize=fontsize,colormap=colormap,colorbarmaxval=colorbarmaxval)
    show()

# generate and plot 2d histogrames from jorek particle distribution
# give a set of coordinate slices
def computes_phase_space_2d_sliced_histograms(groups,key,slice_data_ids=[],\
slice_data_intervals=[],deathvalue=0,bins=[100,100,100],\
fontsize=18,colormap='inferno',colorbarmaxval=[1.,1.,1.,1.]):
  from numpy import array,reshape
  from matplotlib.pyplot import show
  slice_data_intervals = reshape(slice_data_intervals,(len(slice_data_ids),2))
  for group in groups:
    # extract titles, xlabels, ylabels and plot aspect ratio
    if(slice_data_ids==[2,1,0]):
      titles,xlabels,ylabels,aspectratio = define_histogram_labels_titles_aspectratio_2d(key,group['type']) 
    else:
      titles = []; xlabels = []; ylabels = []; aspectratio = [True];
    # compute histograms physical space slacing values and plot it
    hists2d,n_histos = create_phase_space_2d_sliced_histograms(group[key],group['weight'],\
    array(slice_data_ids),slice_data_intervals,group['i_elm'],deathvalue,bins2d=bins)
    fig,axs = plot_2d_histograms(hists2d,n_histos,titles,xlabels,ylabels,aspectratio,\
    fontsize=fontsize,colormap=colormap,colorbarmaxval=colorbarmaxval)
    show()

# main function
def read_analyse_plot_jorek_restart(filename,filepath,separator,deathvalue=0,bins1d=[100,100,100],
bins2d=[100,100,100],fontsize=18,colormap='inferno',spatial_slice_ids=[],\
spatial_slice_intervals=[],velocity_slice_ids=[],velocity_slice_intervals=[],\
colorbarmaxval=[1.,1.,1.,1.]):
  # read the jorek particle restart data
  p_groups,sim_time = read_jorek_particle_restart_file(filename,filepath,separator)
  p_groups = homogenize_particle_structure(p_groups)
  # 0d analysis
  compute_uninitialised_particles(p_groups,['i_elm','weight'],[0,0.])
  # plot 1d histograms
  computes_1d_histograms_jorek_particles(p_groups,'x',deathvalue,bins=bins1d,\
  ylabels=['Nphys','Nphys','Nphys'],fontsize=fontsize)
  computes_1d_histograms_jorek_particles(p_groups,'v',deathvalue,bins=bins1d,\
  ylabels=['Nphys','Nphys','Nphys'],fontsize=fontsize)
  # plot 2d histograms
  computes_2d_histograms_jorek_particles(p_groups,'x',deathvalue,bins=bins2d,\
  fontsize=fontsize,colormap=colormap,colorbarmaxval=colorbarmaxval)
  computes_2d_histograms_jorek_particles(p_groups,'v',deathvalue,bins=bins2d,\
  fontsize=fontsize,colormap=colormap,colorbarmaxval=colorbarmaxval)
  # plot 2d histograms with slicing
  if(len(spatial_slice_ids)!=0 and len(spatial_slice_intervals)==2*len(spatial_slice_ids)):
    computes_phase_space_2d_sliced_histograms(p_groups,'x',slice_data_ids=spatial_slice_ids,\
    slice_data_intervals=spatial_slice_intervals,deathvalue=deathvalue,bins=bins2d,\
    fontsize=fontsize,colormap=colormap,colorbarmaxval=colorbarmaxval)
  else:
    print("spatial slice lists not compatible or empty: skip spatial slices histograms")
  if(len(velocity_slice_ids)!=0 and len(velocity_slice_intervals)==2*len(velocity_slice_ids)):
    computes_phase_space_2d_sliced_histograms(p_groups,'v',slice_data_ids=spatial_slice_ids,\
    slice_data_intervals=spatial_slice_intervals,deathvalue=deathvalue,bins=bins2d,\
    fontsize=fontsize,colormap=colormap,colorbarmaxval=colorbarmaxval)
  else:
    print("velocity slice lists not compatible or empty: skip velocity slices histograms") 

# argument parser 
def generate_argument_parser():
  from argparse import ArgumentParser
  parser = ArgumentParser(description='read and plot a JOREK particle restart')
  parser.add_argument('--filename','-f',type=str,action='store',required=False,\
  dest='filename',default='part_restart.h5',\
  help='name of the particle restart file, default: part_restart.h5')
  parser.add_argument('--filepath','-fpath',type=str,action='store',required=False,\
  dest='filepath',default='.',help='path of the file to be read, default: .')
  parser.add_argument('--separator','-sep',type=str,action='store',required=False,\
  dest='separator',default='/',help='file separator, default: /')
  parser.add_argument('--bins1d','-pbins1d',type=int,nargs='*',action='store',\
  required=False,default=[1000,1000,1000],dest='bins1d',\
  help='number of or method for computing the bins for 1D histograms, default: 100')
  parser.add_argument('--bins2d','-pbins2d',type=int,nargs='*',action='store',\
  required=False,default=[1000,1000,1000],dest='bins2d',\
  help='number of or method for computing the bins for 2D histograms, default: 100')
  parser.add_argument('--fontsize','-fsize',type=int,action='store',required=False,\
  dest='fontsize',default=18,help='plot font size, default: 18')
  parser.add_argument('--colormap','-cmap',type=str,action='store',required=False,\
  dest='colormap',default='inferno',help='colormap of 2D histograms, default: inferno')
  parser.add_argument('--colorbarmaxval','-cbmv',type=float,nargs='*',action='store',\
  required=False,dest='colorbarmaxval',default=[1e0,1e0,1e0,1e0],\
  help='colorbar scaling values of 2D histograms, default: [1.,1.,1.,1.]')
  parser.add_argument('--deathvalue','-dval',type=int,action='store',required=False,\
  dest='deathvalue',default=0,\
  help='flag value below which a particle is set to inactive, default: 0') 
  parser.add_argument('--spatialsliceids','-slids',type=int,nargs='*',action='store',\
  required=False,default=[],dest='spatial_slice_ids',\
  help='indexes of the spatial data slices for computing sliced spatial 2D histograms, default: []')
  parser.add_argument('--velocitysliceids','-vlids',type=int,nargs='*',action='store',\
  required=False,default=[],dest='velocity_slice_ids',\
  help='indexes of the velocity data slices for computing sliced velocity 2D histograms, default: []')
  parser.add_argument('--spatialsliceintervals','-slint',type=float,nargs='*',action='store',\
  required=False,default=[],dest='spatial_slice_intervals',\
  help='intervals of the spatial data slices for computing sliced spatial 2D histograms, default: []')
  parser.add_argument('--velocitysliceintervals','-vlint',type=float,nargs='*',action='store',\
  required=False,default=[],dest='velocity_slice_intervals',\
  help='intervals of the velocity data slices for computing sliced velocity 2D histograms, default: []')
  return parser.parse_args()

# Run main -------------------------------------------- #
if __name__ == "__main__":
  args = generate_argument_parser()
  read_analyse_plot_jorek_restart(args.filename,args.filepath,args.separator,\
  deathvalue=args.deathvalue,bins1d=args.bins1d,bins2d=args.bins2d,\
  fontsize=args.fontsize,colormap=args.colormap,\
  spatial_slice_ids=args.spatial_slice_ids,\
  spatial_slice_intervals=args.spatial_slice_intervals,\
  velocity_slice_ids=args.velocity_slice_ids,\
  velocity_slice_intervals=args.velocity_slice_intervals,\
  colorbarmaxval=args.colorbarmaxval)

# ----------------------------------------------------- #
