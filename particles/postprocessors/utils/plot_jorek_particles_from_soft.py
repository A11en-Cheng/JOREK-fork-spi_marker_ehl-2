## Small tool used for plotting JOREK particles initiliased
## from SOFT orbit output

# read hdf5
def read_jorek_particle_data(filename,filepath,datasets,separator):
  import h5py
  import numpy as np
  data = {}
  fhandler = h5py.File("".join([filepath,separator,filename]))
  for setname in datasets:
    data[setname] = np.array(fhandler[setname])
  fhandler.close()
  return data

# compute particle RZ positions
def compute_particle_RZ_positions(x_pos):
  import numpy as np
  RZ_pos = np.zeros((x_pos.shape[0],2))
  RZ_pos[:,0] = np.sqrt(np.power(x_pos[:,0],2)+np.power(x_pos[:,1],2))
  RZ_pos[:,1] = x_pos[:,2]
  return RZ_pos

# compute the particle toroidal angle
def compute_particle_phi_positions(x_pos):
 from numpy import arctan2,pi
 phi = arctan2(x_pos[:,1],x_pos[:,0])
 phi[phi<0.] = 2.*pi+phi[phi<0.]
 return phi

# plot scatter 2d
def plot_scatter_2d(ax,xy,title='',xlab='x',ylab='y',fontsize=12,\
markertype='.',markersize='1',markercolor='b',aspect='equal'):
  mask = ((xy[:,0]==0.)&(xy[:,1]==0.))==False
  ax.scatter(xy[mask,0],xy[mask,1],s=markersize,\
  marker=markertype,c=markercolor)
  if(len(aspect)!=0):
    ax.set_aspect(aspect)
  ax.set_title(title,fontsize=fontsize)
  ax.set_xlabel(xlab,fontsize=fontsize)
  ax.set_ylabel(ylab,fontsize=fontsize)
  ax.tick_params(axis='x',labelsize=fontsize)
  ax.tick_params(axis='y',labelsize=fontsize)
  ax.grid()

# plot scatter 3d
def plot_scatter_3d(ax,xyz,title='',xlab='x',ylab='y',zlab='z',\
fontsize=12,markertype='.',markersize='1',markercolor='b',\
aspect='equal',rotate=False):
  mask = ((xyz[:,0]==0.)&(xyz[:,1]==0.)&(xyz[:,2]==0.))==False
  ax.scatter3D(xyz[mask,0],xyz[mask,1],xyz[mask,2],\
  s=markersize,marker=markertype,c=markercolor)
  if(len(aspect)!=0):
    ax.set_aspect(aspect)
  ax.set_title(title,fontsize=fontsize)
  ax.set_xlabel(xlab,fontsize=fontsize)
  ax.set_ylabel(ylab,fontsize=fontsize)
  ax.set_zlabel(zlab,fontsize=fontsize)
  ax.tick_params(axis='x',labelsize=fontsize)
  ax.tick_params(axis='y',labelsize=fontsize)
  ax.tick_params(axis='z',labelsize=fontsize)
  ax.xaxis.set_rotate_label(rotate)
  ax.yaxis.set_rotate_label(rotate)
  ax.zaxis.set_rotate_label(rotate)

# plot particle positions
def plot_particle_positions(x_pos,markersize=1,fontsize=16):
  from mpl_toolkits.mplot3d import Axes3D
  from matplotlib import pyplot as plt
  # compute the position in R,Z coordinates
  RZ_pos = compute_particle_RZ_positions(x_pos)
  # plot the figure
  fig = plt.figure(facecolor='white',edgecolor='white')
  axs = [fig.add_subplot(131),fig.add_subplot(132),fig.add_subplot(133,projection='3d')]
  plot_scatter_2d(axs[0],RZ_pos,title='JOREK particle RZ positions from SOFT',\
  xlab='R [m]',ylab='Z [m]',fontsize=fontsize,markertype='.',\
  markersize=markersize,markercolor='b',aspect='equal')
  plot_scatter_2d(axs[1],x_pos[:,0:2],title='JOREK particle xy positions from SOFT',\
  xlab='x [m]',ylab='y [m]',fontsize=fontsize,markertype='.',\
  markersize=markersize,markercolor='b',aspect='equal')
  plot_scatter_3d(axs[2],x_pos,title='JOREK particle positions from SOFT',\
  xlab='x [m]',ylab='y [m]',zlab='z [m]',fontsize=fontsize,markertype='.',\
  markersize=markersize,markercolor='b',aspect='equal',rotate=False)
  plt.show()

# plot particle poincare plot
def plot_poincare_plot(x_pos,phi_target=9e1,dphi=1.,fontsize=16,markersize=1):
  from numpy import maximum,minimum,pi
  from matplotlib import pyplot as plt
  # compute poincare plot
  dphi = float(dphi)
  phi_target = float(phi_target)
  RZ_pos = compute_particle_RZ_positions(x_pos)
  phi = compute_particle_phi_positions(x_pos)
  phi_target = phi_target*pi/1.8e2
  dphi = dphi*pi/1.8e2
  mask = ((phi >= maximum(phi_target-dphi,0.))&(phi <= minimum(phi_target+dphi,2.*pi)))
  RZ_poincare = RZ_pos[mask,:]
  # plot the figure
  fig = plt.figure(facecolor='white',edgecolor='white')
  ax = fig.add_subplot(111)
  plot_scatter_2d(ax,RZ_poincare,title="".join([\
  'JOREK particle RZ poincare from SOFT, phi: ',str(phi_target),' dphi: ',str(dphi)]),\
  xlab='R [m]',ylab='Z [m]',fontsize=fontsize,markertype='.',\
  markersize=markersize,markercolor='b',aspect='equal')
  plt.show()

# plot particle momentum
def plot_particle_momentum(ppar,pperp,fontsize=16,markersize=1):
  from numpy import array,zeros,float64
  from matplotlib import pyplot as plt
  momentum = zeros((ppar.shape[0],2),dtype=float64)
  momentum[:,0] = ppar
  momentum[:,1] = pperp
  fig = plt.figure(facecolor='white',edgecolor='white')
  ax = fig.add_subplot(111)
  plot_scatter_2d(ax,momentum,title='JOREK particle momentum from SOFT',\
  xlab='ppar/mc',ylab='pperp/mc',fontsize=fontsize,markertype='.',\
  markersize=markersize,markercolor='b',aspect='')
  plt.show() 

# function main
def read_and_plot_jorek_particles(filename='',filepath='.',\
datasets=['x','ppar','pperp'],separator= '/',phi_poincare=9e1,\
dphi_poincare=1.,fontsize=16,markersize=1):
  # read data from file
  particle_data = read_jorek_particle_data(filename,filepath,datasets,separator)
  plot_particle_momentum(particle_data['ppar'],particle_data['pperp'],\
  fontsize=fontsize,markersize=markersize)
  plot_poincare_plot(particle_data['x'],phi_target=phi_poincare,dphi=dphi_poincare,\
  fontsize=fontsize,markersize=markersize)
  plot_particle_positions(particle_data['x'],fontsize=fontsize,markersize=markersize)

# input parser
def generate_argument_parser():
  import argparse
  parser = argparse.ArgumentParser(description='file to print')
  parser.add_argument('filename',type=str,action='store',\
  help='name of the hdf5 to be read, required')
  parser.add_argument('--filepath','-fp',type=str,required=False,\
  dest='filepath',action='store',default='.',help='path to the hdf5 to be loaded, default: .')
  parser.add_argument('--separator','-sep',type=str,required=False,\
  dest='separator',action='store',default='/',help='file separator, default: /')
  parser.add_argument('--fontsize','-fsize',type=int,required=False,\
  dest='fontsize',action='store',default=16,help='fontsize for videos, default: 16')
  parser.add_argument('--markersize','-msize',type=int,required=False,\
  dest='markersize',action='store',default=1,help='marker size, default: 1')
  parser.add_argument('--datasets','-dset',type=str,nargs='*',\
  required=False,dest='datasets',action='store',default=['x','ppar','pperp'],\
  help='list contaning the name of the datasets to be loaded, default: [x,ppar,pperp]')
  parser.add_argument('--phi_poincare','-phi',type=float,\
  required=False,dest='phi_poincare',action='store',default=9e1,\
  help='toroidal angle for poincare plot in [deg], default=1.5')
  parser.add_argument('--dphi_poincare','-dphi',type=float,\
  required=False,dest='dphi_poincare',action='store',default=1.,\
  help='toroidal angle interval for poincare plot in [deg], default: 1')
  return parser.parse_args()

if __name__ == "__main__":
  args = generate_argument_parser()
  read_and_plot_jorek_particles(args.filename,filepath=args.filepath,separator=args.separator,\
  datasets=args.datasets,phi_poincare=args.phi_poincare,dphi_poincare=args.dphi_poincare,
  markersize=args.markersize,fontsize=args.fontsize)
