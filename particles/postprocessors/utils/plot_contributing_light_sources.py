# ----------------------------------------------------- #
# Program for reading and plots contributing light      #
# light sources to each point on a lens toghether with  #
# the image plane and their viewing directions          #
# ----------------------------------------------------- #
# Program functions ----------------------------------- #
# Read datasets from HDF5 files
def comput_major_radius(x,y):
  from numpy import power,sqrt
  return sqrt(power(x,2)+power(y,2))

def read_datasets_from_hdf5(filenames,filepath,datasets,separator):
  import h5py
  import numpy as np
  data_list = []
  for filename in filenames:
    data = []
    fhandler = h5py.File("".join([filepath,separator,filename]))
    for setname in datasets:
      data.append(np.array(fhandler[setname]))
    data_list.append(data)  
    fhandler.close()
  return data_list

# Check the consistency of the light data:
# the number of spectra is the same for
# all datasets. Structure of the light data 
# for each light data:
# 0: position of the contributing light sources
# 1: spectral intensities of the contributing light sources
# 2: tokamak limiter (first wall) major radius
# 3: tokamak limiter (first wall) vertical position
def check_equal_light_data(light_data):
  from sys import exit
  num_spectra = light_data[0][1].shape[0]
  test_spectra = [dataset[1].shape[0]==num_spectra for dataset in light_data]
  if(not all(test_spectra)):
    exit("".join(['Number of spectra must be the same for all light datasets ref N# spectra: ',\
    str(num_spectra),' test N# spectra: ',str(test_spectra)]))

# Plot light and camera data
# structure of the light data for each light data
# 0: position of the contributing light sources
# 1: spectral intensities of the contributing light sources
# 2: tokamak limiter (first wall) major radius
# 3: tokamak limiter (first wall) vertical position
# structure of the camera data for each camera data:
# 0: points on the lens
# 1: vertices of the image planes
# 2: view directions of the image planes
def plot_light_and_camera_data(light_data,camera_data,n_tor=100,markersize=1,\
linewidth=3,fontsize=16,colormap_scaling=1e0):
  import numpy as np
  from mpl_toolkits.mplot3d import Axes3D
  from matplotlib import pyplot as plt
  # extract the number of spectra
  n_spectra = light_data[0][1].shape[0]
  # extract limiter coordinates considered equals for all dataset
  R_limiter = light_data[0][2]
  Z_limiter = light_data[0][3]
  cosphi = np.cos(np.linspace(0,2*np.pi,num=n_tor,dtype=np.float64))
  sinphi = np.sin(np.linspace(0,2*np.pi,num=n_tor,dtype=np.float64))
  R_max_limiter = np.amax(R_limiter)
  R_min_limiter = np.amin(R_limiter)
  x_limiter_min = R_min_limiter*cosphi
  x_limiter_max = R_max_limiter*cosphi
  y_limiter_min = R_min_limiter*sinphi
  y_limiter_max = R_max_limiter*sinphi
  del cosphi,sinphi,R_max_limiter,R_min_limiter
  # find the maximum of the specta
  max_spectrum = 0
  for dataset in light_data:
    max_spectrum = colormap_scaling*max(max_spectrum,np.amax(dataset[1]))
  # do scatter plots
  for spectra_id in range(n_spectra):
    fig = plt.figure(facecolor='white',edgecolor='white')
    axs = [fig.add_subplot(131),fig.add_subplot(132),fig.add_subplot(133,projection='3d')]
    n_particles = 0
    for dataset in light_data:
      positions = dataset[0]
      spectra = dataset[1]
      n_particles = n_particles + positions.shape[1]
      for time_id,spectrum in enumerate(np.transpose(spectra[spectra_id,:,:],[1,0])):
        mask = np.where(spectrum>0)
        # aggregated top view
        axs[0].scatter(positions[time_id,mask,0],positions[time_id,mask,1],marker='.',\
        s=markersize,c=spectrum[mask],cmap='inferno',\
        vmin=0,vmax=max_spectrum)
        axs[0].set_aspect('equal')
        axs[0].set_facecolor([0,0,0])
        axs[0].set_title('Aggregated top view',fontsize=fontsize,color='red')
        axs[0].set_xlabel('x [m]',fontsize=fontsize,color='red')
        axs[0].set_ylabel('y [m]',fontsize=fontsize,color='red')
        axs[0].tick_params(axis='x',labelsize=fontsize,colors='red')
        axs[0].tick_params(axis='y',labelsize=fontsize,colors='red')
        # aggregated front view
        axs[1].scatter(comput_major_radius(positions[time_id,mask,0],positions[time_id,mask,1]),\
        positions[time_id,mask,2],s=markersize,c=spectrum[mask],marker='.',cmap='inferno',\
        vmin=0,vmax=max_spectrum)
        axs[1].set_aspect('equal')
        axs[1].set_facecolor([0,0,0])
        axs[1].set_title('Aggregated frontal view',fontsize=fontsize,color='red')
        axs[1].set_xlabel('R [m]',fontsize=fontsize,color='red')
        axs[1].set_ylabel('Z [m]',fontsize=fontsize,color='red')
        axs[1].tick_params(axis='x',labelsize=fontsize,colors='red')
        axs[1].tick_params(axis='y',labelsize=fontsize,colors='red')
        # 3D view
        axs[2].scatter3D(positions[time_id,mask,0],positions[time_id,mask,1],positions[time_id,mask,2],\
        s=markersize,c=spectrum[mask],marker='.',cmap='inferno',vmin=0,vmax=max_spectrum)
        axs[2].set_facecolor([0,0,0])
        axs[2].set_aspect('equal')
        axs[2].set_title('3D view',fontsize=fontsize,color='red')
        axs[2].set_xlabel('x [m]',fontsize=fontsize,color='red')
        axs[2].set_ylabel('y [m]',fontsize=fontsize,color='red')
        axs[2].set_zlabel('z [m]',fontsize=fontsize,color='red')
        axs[2].tick_params(axis='x',labelsize=fontsize,colors='red')
        axs[2].tick_params(axis='y',labelsize=fontsize,colors='red')
        axs[2].tick_params(axis='z',labelsize=fontsize,colors='red')
        axs[2].xaxis.set_rotate_label(False)
        axs[2].yaxis.set_rotate_label(False)
        axs[2].zaxis.set_rotate_label(False)
  # Plotting limiter
  axs[0].plot(x_limiter_min,y_limiter_min,color='red',linewidth=linewidth)
  axs[0].plot(x_limiter_max,y_limiter_max,color='red',linewidth=linewidth)
  axs[1].plot(R_limiter,Z_limiter,color='red',linewidth=linewidth)
  del x_limiter_min,y_limiter_min,x_limiter_max,y_limiter_max,R_limiter,Z_limiter
  # Plotting the camera properties
  for dataset in camera_data:
    # Plot points on lens
    lens_point_position_avg = []
    for points in dataset[0]:
      lens_point_position_avg.append(np.array([np.mean(points[:,0]),np.mean(points[:,1]),np.mean(points[:,2])]))
      axs[0].scatter(points[:,0],points[:,1],marker='.',s=markersize,c='green')
      axs[2].scatter3D(points[:,0],points[:,1],points[:,2],marker='.',s=markersize,c='green')
    # Plot image planes
    for plane in dataset[1]:
      point4 = plane[1,:]-plane[0,:]+plane[2,:]
      axs[0].plot([plane[0,0],plane[1,0]],[plane[0,1],plane[1,1]],color='green',linewidth=linewidth)
      axs[0].plot([plane[0,0],plane[2,0]],[plane[0,1],plane[2,1]],color='green',linewidth=linewidth)
      axs[0].plot([plane[1,0],point4[0]],[plane[1,1],point4[1]],color='green',linewidth=linewidth)
      axs[0].plot([plane[2,0],point4[0]],[plane[2,1],point4[1]],color='green',linewidth=linewidth)
      axs[2].plot3D([plane[0,0],plane[1,0]],[plane[0,1],plane[1,1]],[plane[0,2],plane[1,2]],\
      color='green',linewidth=linewidth)
      axs[2].plot3D([plane[0,0],plane[2,0]],[plane[0,1],plane[2,1]],[plane[0,2],plane[2,2]],\
      color='green',linewidth=linewidth)
      axs[2].plot3D([plane[1,0],point4[0]],[plane[1,1],point4[1]],[plane[1,2],point4[2]],\
      color='green',linewidth=linewidth)
      axs[2].plot3D([plane[2,0],point4[0]],[plane[2,1],point4[1]],[plane[2,2],point4[2]],\
      color='green',linewidth=linewidth)
    del point4
    # Plot view directions
    for time_id,direction in enumerate(dataset[2]):
      point = lens_point_position_avg[time_id]
      axs[0].quiver(point[0],point[1],direction[0],direction[1],color='green',\
      linewidth=linewidth)
      axs[2].quiver(point[0],point[1],point[2],direction[0],direction[1],direction[2],\
      color='green',linewidth=linewidth)
    del lens_point_position_avg
  # generating image
  plt.suptitle("".join(['Point light source intensities, spectrum N#:',str(spectra_id)]),\
  fontsize=fontsize)
  plt.show() 

# Main function
# The light source datasets are:
#   contributing_light_positions: position of the contributing
#     lights per lens point
#   contributing_light_intensities: spectral intensity of the
#     contributing lights per lens point
#   limiter_major_radius
#     tokamak limiter (first wall) major radius
#   limiter_vertical_coordinate
#     tokamak limiter (first wall) vertical position
# The camera datasets are:
#   point_on_lens_positions: positions of the points on 
#     the camera lens
#   image_plane_vertices: vertices of the image planes
#   image_plane_directions: viewing directions of each
#     image plane
def load_and_plot_contributing_light_sources(\
light_filenames=[],camera_filenames=[],filepath=".",\
light_dataset=["contributing_light_positions","contributing_light_intensities",\
"limiter_major_radius","limiter_vertical_coordinate"],\
camera_dataset=["point_on_lens_positions","image_plane_vertices",\
"image_plane_directions"],separator="/",n_tor_mesh=100,marker_size=1,\
line_width=3,font_size=16,color_scaling=1e0):
  # Read data from files
  light_data = read_datasets_from_hdf5(light_filenames,filepath,\
  light_dataset,separator)  
  camera_data = read_datasets_from_hdf5(camera_filenames,filepath,\
  camera_dataset,separator)
  # Check consistency of the data
  check_equal_light_data(light_data)
  # Plot data
  plot_light_and_camera_data(light_data,camera_data,\
  n_tor=n_tor_mesh,markersize=marker_size,linewidth=line_width,\
  fontsize=font_size,colormap_scaling=color_scaling)

def generate_argument_parser():
  import argparse
  parser = argparse.ArgumentParser(\
  description='read and plot particles contributing to an image')
  parser.add_argument('--light_names','-ln',type=str,nargs='*',\
  action='store',required=True,help='list of ligth files to be read')
  parser.add_argument('--camera_names','-cn',type=str,nargs='*',\
  action='store',required=True,help='list of camera files to be read')
  parser.add_argument('--filepath','-fp',type=str,required=False,\
  dest='filepath',action='store',default='.',\
  help='path to light and camera hdf5 files to be loaded, default: .')
  parser.add_argument('--separator','-sep',type=str,required=False,\
  dest='separator',action='store',default='/',help='file separator, default: /')
  parser.add_argument('--light_dataset','-ldset',type=str,nargs='*',required=False,\
  dest='light_dataset',action='store',default=["contributing_light_positions",\
  "contributing_light_intensities","limiter_major_radius",\
  "limiter_vertical_coordinate"],help='light datasets to be loaded')
  parser.add_argument('--camera_dataset','-cdset',type=str,nargs='*',required=False,\
  dest='camera_dataset',action='store',default=["point_on_lens_positions","image_plane_vertices",\
  "image_plane_directions"],help='camera datasets to be loaded')
  parser.add_argument('--fontsize','-font',type=int,required=False,\
  dest='fontsize',action='store',default=16,help='plot font size, default: 16')
  parser.add_argument('--markersize','-msize',type=int,required=False,\
  dest='markersize',action='store',default=1,help='size of the markers, default: 1')
  parser.add_argument('--linewidth','-lwidth',type=int,required=False,\
  dest='linewidth',action='store',default=3,help='plot line width, default: 3')
  parser.add_argument('--colorscalefactor','-cfact',type=float,required=False,\
  dest='colorscalefactor',action='store',default=1e0,help='scaling factor of the colorbar, default: 1e0')
  parser.add_argument('--n_tor_mesh','-ntor',type=int,required=False,\
  dest='n_tor_mesh',action='store',default=100,help='number of poloidal points for plotting the wall, default: 100')
  return parser.parse_args()

# Run main -------------------------------------------- #
if __name__ == "__main__":
  args = generate_argument_parser()
  fig=load_and_plot_contributing_light_sources(light_filenames=args.light_names,\
  camera_filenames=args.camera_names,filepath=args.filepath,separator=args.separator,\
  light_dataset=args.light_dataset,camera_dataset=args.camera_dataset,
  n_tor_mesh=args.n_tor_mesh,marker_size=args.markersize,line_width=args.linewidth,\
  font_size=args.fontsize,color_scaling=args.colorscalefactor)
# ----------------------------------------------------- #

