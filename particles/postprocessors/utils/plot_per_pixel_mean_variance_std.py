# ----------------------------------------------------- #
# Program for reading and plots contributing light      #
# light sources per pixel average, variance and std     #
# from HDF5_array2D_saving                              #
# ----------------------------------------------------- #
# Program functions ----------------------------------- #
# Read datasets from HDF5 files
def comput_major_radius(x,y):
  from numpy import power,sqrt
  return sqrt(power(x,2)+power(y,2))

def read_datasets_from_hdf5(filenames,filepath,datasets,separator):
  from h5py import File
  from numpy import array
  data_list = []
  for filename in filenames:
    data = dict()
    fhandler = File("".join([filepath,separator,filename]),'r')
    for setname in datasets:
      data[setname] = array(fhandler[setname])
    data_list.append(data)  
    fhandler.close()
  return data_list

# write statistics in hdf5 file
def write_statistics_in_hdf5(average,variance,stddev,filename,filepath,separator):
  from h5py import File
  from numpy import float64 as npfloat64
  fhandler = File("".join([filepath,separator,filename]),'w')
  dset = fhandler.create_dataset('stat_size',len(average.shape),dtype=int,data=average.shape)
  dset = fhandler.create_dataset('average',average.shape,dtype=average.dtype,data=average)
  dset = fhandler.create_dataset('variance',variance.shape,dtype=variance.dtype,data=variance)
  dset = fhandler.create_dataset('standard_deviation',stddev.shape,dtype=stddev.dtype,data=stddev)
  fhandler.close()

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
  num_spectra = light_data[0]["contributing_light_intensities"].shape[0]
  test_spectra = [dataset["contributing_light_intensities"].shape[0]==num_spectra for dataset in light_data]
  if(not all(test_spectra)):
    exit("".join(['Number of spectra must be the same for all light datasets ref N# spectra: ',\
    str(num_spectra),' test N# spectra: ',str(test_spectra)]))

# compute the average light intensity for each pixel
def compute_per_pixel_average(n_pixels_spectra,light_data):
  from numpy import zeros
  from numpy import float64 as npfloat64
  n_particles = 0
  average_intensity = zeros((n_pixels_spectra[1],n_pixels_spectra[2],\
  n_pixels_spectra[0]),dtype=npfloat64)
  for light_dataset in light_data:
    pixel_ids = light_dataset["contributing_light_pixel_indexes"]
    light_intensities = light_dataset["contributing_light_intensities"]
    for timesId,times in enumerate(pixel_ids):
      for pixelId,pixel in enumerate(times):
        average_intensity[pixel[0]-1,pixel[1]-1,:] = average_intensity[pixel[0]-1,pixel[1]-1,:] + \
        light_intensities[timesId,pixelId,:]
    n_particles = n_particles + light_intensities.shape[1]
  return average_intensity/n_particles
    
# compute the light variance and standard deviation
def compute_per_pixel_variance_stddev(n_pixels_spectra,average,light_data):
  from numpy import zeros,power,sqrt
  from numpy import float64 as npfloat64
  n_particles = 0
  variance_intensity = zeros((n_pixels_spectra[1],n_pixels_spectra[2],\
  n_pixels_spectra[0]),dtype=npfloat64)
  for light_dataset in light_data:
    pixel_ids = light_dataset["contributing_light_pixel_indexes"]
    light_intensities = light_dataset["contributing_light_intensities"]
    for timesId,times in enumerate(pixel_ids):
      for pixelId,pixel in enumerate(times):
        variance_intensity[pixel[0]-1,pixel[1]-1,:] = variance_intensity[pixel[0]-1,pixel[1]-1,:] + \
        power((light_intensities[timesId,pixelId,:]-average[pixel[0]-1,pixel[1]-1,:]),2)
    n_particles = n_particles + light_intensities.shape[1]
  variance_intensity = variance_intensity/(n_particles-1)
  return variance_intensity,sqrt(variance_intensity) 

# Method used for plotting 4d arrays as 2d images
# inputs:
#   frames_spectra: (nx,ny,n_spectra) array to plot
#   x_positions:    (nx)(optional) x positions for scaling
#                   default: no scaling is used
#   y_positions:    (ny)(optional) y positions for scaling
#                   default: no scaling is used
#   title:          (string) figure title, default: empty
def imshow_3d(frame,x_positions=[],y_positions=[],title=""):
  from numpy import newaxis
  from matplotlib.pyplot import figure,imshow,colorbar,show
  from matplotlib.pyplot import title as tit
  # plot with extension
  if(not ((len(x_positions)==0) and (len(y_positions)==0))):
    dx = 0.5*(x_positions[1]-x_positions[0])
    dy = 0.5*(y_positions[1]-y_positions[0])
    ext = [x_positions[0]-dx, x_positions[-1]+dx,\
    y_positions[0]-dy,y_positions[-1]+dy]
  else:
    ext = None
  for spectrum_id,spectrum in enumerate(frame):
    figure()
    imshow(spectrum,extent=ext)
    colorbar()
    tit("".join([title,' spectrum N# ',str(spectrum_id+1)]))
  show()
  return

# Main function
# The light source datasets are:
#   contributing_light_positions: position of the contributing
#     lights per lens point
#   contributing_light_intensities: spectral intensity of the
#     contributing lights per lens point
#   limiter_major_radius
#     tokamak limiter (first wall) major radius
#   limiter_vertical_coordinate:
#     tokamak limiter (first wall) vertical position
#   contributing_light_pixel_indexes:
#     pixel indexes of each light source
#   contributing_light_local_pixel_coordinates:
#     local pixel coordinates of each light source
# The camera datasets are:
#   point_on_lens_positions: positions of the points on 
#     the camera lens
#   image_plane_vertices: vertices of the image planes
#   image_plane_directions: viewing directions of each
#     image plane
#   n_pixels_n_spectra: n_spectra,n_x_pixel,n_y_pixel
def load_and_plot_contributing_light_sources(\
light_filenames=[],camera_filenames=[],filepath=".",\
light_dataset=["contributing_light_positions","contributing_light_intensities",\
"limiter_major_radius","limiter_vertical_coordinate",\
"contributing_light_pixel_indexes","contributing_light_local_pixel_coordinates"],\
camera_dataset=["point_on_lens_positions","image_plane_vertices",\
"image_plane_directions","n_pixels_n_spectra"],stats_filename="image_statistics.h5",\
separator="/",font_size=16):
  from numpy import transpose
  # Read data from files
  light_data = read_datasets_from_hdf5(light_filenames,filepath,\
  light_dataset,separator)  
  camera_data = read_datasets_from_hdf5(camera_filenames,filepath,\
  camera_dataset,separator)
  # Check consistency of the data
  check_equal_light_data(light_data)
  # compute statistics on the image, generate the pixel mesh assuming constant plane
  print("computing light intensity average ...")
  average = compute_per_pixel_average(camera_data[-1]["n_pixels_n_spectra"],light_data)
  print("computing light intensity average completed!")
  print("computing light intensity variance and standard deviation ...")
  variance,stddev = compute_per_pixel_variance_stddev(camera_data[-1]["n_pixels_n_spectra"],\
  average,light_data)
  print("computing light intensity variance and standard deviation completed!")
  print("write data in hdf5 file")
  write_statistics_in_hdf5(average,variance,stddev,stats_filename,filepath,separator)
  print("Plotting data")
  imshow_3d(transpose(average,axes=[2,1,0]),title="Per pixel light intensity average for ")
  imshow_3d(transpose(variance,axes=[2,1,0]),title="Per pixel light intensity variance for ")
  imshow_3d(transpose(stddev,axes=[2,1,0]),title="Per pixel light intensity std deviation for ")

def generate_argument_parser():
  import argparse
  parser = argparse.ArgumentParser(\
  description='plot the mean value, the variance and the standard deviation of the light intensity of each pixel')
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
  "limiter_vertical_coordinate","contributing_light_pixel_indexes",\
  "contributing_light_local_pixel_coordinates"],help='light datasets to be loaded')
  parser.add_argument('--camera_dataset','-cdset',type=str,nargs='*',required=False,\
  dest='camera_dataset',action='store',default=["point_on_lens_positions","image_plane_vertices",\
  "image_plane_directions","n_pixels_n_spectra"],help='camera datasets to be loaded')
  parser.add_argument('--fontsize','-font',type=int,required=False,\
  dest='fontsize',action='store',default=16,help='plot font size, default: 16')
  parser.add_argument('--stats_filename','-sf',type=str,dest='stats_filename',\
  default = "image_statistics.h5",action='store',required=False,\
  help='filename in which the image statistics are written, default: image_statitics.h5')
  return parser.parse_args()

# Run main -------------------------------------------- #
if __name__ == "__main__":
  args = generate_argument_parser()
  fig=load_and_plot_contributing_light_sources(light_filenames=args.light_names,\
  camera_filenames=args.camera_names,filepath=args.filepath,\
  stats_filename=args.stats_filename,separator=args.separator,\
  light_dataset=args.light_dataset,camera_dataset=args.camera_dataset,\
  font_size=args.fontsize)
# ----------------------------------------------------- #

