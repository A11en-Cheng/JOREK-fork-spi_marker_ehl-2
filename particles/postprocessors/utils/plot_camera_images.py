# --------------------------------------------------------------- #
# Read and plot the camera pixel and filter intensities HDF5 file
# --------------------------------------------------------------- #
# Method used for adapting the sizes of position coordinates
def adapt_positions_for_plot(positions):
  from numpy import newaxis,transpose
  if(positions.size==0):
    return positions
  if(len(positions.shape)<2):
    positions = positions[newaxis,:]
  positions = transpose(positions,axes=[1,0])
  return positions

# Method used for adapitng images for plot
def adapt_images_for_plot(images):
  from numpy import newaxis,transpose
  if(images.size==0):
    return images
  if(len(images.shape)<3):  
    images  = images[:,:,newaxis]
  if(len(images.shape)<4):
    images = images[:,:,:,newaxis]
  images = transpose(images,axes=[3,2,1,0])
  return images

# Method used for plotting 4d arrays as 2d images
# inputs:
#   frames_spectra: (nx,ny,n_spectra,n_times) array to plot
#   x_positions:    (nx,n_times)(optional) x positions for scaling
#                   default: no scaling is used
#   y_positions:    (ny,n_times)(optional) y positions for scaling
#                   default: no scaling is used
#   title:          (string) figure title, default: empty
def imshow_4d(frames_spectra,x_positions=[],y_positions=[],title=""):
  from numpy import newaxis
  from matplotlib.pyplot import figure,imshow,colorbar
  from matplotlib.pyplot import title as tit
  # plot with extension
  for frame_id,frame in enumerate(frames_spectra):
    if(not ((len(x_positions)==0) and (len(y_positions)==0))):
      dx = 0.5*(x_positions[1,frame_id]-x_positions[0,frame_id])
      dy = 0.5*(y_positions[1,frame_id]-y_positions[0,frame_id])
      ext = [x_positions[0,frame_id]-dx, x_positions[-1,frame_id]+dx,\
      y_positions[0,frame_id]-dy,y_positions[-1,frame_id]+dy]
    else:
      ext = None
    for spectrum_id,spectrum in enumerate(frame):
      figure()
      imshow(spectrum,extent=ext)
      colorbar()
      tit("".join([title,' spectrum N# ',str(spectrum_id+1),\
      ' time N#: ',str(frame_id+1)]))   
  return  

# load and plot a jorek image set
# inputs:
#   filenames:                    (list) list of names of the hdf5 files containing the images
#   filepath:                     (string) path to the file containing the images
#   image_datasetname:            (string) name of the image hdf5 dataset
#   pixel_x_position_datasetname: (string) name of the x coordinate map hdf5 dataset
#   pixel_y_position_datasetname: (string) name of the y coordinate map hdf5 dataset
#   separator:                    (string) filename-filepath separator
# outputs:
#   images: (n_images,ntimes,nspectrum,ny,nx) set of jorek images
def load_and_plot_jorek_images(filenames=[],filepath=".",image_datasetname="",\
pixel_x_postion_datasetname="",pixel_y_postion_datasetname="",\
separator="/"):
  import h5py
  from numpy import array,isnan,transpose
  from matplotlib.pyplot import show
  # Initialisation
  pixel_intensities = []
  # Read image
  for filename in filenames:
    fhandler = h5py.File("".join([filepath,separator,filename]),'r')
    pixel_filter_intensities = array(fhandler[image_datasetname])
    x_positions = array(fhandler[pixel_x_postion_datasetname])
    y_positions = array(fhandler[pixel_y_postion_datasetname])
    fhandler.close()
    # Reorder pixel and filter intensities for plotting
    x_positions = transpose(x_positions,axes=[1,0])
    y_positions = transpose(y_positions,axes=[1,0])
    pixel_filter_intensities = transpose(pixel_filter_intensities,axes=[3,4,0,1,2])
    pixel_intensities.append(pixel_filter_intensities[0])
    filter_intensities = pixel_filter_intensities[1]
    # Plot pixel and filter intensities
    imshow_4d(pixel_intensities[-1],x_positions=x_positions,y_positions=y_positions,title="Image for" )
    imshow_4d(filter_intensities,x_positions=x_positions,y_positions=y_positions,title="Filter for" )
    if((isnan(pixel_intensities[-1])).any()):
      print('Warning: found pixel(s) with NaN intensity')
  show()
  return pixel_intensities

# load and plot a generic image set
# inputs:
#   filename:                     (string) name of the hdf5 containing the images
#   filepath:                     (string) path to the file containing the images
#   image_datasetname:            (string) name of the image hdf5 dataset
#   pixel_x_position_datasetname: (string) name of the x coordinate map hdf5 dataset
#   pixel_y_position_datasetname: (string) name of the y coordinate map hdf5 dataset
#   separator:                    (string) filename-filepath separator
# outputs:
#   image: (ntimes,nspectrum,ny,nx) set of shape adapted generic images
def load_and_plot_generic_images(filename="",filepath=".",image_datasetname="",\
pixel_x_postion_datasetname="",pixel_y_postion_datasetname="",\
separator="/"):
  import h5py
  from numpy import array,isnan
  from matplotlib.pyplot import show
  # Initialisation
  x_positions = array([]); y_positions = array([]);
  # Read image
  fhandler = h5py.File("".join([filepath,separator,filename]),'r')
  image = array(fhandler[image_datasetname])
  if(pixel_x_postion_datasetname in fhandler):
    x_positions = array(fhandler[pixel_x_postion_datasetname])
  if(pixel_y_postion_datasetname in fhandler):
    y_positions = array(fhandler[pixel_y_postion_datasetname])
  fhandler.close()
  # Reorder pixel and filter intensities for plotting
  x_positions = adapt_positions_for_plot(x_positions)
  y_positions = adapt_positions_for_plot(y_positions)
  image = adapt_images_for_plot(image)
  # Plot pixel and filter intensities
  imshow_4d(image,x_positions=x_positions,y_positions=y_positions,title="Comparison image")
  if((isnan(image)).any()):
    print('Warning: found generic image with NaN intensity')
  show()
  return image

# compute the differences between the jorek images and the test images.
# the number of spectra and the number of pixels of the two image set must be the same
# inputs:
#   jorek_image: (n_images,ntimes,nspectra,ny,nx) set of jorek images per time per spectrum
#   test_image:  (ntimes,nspectra,ny,nx) set of test images per time per spectrum
def compute_and_plot_image_differences(jorek_images,test_image):
  from numpy import array,zeros,float64,amax,abs
  from matplotlib.pyplot import show
  # loop on the JOREK images
  for jorek_image_id,jorek_image in enumerate(jorek_images):
    # check image shape compatibility
    if(jorek_image.shape[1:]!=test_image.shape[1:]):
      raise Exception("".join(["JOREK and test image sets have different shapes jorek shape: ",\
      str(jorek_image.shape)," test image shape: ",str(test_image.shape)]))
    # initialisation
    image_error = zeros((jorek_image.shape[0]*test_image.shape[0],\
    jorek_image.shape[1],jorek_image.shape[2],jorek_image.shape[3],),dtype=float64)
    normalised_image_error = zeros(image_error.shape,dtype=float64)
    # compute differences
    for test_id,test_frame in enumerate(test_image):
      for jorek_id,jorek_frame in enumerate(jorek_image):
        for spectra_id,spectra in enumerate(jorek_frame):
          image_error[test_id*jorek_image.shape[0],spectra_id,:,:] = spectra-test_frame[spectra_id]
          normalised_image_error[test_id*jorek_image.shape[0],spectra_id,:,:] = \
          (spectra/amax(abs(spectra))) - (test_frame[spectra_id]/amax(abs(test_frame[spectra_id])))
    # plot images
    imshow_4d(image_error,title="".join(["Differences between the JOREK image N# ",str(jorek_image_id+1),\
    " and the test image"]))
    imshow_4d(normalised_image_error,title="".join(["Normalised differences between the JOREK image N# ",\
    str(jorek_image_id+1)," and the test image"]))
    imshow_4d(abs(image_error),title="".join(["Error between the JOREK image N# ",str(jorek_image_id+1),\
    " and the test image"]))
    imshow_4d(abs(normalised_image_error),title="".join(["Normalised error between the JOREK image N# ",\
    str(jorek_image_id+1)," and test the image"]))
  show()

def generate_argument_parser():
  import argparse
  parser = argparse.ArgumentParser(\
  description='read and plot fast camera image and filter')
  parser.add_argument('-f','--filenames',type=str,nargs='*',\
  action='store',dest='filenames',default=['pixel_filter_intensities.h5'],\
  help='name of the file to be read')
  parser.add_argument('--filepath','-fp',type=str,required=False,\
  dest='filepath',action='store',default='.',\
  help='path to the image/filter hdf5 files to be loaded, default: .')
  parser.add_argument('--separator','-sep',type=str,required=False,\
  dest='separator',action='store',default='/',help='file separator, default: /')
  parser.add_argument('--image_datasetname','-imdset',type=str,required=False,\
  dest='image_datasetname',action='store',default='pixel_filter_intensities',\
  help='image-filter dataset name')
  parser.add_argument('--x_pixel_coordinates_datasetname','-icx',type=str,\
  required=False,dest='pixel_x_postion_datasetname',action='store',\
  default='x_pixel_coordinates',help='x position coordinates dataset name')
  parser.add_argument('--y_pixel_coordinates_datasetname','-icy',type=str,\
  required=False,dest='pixel_y_postion_datasetname',action='store',\
  default='y_pixel_coordinates',help='y position coordinates dataset name')
  parser.add_argument('-fgen','--genericfilename',type=str,action='store',required=False,\
  dest='generic_filename',default='',help='name of the SOFT code filename to be read')
  parser.add_argument('--generic_image_datasetname','-gimdset',type=str,required=False,\
  dest='generic_image_datasetname',action='store',default='image',\
  help='SOFT image-filter dataset name')
  return parser.parse_args()

# Run main ------------------------------------------------------ #
if __name__ == "__main__":
  args = generate_argument_parser()
  if(len(args.filenames)!=0):
    jorek_images=load_and_plot_jorek_images(filenames=args.filenames,\
    filepath=args.filepath,image_datasetname=args.image_datasetname,\
    pixel_x_postion_datasetname=args.pixel_x_postion_datasetname,\
    pixel_y_postion_datasetname=args.pixel_y_postion_datasetname,\
    separator=args.separator)
  if(len(args.generic_filename)!=0):
    generic_image=load_and_plot_generic_images(filename=args.generic_filename,\
    filepath=args.filepath,image_datasetname=args.generic_image_datasetname,\
    separator=args.separator)
  if((len(args.generic_filename)!=0) and (len(args.filenames)!=0)):
    compute_and_plot_image_differences(jorek_images,generic_image)

# --------------------------------------------------------------- #
