# --------------------------------------------- #
# The python script generates a new fast camera #
# example file modifying the lines requested by #
# the used. The main application is to          #
# automatically generate example files with the #
# right input parameters for running non        #
# regression tests.                             #
# --------------------------------------------- #
# Tools -----------------------------------------
# Append path of packages
# inputs:
#   package_path: (list,string) package paths to add
#   position:     (integer) position in the path list  
def insert_packages_paths(package_paths,position=0):
  import os; import sys;
  for package_path in package_paths:
    sys.path.insert(position,os.path.abspath(package_path))

# Modify the particle example file an write it in a new file
# inputs:
#   example_file:          (string) name of the example file
#   example_path:          (string) path to the example directory
#   dest_path:             (string) destination directory of 
#                          the new example file
#   name_modifier:         (string) modifier of the example name
#   separator:             (string) path separator
#   dict_params_to_change: (dict) dictionary containing the
#                          parameter to change
def create_modified_example(example_file,example_path,\
  dest_path,name_modifier,separator,dict_params_to_change):
  from modify_fortran_file import ModifyFortranFile
  file_modifier = ModifyFortranFile(example_file,\
  source_dir=example_path,dest_dir=dest_path,\
  modify_name=name_modifier,separator=separator)
  file_modifier.modify_variables(dict_params_to_change)
  del file_modifier

# Argument parser -------------------------------
# outputs:
#   parser: (namespace) namespace containing the 
#           inputs as attributes
def generate_argument_parser():
  from argparse import ArgumentParser
  from dictionary_parser_class import ParseDictionary
  # define parser
  parser = ArgumentParser(description=\
  "Automatically generate Fast camera example files with user's inputs")
  parser.add_argument('--example-dir','-ed',type=str,required=False,\
  action='store',dest='example_dir',default='./particles/postprocessors/examples',\
  help='Directory path to the examples, default: ./particles/postprocessors/examples')
  parser.add_argument('--example-name','-efn',type=str,required=False,\
  action='store',dest='example_name',default='camera_RE_gyroaverage_synchrotron_example',\
  help='Filename of the example to be modified, default: camera_RE_gyroaverage_synchrotron_example')
  parser.add_argument('--destination-dir','-dd',type=str,required=False,\
  action='store',dest='dest_dir',default='.',\
  help='Destination directory of the modified example file, default: .')
  parser.add_argument('--name-modifier','-nm',type=str,required=False,\
  action='store',dest='name_modifier',default='.',\
  help='Characters recognising a modified example, default: _modified')
  parser.add_argument('--separator','-sep',type=str,required=False,\
  action='store',dest='separator',default='/',\
  help='Path separator, default: /')
  parser.add_argument('-dptc','--dict-params-to-change',dest='params_to_change_dict',\
  action=ParseDictionary,required=None,default={'n_frames':1,'n_times':1,\
  'fields_filename':'jorek_restart','image_filename':'result_intensities',\
  'n_groups':1,'n_spectra':1,'n_wavelengths':40,'n_int_camera_param':5,\
  'n_real_camera_param':9,'write_gc_in_txt':False,'particle_filenames':['particle_restart'],\
  'min_spectra':[3e-6],'max_spectra':[3.5e-6],'pinhole_position':[-8.86e-1,-4.002,-3.32e-1],\
  'int_camera_param':[0,600,600,0,1],'real_camera_param':[5.23e-1,5.23e-1,\
  1.5707963267948966,9.998025e-1,1.5807965,2.09801,-8.86e-1,-4.002,-3.32e-1]},\
  help='Dictionary of the parameter to change, default: camera_RE_gyroaverage_synchrotron_example inputs')
  return parser.parse_args()

# Execute script --------------------------------
if __name__ == "__main__":
  # data
  package_paths = ['./util/python_utils']
  # execute script
  insert_packages_paths(package_paths)
  args = generate_argument_parser()
  create_modified_example(args.example_name,args.example_dir,args.dest_dir,\
  args.name_modifier,args.separator,args.params_to_change_dict)
# -----------------------------------------------
