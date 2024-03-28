# --------------------------------------------------------------- #
# Script aimed to the automatic generation, compilation and       #
# execution of unit tests. Unit test results are synthesized in   #
# a XML file compatible with the Atlassian Bamboo CI/CD installed #
# on the ITER HPC infrastructure.                                 #
# --------------------------------------------------------------- #

# general purpose routines -------------------------------------- #

# Check for launchers in environment variables
# inputs:
#   launchers: (dict(string)) launchers to be invoked for 
#              executing a unit test application
# outputs:
#   launchers: (dict(string)) launchers to be invoked for 
#              executing a unit test application
def get_launchers_from_environment(launchers):
  from os import environ
  from os import getenv
  if("MPIRUN" in environ):
    print("Found variable for MPI launcher: override!")
    launchers["mpi"] = getenv("MPIRUN")+"2 "
    if("SERIALRUN" in environ):
      launchers["mpi"] = launchers["mpi"]+getenv("SERIALRUN")
    else:
      launchers["mpi"] = launchers["mpi"]+"./"
    print(launchers["mpi"])
  if("SERIALRUN" in environ):
    print("Found variable for MPI launcher: override!")
    launchers["serial"] = getenv("SERIALRUN")
  return launchers

# Append path of packages
# inputs:
#   package_path: (list,string) package paths to add
#   position:     (integer) position in the path list  
def insert_packages_paths(package_paths,position=0):
  import os; import sys;
  for package_path in package_paths:
    sys.path.insert(position,os.path.abspath(package_path))

# return a boolean from a string
# inputs:
#   string: (character) string to be converted in bool,
#           it must be: true,True,False,false
def return_bool_from_string(string):
  if((string=='True')or(string=='true')):
    return True
  elif((string=='False')or(string=='false')):
    return False
  else:
    raise ValueError("Error string cannot be converted to bool, not: True,true,False,false")

# remove file if requested
# inputs:
#   filename: (string) filename of file to be removed
#   remove:   (boolean) if true file is removed
def remove_file(filename,remove):
  from os import system
  if(remove):
    system(''.join(['rm -f ',filename]))

# create a dictionary of empty list from keys
# inputs:
#   keys_list:  (list,strings) list of keys
# outputs:
#   list_dict:  (dict,list) dictionary of empty lists
def create_list_dictionary_from_keys(keys_list):
  list_dict = {}
  for key in list(set(keys_list)):
    list_dict[key] = []
  return list_dict

# routines for handling unit test drivers ----------------------- #

# write the driver file for fruit serial unit tests
# inputs:
#   driver_path:        (path)   path to the unit test driver file
#   test_name:          (string) basename of the unit test
#   test_basket_prefix: (string) prefix of the basket subroutine
#   driver_suffix:      (string) suffix of the unit test driver
#   test_prefix:        (string) prefix of the unit test module
#   test_suffix:        (string) suffix of the unit test module
#   log_fruit_summary:  (bool) if true, the fruit summary is logged 
def write_test_driver_serial(driver_path,test_name,test_basket_prefix,\
driver_suffix,test_prefix,test_suffix,log_fruit_summary):
  with driver_path.open(mode='w') as driver:
    driver.write(''.join(['program ',test_name,driver_suffix,'\n']))
    driver.write('use fruit\n')
    driver.write(''.join(['use ',''.join([test_prefix,test_name,\
    test_suffix,', only: ',test_basket_prefix,test_name,'\n'])]))
    driver.write('  implicit none\n')
    driver.write('\n')
    driver.write('  ! init fruit suite\n')
    driver.write('  call init_fruit_xml\n')
    driver.write('  call init_fruit\n')
    driver.write('  call fruit_hide_dots')
    driver.write('\n')
    driver.write(''.join(['  ! run ',test_name,' test basket\n']))
    driver.write(''.join(['  call ',test_basket_prefix,test_name,'\n']))
    driver.write('\n')
    driver.write('  ! write test summary and finalize test suit\n')
    driver.write('  call fruit_summary_xml\n')
    if(log_fruit_summary):
      driver.write('  call fruit_summary\n')
    driver.write('  call fruit_finalize\n')
    driver.write(''.join(['end program ',test_name,driver_suffix,'\n']))

# write the driver file for fruit mpi-enabled unit tests
# inputs:
#   driver_path:        (path)   path to the unit test driver file
#   test_name:          (string) basename of the unit test
#   test_basket_prefix: (string) prefix of the basket subroutine
#   driver_suffix:      (string) suffix of the unit test driver
#   test_prefix:        (string) prefix of the unit test module
#   test_suffix:        (string) suffix of the unit test module
#   log_fruit_summary:  (bool) if true, the fruit summary is logged
def write_test_driver_parallel(driver_path,test_name,test_basket_prefix,\
driver_suffix,test_prefix,test_suffix,log_fruit_summary):
  with driver_path.open(mode='w') as driver:
     driver.write(''.join(['program ',test_name,driver_suffix,'\n']))
     driver.write('use fruit\n')
     driver.write('use fruit_mpi\n')   
     driver.write('use mod_mpi_tools, only: init_mpi_threads\n')
     driver.write('use mod_mpi_tools, only: finalize_mpi_threads\n')
     driver.write(''.join(['use ',''.join([test_prefix,test_name,\
     test_suffix,', only: ',test_basket_prefix,test_name,'\n'])]))
     driver.write('  implicit none\n')
     driver.write('  integer :: rank,n_tasks,ifail\n')
     driver.write('\n')
     driver.write('  ! init the mpi communicator\n')
     driver.write('  call init_mpi_threads(rank,n_tasks,ifail)\n')
     driver.write('\n')
     driver.write('  ! init the fruit suit\n')
     driver.write('  call init_fruit\n')
     driver.write('  call fruit_init_mpi_xml(rank)\n')
     driver.write('  call fruit_hide_dots')
     driver.write('\n')
     driver.write(''.join(['  ! run ',test_name,' test basket\n']))
     driver.write(''.join(['  call ',test_basket_prefix,test_name,\
     '(rank,n_tasks,ifail)','\n']))
     driver.write('\n')  
     driver.write('  ! write test summary and finalize test suit\n')
     if(log_fruit_summary):
       driver.write('  call fruit_summary_mpi(n_tasks,rank)\n')
     driver.write('  call fruit_summary_mpi_xml(n_tasks,rank)\n')
     driver.write('  call fruit_finalize_mpi(n_tasks,rank)\n')
     driver.write('\n')
     driver.write('  ! finalize MPI communicator\n')
     driver.write('  call finalize_mpi_threads(ifail)\n')
     driver.write(''.join(['end program ',test_name,driver_suffix,'\n']))

# reorder and store path to be tested in dictionary
# inputs:
#   path_list: (list(string)) string of path to be reordered
#   test_parallel: (list,string) list of the parallelization types 
#   path_dict: (dict) initialised dictionary to be filled with path
# output:
#   path_dict: (dict) filled dictionary to be filled with path
def reorder_store_test_paths_in_dict(path_list,test_parallel,path_dict):
  for posix_path in path_list:
    if(posix_path.is_file()):
      for key in test_parallel:
        if(key in posix_path.name):
          if(posix_path not in path_dict[key]):
            path_dict[key].append(posix_path)
        else:
          if(posix_path not in path_dict['serial']):
            path_dict['serial'].append(posix_path)
    else:
      print("".join(['Warning: ',posix_path.name,' is not a file: skip!']))
  return path_dict

# reorder and store user selected unit tests to be run 
# inputs:
#   test_list:     (list,string) list relative paths of unit test
#                                to be run from root folder 
#   test_parallel: (list,string) list of the parallelization types 
# outputs:
#   test_modules: (dict) dictionary containing the path associated
#                 to the test modules
def reorder_store_unit_tests(test_list,test_parallel):
  from pathlib import Path
  # generate path list
  path_list = [Path("".join(['./',filepath])) for filepath in test_list]
  # initialize test module dictionary
  test_parallel.append('serial')
  test_modules = create_list_dictionary_from_keys(test_parallel)
  test_parallel.remove('serial')
  # loop on the path to be stored
  for test in path_list:
    test_modules = reorder_store_test_paths_in_dict(path_list,\
    test_parallel,test_modules) 
  return test_modules

# Find unit test modules in a given folder. Unit test modules are
# identified via the prefix mod_, the suffix _test and the
# extension .f90
# inputs:
#   test_dirs:     (list)(string) list paths to the search directories
#   test_prefix:   (string) prefix identifying a test module
#   test_suffix:   (string) suffix identifying a test module
#   test_ext:      (string) file extension identifying a test module
#   test_parallel: (list,string) list of the parallelization types 
# outputs:
#   test_modules: (dict) dictionary containing the path associated
#                 to the test modules
def find_unit_test_modules(test_dirs,test_prefix,test_suffix,\
test_ext,test_parallel):
  from pathlib import Path
  # initialize test module dictionary
  test_parallel.append('serial')
  test_modules = create_list_dictionary_from_keys(test_parallel)
  test_parallel.remove('serial')
  # loop on all directories in which perform the search
  for test_dir in test_dirs:
    # loop on the paths for all test modules in test_dir 
    # and store them in the test_modules dictionary as a 
    # function of their parallelism
    test_dir_path = Path(test_dir)
    test_modules = reorder_store_test_paths_in_dict(\
    test_dir_path.glob("".join([test_prefix,'*',test_suffix,\
    '.',test_ext])),test_parallel,test_modules)
  # return the module dictionary   
  return test_modules

# generate a unit test driver from a test module path
# inputs:
#   test_path:          (path) path posix of the test module
#   test_dir:           (string) unit test directory
#   test_basket_prefix: (string) prefix of the unit test basket
#   test_ext:           (string) test extension
#   driver_suffix:      (string) suffix of the unit test driver
#   log_fruit_summary:  (bool) if true, the fruit summary is logged 
# outputs:
#   driver_path:        (path) path posix of the test driver
def generate_unit_test_driver(test_path,test_basket_prefix,\
test_prefix,test_suffix,test_ext,driver_suffix,log_fruit_summary):
  from pathlib import Path
  # create the unit test driver path
  test_name = test_path.name.replace(test_prefix,'').replace(\
  ''.join([test_suffix,'.',test_ext]),'')
  driver_path = Path(''.join([test_name,driver_suffix,'.',test_ext]))
  # force remove old file and create new one
  driver_path.unlink(missing_ok=True)
  driver_path.touch()
  # write unit test driver as a function of the parallelism
  if('mpi' in test_name):
    write_test_driver_parallel(driver_path,test_name,test_basket_prefix,\
    driver_suffix,test_prefix,test_suffix,log_fruit_summary)
  else:
    write_test_driver_serial(driver_path,test_name,test_basket_prefix,\
    driver_suffix,test_prefix,test_suffix,log_fruit_summary)
  return driver_path

# cean all jorek compilation files:
def make_cleanall():
  from os import system
  system('make cleanall')

# compile the unit test driver using 8 threads
# inputs:
#   driver_path:   (path) path posix of the test driver
#   debug_options: (string) define which set of pre-defined
#                  debug options to use for compilation
def compile_unit_test_driver(driver_path,debug_options):
  from os import system
  error = 0
  exec_name = driver_path.name.replace(driver_path.suffix,'')
  # check if the driver is valid
  if(driver_path.is_file()):
    system(''.join(['rm -f ',exec_name]))
    system(''.join(['make -j8 ',debug_options,\
    ' ',exec_name]))
  else:
    print('Warning: driver is not a valid file!')
    error = 0
  return error

# find the right launcher for the application
# inputs:
#   exec_name: (string) name of the application
def find_launcher_type(exec_name):
  if('mpi' in exec_name):
    return 'mpi'
  else:
    return 'serial'

# run the unit test driver using 2 mpi tasks maximum 
# and 2 omp threads maximum
# inputs:
#   driver_path:   (path) path posix of the test driver
#   launchers:     (dict(string)) launchers to be invoked for 
#                  executing a unit test application
def run_unit_test_driver(driver_path,launchers):
  from pathlib import Path
  from os import system,environ
  # set error and the number of OMP threads
  error = 0; omp_num_threads_old = str(1);
  if('OMP_NUM_THREADS' in environ):
    omp_num_threads_old = environ['OMP_NUM_THREADS']
  system('export OMP_NUM_THREADS=2')
  exec_name = driver_path.name.replace(driver_path.suffix,'')
  # check if the exec exists and is a file
  exec_path = Path(exec_name)
  if(exec_path.is_file()):
    system(''.join([launchers[find_launcher_type(exec_name)],exec_name]))
    # remove executable
    system(''.join(['rm -f ',exec_name]))
  else:
    print('Warning: executable is not a valid file!')
    error = 1
  # restore original number of omp threads
  system(''.join(['export OMP_NUM_THREADS=',omp_num_threads_old]))
  return error

# routines for handling FRUIT XML files ------------------------- #

# Read multiple fruit files and store their path in list. Discard
# temporary result file.
# inputs:
#   result_dir:    (string) directory containing the FRUIT results
#   result_prefix: (string) root of the FRUIT result filenames
#   result_ext:    (string) extension of the FRUIT result filenames
#   remove_result: (bool) if True temporary result files are removed  
# outputs:
#   result_paths:  (list,path) path list of all FRUIT result files
def find_fruit_result_file(result_dir,result_prefix,\
result_ext,remove_result):
  from pathlib import Path
  result_dir_path = Path(result_dir)
  result_paths = []
  for result in result_dir_path.glob("".join(\
  [result_prefix,'*','.',result_ext])):
    if('_tmp' not in result.name):
      result_paths.append(result)
    else:
      remove_file(result.name,remove_result)
  return result_paths

# Read the file containing the test results as run by FRUIT.
# The total number of successes and failures is returned.
# inputs:
#   result_path:   (path) path of the FRUIT result file
#   result_map:    (dict) dictionary containing the association
#                         between the fruit errors,tests,failures,id
#                         fields and the index of the integer in list
#   remove_result: (bool) if True result files are removed
# outputs:
#   n_successes: (int) number of test ending in success
#   n_failures:  (int) number of test ending in failures
#   n_errors:    (int) number of test ending in errors
#   test_id:     (int) test id
def read_fruit_result_file(result_path,result_map,remove_result):
  from re import findall
  with result_path.open(mode='r') as result:
    lines = result.readlines()
  # extract the line containing the results
  for line in lines:
    if('failures' in line):
      results = [int(result) for result in findall(r'\d+',line)]
      break
  # remove result file if required
  remove_file(result_path.name,remove_result)
  # return results
  return results[result_map['tests']]-results[result_map['failures']],\
  results[result_map['failures']],results[result_map['errors']],\
  results[result_map['id']]

# Read and reduce all fruit results (more than one result file can
# exist for the same FRUIT test if run using MPI)
# inputs:
#   result_dir:     (string) directory containing the FRUIT results
#   result_prefix:  (string) root of the FRUIT result filenames
#   result_ext:     (string) extension of the FRUIT result filenames
#   result_map:     (dict) dictionary containing the association
#                          between the fruit errors,tests,failures,id
#                          fields and the index of the integer in list
#   remove_results: (boolean) if True result files are removed
# outputs:
#   n_successes_tot: (int) total number of test ending in success
#   n_failures_tot:  (int) total number of test ending in failures
#   n_errors_tot:    (int) total number of test ending in errors
def read_reduce_fruit_results(result_dir,result_prefix,result_ext,\
result_map,remove_results):
  # initialisations
  n_successes_tot=0; n_failures_tot=0; n_errors_tot=0;
  # find all FRUIT result paths
  result_path_list = find_fruit_result_file(result_dir,\
  result_prefix,result_ext,remove_results)
  # read and reduce fruit results, if no result is found, return an error
  if(len(result_path_list)<1):
    n_errors_tot = n_errors_tot + 1
  else:
    for result_path in result_path_list:
      n_successes,n_failures,n_errors,test_id = read_fruit_result_file(\
      result_path,result_map,remove_results)
      n_successes_tot = n_successes_tot + n_successes
      n_failures_tot = n_failures_tot + n_failures
      n_errors_tot = n_errors_tot + n_errors    
  return n_successes_tot,n_failures_tot,n_errors_tot

# global unit test routines ------------------------------------- #

# check results, if there is a non-zero number of failures
# or a non zero number of erros, the script return the 
# error exit code of 1. The success exit code of 0 is 
# returned otherwise
# inputs: 
#   n_failures: (integer) total number of failed tests
#   n_errors:   (integer) total number of errors
def check_results(n_failures,n_errors):
  if(n_failures!=0):
    print(''.join(['Unit test failed: N# failures: ',\
    str(n_failures),' N# errors: ',str(n_errors)]))
    return 1
  if(n_errors!=0):
    print(''.join(['Unit test failed: N# failures: ',\
    str(n_failures),' N# errors: ',str(n_errors)]))
    return 1
  else: 
    print('Unit test successfully completed!')
    return 0

# execute unit test: generate, compile, run unit test driver,
# read and reduce fruit result for the unit test
#   test_path:          (path) path posix of the test module
#   test_basket_prefix: (string) prefix of the unit test basket
#   test_ext:           (string) test extension
#   driver_suffix:      (string) suffix of the unit test driver
#   result_dir:         (string) directory containing the FRUIT results
#   result_prefix:      (string) root of the FRUIT result filenames
#   result_ext:         (string) extension of the FRUIT result filenames
#   result_map:         (dict) dictionary containing the association
#                              between the fruit errors,tests,failures,id
#                              fields and the index of the integer in list
#   remove_driver:      (boolean) if true the driver file
#   remove_results:     (boolean) if True result files are removed
#   log_fruit_summary:  (bool) if true, the fruit summary is logged 
#   debug_options:      (string) define which set of pre-defined
#                       debug options to use for compilation
def execute_unit_test(test_path,test_basket_prefix,test_prefix,\
test_suffix,test_ext,driver_suffix,launchers,result_dir,\
result_prefix,result_ext,result_map,remove_driver,remove_results,\
log_fruit_summary,debug_options):
  # generate the unit test driver
  driver_path = generate_unit_test_driver(test_path,test_basket_prefix,\
  test_prefix,test_suffix,test_ext,driver_suffix,log_fruit_summary) 
  # compile the unit test driver
  error_driver = compile_unit_test_driver(driver_path,debug_options)
  # execute the unit test driver
  error_exec = run_unit_test_driver(driver_path,launchers)
  # remove the driver
  remove_file(driver_path.name,remove_driver)
  # find, read and reduce all FRUIT results and remove result files
  n_successes,n_failures,n_errors = read_reduce_fruit_results(\
  result_dir,result_prefix,result_ext,result_map,remove_results)
  n_errors = n_errors + error_driver + error_exec
  return n_successes,n_failures,n_errors

# log unit tests results
# inputs: 
#   unit_test_log: (list(characters)) path of the unit tests to be logged
#   log_header:    (character) header of the unit test log to be printed
def log_unit_test_results(unit_test_log,log_header):
  if(len(unit_test_log)>0):
    print(log_header)
    for unit_test in unit_test_log:
      print("".join(['  ',unit_test]))

# execute the overall suite of unit tests
# inputs:
#   test_dirs:             (list(string)) list of paths to the search directories
#   test_parallel:         (list(string)) list of the parallelization types 
#   test_basket_prefix:    (string) prefix of the unit test basket
#   test_ext:              (string) test extension
#   driver_suffix:         (string) suffix of the unit test driver
#   result_dir:            (string) directory containing the FRUIT results
#   result_prefix:         (string) root of the FRUIT result filenames
#   result_ext:            (string) extension of the FRUIT result filenames
#   result_map:            (dict) dictionary containing the association
#                                 between the fruit errors,tests,failures,id
#                                 fields and the index of the integer in list
#   remove_driver:         (boolean) if true the driver file
#   remove_results:        (boolean) if True result files are removed
#   log_fruit_summary:     (bool) if true, the fruit summary is logged 
#   test_module_to_be_run: (list(string)) list of test module paths selected
#                                         by the user to be run
#   debug_options:         (string) define which set of pre-defined
#                          debug options to use for compilation
# outputs:
#   exit_code:          (integer) 0 if all tests terminated successfully
#                                 1 otherwise
def execute_all_unit_tests(test_dirs,test_parallel,test_basket_prefix,\
test_prefix,test_suffix,test_ext,driver_suffix,launchers,result_dir,\
result_prefix,result_ext,result_map,remove_driver,remove_results,\
log_fruit_summary,test_modules_to_be_run=[],debug_options=0):
  # check for launchers defined by environment variables
  launchers = get_launchers_from_environment(launchers)
  # initialise the failure and error counters
  n_failures=0; n_errors=0; failed_tests=[]; error_tests=[];
  if(len(test_modules_to_be_run)>0):
    # re-order and store unit tests to be run
    test_modules = reorder_store_unit_tests(test_modules_to_be_run,test_parallel)
  else:
    # find all unit test modules in test_dir
    test_modules =  find_unit_test_modules(test_dirs,test_prefix,\
    test_suffix,test_ext,test_parallel)
  # remove all jorek compilation files
  make_cleanall()
  # execute all unit tests in the test modules
  for tests in test_modules.values():
    for test in tests:
      # execute a unit test
      n_successes_loc,n_failures_loc,n_errors_loc = \
      execute_unit_test(test,test_basket_prefix,test_prefix,\
      test_suffix,test_ext,driver_suffix,launchers,result_dir,\
      result_prefix,result_ext,result_map,remove_driver,\
      remove_results,log_fruit_summary,debug_options)
      # store the name of failed tests
      if(n_failures_loc!=0): 
        failed_tests.append(test.name)
      # store the name of tests with errors
      if(n_errors_loc!=0):
        error_tests.append(test.name)
      # reduce the total number of successes, failures and errors
      n_failures  = n_failures + n_failures_loc
      n_errors    = n_errors + n_errors_loc
  # log tests with failures and/or errors
  log_unit_test_results(failed_tests,'Unit test module with failed tests:')
  log_unit_test_results(error_tests,'Unit test module with errors:')
  # check the validity of the overall result
  return check_results(n_failures,n_errors)
     
# Argument parsers ---------------------------------------------- #
# outpus:
#   parser: (namespace) namespace having the inputs as attributes
def generate_argument_parser(dict_path='./util/python_utils' ):
  from argparse import ArgumentParser 
  # add dictionary path
  insert_packages_paths([dict_path])
  from dictionary_parser_class import ParseDictionary
  # define parser
  parser = ArgumentParser(\
  description='generate, compile and execute JOREK unit tests')
  parser.add_argument('--directories','-d',type=str,nargs='*',\
  required=False,action='store',dest='test_dirs',\
  default=['./particles/tests','./elements/tests','./diagnostics/tests',\
  './communication/tests','./core/tests','./grids/tests',\
  './tools/tests','./particles/postprocessors/tests'],\
  help='relative paths of the directories containing unit tests')
  parser.add_argument('--parallelisms','-p',type=str,nargs='*',\
  required=False,action='store',dest='test_parallel',\
  default=['mpi'],help='type of parallelism of the unit tests,default: mpi')
  parser.add_argument('--fruit-basket-prefix','-fbp',type=str,\
  required=False,action='store',dest='test_basket_prefix',default='run_fruit_',\
  help='prefix of the fruit basket to be run,default: run_fruit_')
  parser.add_argument('--test-prefix','-tp',type=str,required=False,\
  action='store',dest='test_prefix',default='mod_',\
  help='prefix of the unit test module file, default: mod_')
  parser.add_argument('--test-suffix','-ts',type=str,required=False,\
  action='store',dest='test_suffix',default='_test',\
  help='suffix of the unit test module file, default: _test') 
  parser.add_argument('--test-extension','-te',type=str,required=False,\
  action='store',dest='test_ext',default='f90',\
  help='extension of the unit test module file, default: f90')
  parser.add_argument('--test-driver-suffix','-ds',type=str,required=False,\
  action='store',dest='driver_suffix',default='_test_driver',\
  help='suffix of unit test driver file, default: _test_driver')
  parser.add_argument('--result-dir','-rd',type=str,required=False,\
  action='store',dest='result_dir',default='.',\
  help='folder of the unit test results, default: .')
  parser.add_argument('--result-prefix','-rp',type=str,required=False,\
  action='store',dest='result_prefix',default='result',\
  help='prefix of the test result file, default: result')
  parser.add_argument('--result-extension','-re',type=str,required=False,\
  action='store',dest='result_ext',default='xml',\
  help='extension of the test result file, default: xml')
  parser.add_argument('--remove-driver','-rmd',type=str,required=False,\
  action='store',dest='remove_drivers',default='True',\
  help='if true the test drivers are removed after execution, default: true')
  parser.add_argument('--remove-results','-rmr',type=str,required=False,\
  action='store',dest='remove_results',default='True',\
  help='if true the test results are removed after execution, default: true')
  parser.add_argument('--log-fruit-summary','-lfs',type=str,required=False,\
  action='store',dest='log_fruit_summary',default='False',\
  help='if true the fruit summary is logged, default: false')
  parser.add_argument('--fruit-result-map','-frm',required=False,action=ParseDictionary,\
  dest='fruit_result_map',default={'errors':0,'tests':1,'failures':2,'id':3},
  help='relative position of fruit results as read by result file, default: errors:0,tests:1,failures:2,id:3')
  parser.add_argument('--launchers','-l',required=False,action=ParseDictionary,
  dest='launchers',default={'serial':'./','mpi':'mpirun -np 2 ./'},\
  help='launcher to be used for executing a test in the order: default: serial:./,mpi:mpirun -np 2]')
  parser.add_argument('--test_to_run','-ttr',type=str,nargs='*',\
  required=False,action='store',dest='tests_to_run',default=[],\
  help='path from root flder to the tests to be run, default: []')
  parser.add_argument('--debug-options','-debug',type=str,\
  required=False,action='store',dest='debug_options',default='',\
  help='integer defining the compilation debug option, default: empty string')
  return parser.parse_args()

# Execute script ------------------------------------------------ #
if __name__ == '__main__':
  from sys import exit as sysexit
  # parse the inputs
  args = generate_argument_parser()
  # run all unit test suites
  exit_code = execute_all_unit_tests(args.test_dirs,args.test_parallel,\
  args.test_basket_prefix,args.test_prefix,args.test_suffix,args.test_ext,\
  args.driver_suffix,args.launchers,args.result_dir,args.result_prefix,\
  args.result_ext,args.fruit_result_map,return_bool_from_string(args.remove_drivers),\
  return_bool_from_string(args.remove_results),return_bool_from_string(args.log_fruit_summary),\
  test_modules_to_be_run=args.tests_to_run,debug_options=args.debug_options) 
  # exit with the appropriate exit code
  sysexit(exit_code)

# End-of-the-scripts -------------------------------------------- #
