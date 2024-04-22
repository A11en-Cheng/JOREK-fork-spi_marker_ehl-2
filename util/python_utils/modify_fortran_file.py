# --------------------------------------------------- #
# The procedures contained in the module are used     #
# for modifying a fortran file from a input           #
# dictionary line by line. Example of its usage is    #
# the modification of particle example files for      #
# non regression testing. The method is write safe    #
# in the sens that modifications are applied only     #
# to a copy of the source file but not on the source. #
# --------------------------------------------------- #

# add ' as characters to string
# inputs: 
#   string: (string) string to be modified
# outputs:
#   string: (string) string beginning and ending with '
def fortran_compatible_string(string):
  return "'"+string+"'"

# Class containing procedures used for modifying
# FORTRAN files
class ModifyFortranFile():
  # Class constructor
  # inputs:
  #   filename: (string) the name of the file
  #   filedir:  (string) directory containing the file 
  def __init__(self,source_file,source_dir='.', \
    dest_dir='.',modify_name='_modified',separator='/'): 
    from pathlib import Path
    self.separator = separator
    self.source_path = Path("".join([source_dir,separator,source_file]))
    self.dest_path   = Path("".join([dest_dir,separator,\
    self.source_path.stem,modify_name,self.source_path.suffix]))

  # Class destructor, reset attributes to zero for safety
  def __del__(self):
    self.separator    = None
    self.source_path  = None
    self.dest_path    = None

  # Override class printer
  def __str__(self):
    return "".join(['Class for modiying FORTRAN ocdes source: ',\
    self.source_path.as_posix(),', destination: ',self.dest_path.as_posix()])

  # convert a variable in a fortran string
  # inputs:
  #   variable: (*) python variable to be converted in string
  # output: 
  #   fortran_string: (string) variable converted in string
  def convert_variable_fortran_string(self,variable):
    if(variable is dict):
      print('Warning dictionaries are not converted in fortran string!')
      return ''
    if(type(variable) is str):
      fortran_string = fortran_compatible_string(variable)
    else:
      fortran_string = variable
    fortran_string = str(fortran_string)
    fortran_string = fortran_string.replace('False','.false.')
    fortran_string = fortran_string.replace('True','.true.')
    fortran_string = fortran_string.replace('(','[')
    fortran_string = fortran_string.replace(')',']')
    fortran_string = fortran_string.replace('{','[')
    fortran_string = fortran_string.replace('}',']')
    return fortran_string

  # modify a string using regex
  # inputs:
  #   string:          (string) string to modify
  #   dict_substitute: (dict) dictionary containing the
  #                    substitution strings
  #   dict_regex:      (dict) dictionary containing the
  #                    compiled regexes
  # outputs:
  #   string: (string) modified string
  def modify_string_regex(self,string,remove_next_line,dict_substitutes,\
  dict_regex_sub,dict_regex_return,regex_return_break,regex_double_return):
    import re 
    # check if the line or part of it must be remobed
    if(remove_next_line):
      # check if there is a separator from the beginning of the line
      if(regex_return_break.match(string)):
        # remove all characters from the beginning of the line up to
        # the separaor ; and do not remove the next line line
        string = regex_return_break.sub('',string); 
        remove_next_line = False;
      # check if there is a retun at the end of the line without separators
      elif(regex_double_return.match(string)):
        string = '' # remove all line
      # there is no separator but the line is not returned at the end
      # remove all line and do not remove the next line
      else:
        string = ''; remove_next_line = False;
    # check if the line contains a variable to be modified
    if(any(key in string for key in dict_substitutes.keys())):
      # check if there is a variable to be modified at the end
      # of the string and the string is returned, if yes remove
      # the next line of the beginning of it
      remove_next_line = any(regex_return.match(string) \
      for regex_return in dict_regex_return.values())
      # modify all selected variables with the new ones
      for key,substitute in dict_substitutes.items():
        string = dict_regex_sub[key].sub(substitute,string)
    return string,remove_next_line

  # Modify variables in destination files from dictionary.
  # variables are identified looking for 'key =' substrings
  # while reading the file. Warning: the substitution 
  # operation is performed for each matching line!
  # inputs: 
  #   dict_variables: (dict) dictionary containing the variables
  #                   to be substituted
  def modify_variables(self,dict_variables):
    from string import whitespace
    import re
    remove_next_line = False
    # create dictionary of compiled regular expression and substitution strings
    regex_return_break = re.compile(r"^.*?;")
    regex_double_return = re.compile(r"^.*,&$")
    dict_compiled_regex_return = dict(("(;|\s+|^)"+key,re.compile(re.escape(key)+\
    "(\s+|)=.*&$")) for key in dict_variables.keys())
    dict_compiled_regex_sub    = dict((key,re.compile("(;|\s+|^)"+re.escape(key)+\
    "(\s+|)=.*?(;|$)")) for key in dict_variables.keys())
    dict_substitute_str = dict((key,"".join([key,'=',self.convert_variable_fortran_string(\
    variable),';'])) for key,variable in dict_variables.items()) 
    # read the line from srouce, apply the regex and write the line to destination
    # create destination if does not exit
    with self.source_path.open(mode='r') as source:
      with(self.dest_path.open(mode='w+')) as destination:
        for line in source:
          line,remove_next_line = self.modify_string_regex(line,remove_next_line,\
          dict_substitute_str,dict_compiled_regex_sub,dict_compiled_regex_return,\
          regex_return_break,regex_double_return); destination.write(line)

# Test main
if __name__ == '__main__':

  # test data
  source_file    = 'ex1.f90'
  source_dir     = './particles/examples'
  dest_dir       = './particles'
  modify_name    = 'mode'
  separator      = '+'
  dict_variables = {'p%x':1,'p%v':[1,2,3]}

  # Constructor 
  fortran_modifier = ModifyFortranFile(source_file,source_dir=source_dir,\
  dest_dir=dest_dir,modify_name=modify_name,separator=separator)  
  print(fortran_modifier)
  fortran_modifier.modify_variables(dict_variables)
  del fortran_modifier 
