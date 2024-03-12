# ----------------------------------------------------- #
# Class used for parsing a dictionary as argument       #
# using the pyhton argument parser action.              #
# The syntax of the arguments must be:                  #
# key1:type1:value11,value12,value13+key2:type2:value21 #
# For using the action, the call of the                 #
# add_argument function must contain the argument:      #
# (...,action=ParseDictionary,...)                      #
# ----------------------------------------------------- #
from argparse import Action

# Cast a list of string in a list of desired types 
# inputs:
#   type_str: (string) string defining the casting type
#   value:    (string) string containing the values
# outputs:
#   value: value with the type defined by type string
def cast_type(type_str,value):
  if(type_str=='int'):
    return int(value)
  elif(type_str=='float'):
    return float(value)
  elif(type_str=='ord'):
    return ord(value)
  elif(type_str=='hex'):
    return hex(value)
  elif(type_str=='oct'):
    return oct(value)
  elif(type_str=='bin'):
    return bin(value)
  elif(type_str=='bool' or type_str=='logical'):
    if(value=='True' or value=='true'):
      return True
    else:
      return False
  elif(type_str=='str'):
    return str(value)
  else:
    print(''.join(['Unrecognised type ',type_str,\
    ' for casting for value: ',value,'!']))
    return value

# Cast a list of string in a list of desired types 
# inputs:
#   type_str:  (string) string defining the casting type
#   value:     (string) string containing the values
#   delimiter: (string) delimiter for string splitting
# outputs:
#   value: value with the type defined by type string
def type_cast_list(type_str,values,delimiter):
  return [cast_type(type_str,value) for value in \
  values.split(delimiter)]

# Class defininf the action for parsing dictionaries
class ParseDictionary(Action):
  # defining the initialisation of the action for parsing 
  # a dictionary
  def __init__(self,option_strings,dest,nargs=None,\
    dict_delimiter='+',entries_delimiter=':',\
    values_delimiter=',',*args,**kwargs):
    if(nargs is not None):
      print('Warning dictionary argument parser: only one dictionary is parsed!')
    super().__init__(option_strings,dest,**kwargs)
    self.dict_delimiter=dict_delimiter
    self.values_delimiter=values_delimiter
    self.entries_delimiter = entries_delimiter

  # Defining the action for parsing a dictionary
  # its extension with multiple inputs requires to add
  # the args argument into the call inputs
  def __call__(self,parser,namespace,entries,option_string=None):
    setattr(namespace,self.dest,dict())
    if(self.dict_delimiter in entries):
      entries = entries.split(self.dict_delimiter)
    else:
      entries = [entries]
    for entry in entries:
      key,type_str,values = entry.split(self.entries_delimiter)
      if(self.values_delimiter in values): 
        getattr(namespace,self.dest)[key] = type_cast_list(\
        type_str,values,self.values_delimiter)
      else:
        getattr(namespace,self.dest)[key] = cast_type(type_str,values)

 # Main for testing      
if __name__ == '__main__':

  from argparse import ArgumentParser

  parser = ArgumentParser()
  parser.add_argument('-d','--dictionary',dest='dictionary',nargs='*',\
  action=ParseDictionary,required=None,default={'values':[10,4]})
  args = parser.parse_args()
  print(args.dictionary)
