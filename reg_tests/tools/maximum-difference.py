# method used for comparing and print the two datasets
def difference_norm_infinity(filename1,filename2,dsetname):
  from h5py import File as H5File
  from numpy import array,amax
  from numpy import abs as npabs
  with H5File(filename1) as h5f1:
    with H5File(filename2) as h5f2:
      dset1 = array(h5f1[dsetname])
      print(amax(npabs(array(h5f2[dsetname])-dset1),\
      axis=tuple([id for id in range(dset1.ndim)])))

# Argument parser
def generate_argument_parser():
  from argparse import ArgumentParser
  parser = ArgumentParser(description=\
  "Method for comparing the maximum difference of the same dataset between two different HDF5 files")
  parser.add_argument('--first-file-name','-fn1',type=str,required=False,\
  action="store",dest="file1_name",default="end.h5",help="Name of the first file")
  parser.add_argument('--second-file-name','-fn2',type=str,required=False,\
  action="store",dest="file2_name",default="jorek_restart.h5",help="Name of the second file")
  parser.add_argument('--dataset-name','-dn',type=str,required=False,\
  action="store",dest="dset_name",default="values",help="Name of the dataset to be compared")
  return parser.parse_args()

# script name
if __name__ == "__main__":
  args = generate_argument_parser()
  difference_norm_infinity(args.file1_name,\
  args.file2_name,args.dset_name)  
