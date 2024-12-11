# --- General settings
jorekmodel="600"
options="with_vpar=.true. with_TiTe=.false. with_neutrals=.false. with_impurities=.false. with_refluid=.false."
description="Exports restart file to IMAS, then reads it from IMAS to create a JOREK input file and compares it with a reference."
mpitasks=1
binaries="jorek2_IDS" 
binaries_initial="jorek2_IDS" 
python_scripts="./communication/IMAS/imas2jorek.py" 
requiredfiles="input0 imas.nml jorek_namelist_ref"
extra_remote_files=""

# --- Compile the code for the test case
function compile_jorek () {
  ./util/config.sh model=$jorekmodel n_tor=1 n_plane=1 n_period=1 $options           || exit 1
  make $compilopt $debugoptions jorek2_IDS                                           || exit 1
}

# --- Initial run only required when preparing or updating the test case
function initial_run () {
  echo "Initial run is not necessary for this case"                                  || exit 1
}

# --- Carry out the test case
function restart_run () {
  database="imas_regtest_db"
  if [ -d ~/public/imasdb/${database}/4/111111/1/ ]; then
      echo "Directory exists."
  else
      pwd
      echo "directory does not exist"
  fi
 
  cp jorek_restart.h5 jorek00000.h5
  ./jorek2_IDS < input0                                                              || exit 1
  python imas2jorek.py -d ${database} -p 111111 -r 1 -dd 4 -tk inxflow               || exit 1
  if [ -d ~/public/imasdb/${database}/4/111111/1/ ]; then
      rm -rf ~/public/imasdb/${database}/4/111111/1/
  fi
      
}

# --- Compare the results of the test case to the reference solution
function compare_results () {
    sed '/^\s*\*/d' jorek_namelist > file1
    sed '/^\s*\*/d' jorek_namelist_ref > file2  
  diff file1 file2                                             || exit 1
}
