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
  yes | imasdb -RR ${database}  # Remove pre-existing database with that name
  imasdb ${database}      # 
  cp jorek_restart.h5 jorek00000.h5
  ./jorek2_IDS < input0                                                              || exit 1
  python imas2jorek.py -d ${database} -s 111111 -r 1 -tk inxflow                     || exit 1
}

# --- Compare the results of the test case to the reference solution
function compare_results () {
  diff jorek_namelist jorek_namelist_ref                                             || exit 1
}
