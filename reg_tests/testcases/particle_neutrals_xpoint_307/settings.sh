# --- General settings
jorekmodel="307"
description="kinetic_neutrals with model$jorekmodel using xpoint plasma"
mpitasks=2
num_threads=8
binaries="kinetic_neutrals"
binaries_initial="jorek_model${jorekmodel}_1"
extra_restart="part_restart.h5"
requiredfiles="input acd12_h.dat  acd96_h.dat  ccd12_h.dat  ccd96_h.dat plt12_h.dat  plt96_h.dat  prb12_h.dat  prb96_h.dat  prc12_h.dat  prc96_h.dat  scd12_h.dat scd96_h.dat y_DD.dat  ye_DD.dat"
extra_remote_files="acd12_h.dat  acd96_h.dat  ccd12_h.dat  ccd96_h.dat plt12_h.dat  plt96_h.dat  prb12_h.dat  prb96_h.dat  prc12_h.dat  prc96_h.dat  scd12_h.dat scd96_h.dat y_DD.dat  ye_DD.dat part_restart.h5"
particle_example="kinetic_neutrals"
particle_example_dir="particles/examples"


# --- Compile the code for the test case
function compile_jorek () {
  if [ "$initialrun" == "yes" ]; then
    ./util/config.sh model=$jorekmodel n_tor=1 n_coord_tor=1 l_pol_domm=0 n_plane=1 n_period=1 n_coord_period=1 || exit 1
    make $compilopt $debugoptions jorek_model${jorekmodel}                                                      || exit 1
    mv jorek_model${jorekmodel} jorek_model${jorekmodel}_1                                                      || exit 1
    make cleanall                                                                                               || exit 1
  fi
  ./util/config.sh model=$jorekmodel n_tor=1 n_coord_tor=1 l_pol_domm=0 n_plane=1 n_period=1 n_coord_period=1   || exit 1
  make $compilopt $debugoptions kinetic_neutrals                                                                || exit 1
}


# --- Initial run only required when preparing or updating the test case
function initial_run () {
  ${codedir}/util/setinput.sh input restart_particles=.f. nstep_n=10,10,10 tstep_n=10.,100.,1000. use_ncs=.f.   || exit 1
  $MPIRUN 1 ./jorek_model${jorekmodel}_1 < input | tee logfile_initial                                          || exit 1
  ${codedir}/util/setinput.sh input restart_particles=.f. nstep_n=10 tstep_n=10. use_ncs=.t.                    || exit 1
  export OMP_NUM_THREADS=$num_threads                                                                           || exit 1
  echo "setting OMP_NUM_THREADS=$num_threads, due to the requirements of the test"                              || exit 1
  $MPIRUN $mpitasks ./kinetic_neutrals < input | tee logfile_initial2                                           || exit 1
}


# --- Carry out the test case
function restart_run () {
  ${codedir}/util/setinput.sh input restart_particles=.t. nstep_n=1 tstep_n=10.                                 || exit 1
  export OMP_NUM_THREADS=$num_threads                                                                           || exit 1
  echo "setting OMP_NUM_THREADS=$num_threads, due to the requirements of the test"                              || exit 1
  $MPIRUN $mpitasks ./kinetic_neutrals < input | tee logfile                                                    || exit 1
}


# --- Compare the results of the test case to the reference solution
function compare_results () {
  compare_results_generic 5.e-7                                                                                 || exit 1
}
