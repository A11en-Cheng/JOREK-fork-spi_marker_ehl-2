# --- General settings
jorekmodel="303"
tor_harmonics=3
fourier_planes=4
period_harmonics=1
order=3
n_nodes=60001
n_elements=60001
n_boundary=1001
n_pieces=6001
debugoptions=""
options=""
description="Gyroaverage synchrotron image from gc distribution, model$jorekmodel, n_tor=$tor_harmonics."
mpitasks=2
requiredfiles="input jorek_restart.h5"
extra_remote_files="jorek_restart.h5"
restart_file="particle_restart.h5"
result_file="result_intensities.h5"
threshold=1.e-15
test_dataset="pixel_filter_intensities"
particle_example="camera_RE_gyroaverage_synchrotron_example"
particle_example_dir="particles/postprocessors/examples"
example_name_modifier="_modified"
example_inputs="n_frames:int:1+n_times:int:1+fields_filename:str:jorek+image_filename:str:result_intensities+n_groups:int:1+n_spectra:int:1+n_wavelengths:int:40+n_int_camera_param:int:5+n_real_camera_param:int:9+write_gc_in_txt:bool:False+particle_filenames:str:particle_restart+min_spectra:float:3e-6+max_spectra:float:3.5e-6+pinhole_positions:float:-1e0,-4e0,-3.32e-1+int_camera_param:int:1,60,60,0,1+real_camera_param:float:5.23e-1,5.23e-1,1.5707963267948966,9.998025e-1,1.5807965,2.09801,-1e0,-4e0,-3.32e-1"
binaries="${particle_example}${example_name_modifier}"
binaries_initial=""

# --- Compile the code for the test case
function compile_jorek () {
    ./util/config.sh model=$jorekmodel n_tor=$tor_harmonics \
    n_order=$order n_plane=$fourier_planes n_period=$period_harmonics \
    n_nodes_max=$n_nodes n_elements_max=$n_elements n_boundary_max=$n_boundary \
    n_pieces_max=$n_pieces $options                                             || exit 1
    python3 "${codedir}/particles/utils/create_temporary_example.py" \
    -ed "${codedir}/$particle_example_dir" -efn "${particle_example}.f90" \
    -dd $codedir -nm $example_name_modifier -dptc $example_inputs               || exit 1
    make $compilopt $debugoptions ${binaries}                                   || exit 1
    rm "${codedir}/${particle_example}${example_name_modifier}.f90"             || exit 1
}

# --- Dummy initial run
function initial_run() {
  dummy_initial_run_particles                                                   || exit 1
}

# --- Carry out the test case
function restart_run () {
  $MPIRUN $mpitasks ./$binaries < input | tee logfile                                 || exit 1
}


# --- Compare the results of the test case to the reference solution
function compare_results () {
  compare_results_generic $threshold $result_file $test_dataset || exit 1
}
