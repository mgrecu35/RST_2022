gfortran -c -fPIC -O3 src_f90/band_dble.f90
f2py -c -m lidarSim lidar_simulator.F90 src_f90/bisection.f90 \
     bhmie.f mieRoutines.f90 src_f90/radtran_tau_dble.f src_f90/rosen.f band_dble.o
