gfortran -o reweight.x reweight_2d.f90
ln -s tass_0.5/COLVAR COLVAR_1
ln -s tass_0.5/HILLS HILLS_1
./reweight.x -T0 300.0 -T 2700 -dT 2700.0  -grid 0.5  10.0 0.5 5.0 13.0 0.08 -pfrqMD 50 -dtMTD 500 -nr 20  
