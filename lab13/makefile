main.o: square_well_1d.o lennard_jones_1d.o lab13.cxx
	g++ square_well_1d.o lennard_jones_1d.o lab13.cxx -o main.o `root-config --cflags --glibs`

square_well_1d.o: square_well_1d.h square_well_1d.cxx
	g++ -c square_well_1d.h square_well_1d.cxx

lennard_jones_1d.o: lennard_jones_1d.h lennard_jones_1d.cxx
	g++ -c lennard_jones_1d.h lennard_jones_1d.cxx
