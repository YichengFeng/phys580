main.o: square_well_1d.o match_square_well_1d.o qm_pow4.o lennard_jones_1d.o hw6.cxx
	g++ square_well_1d.o match_square_well_1d.o qm_pow4.o lennard_jones_1d.o hw6.cxx -o main.o `root-config --cflags --glibs`

square_well_1d.o: square_well_1d.h square_well_1d.cxx
	g++ -c square_well_1d.h square_well_1d.cxx

match_square_well_1d.o: match_square_well_1d.h match_square_well_1d.cxx
	g++ -c match_square_well_1d.h match_square_well_1d.cxx

qm_pow4.o: qm_pow4.h qm_pow4.cxx
	g++ -c qm_pow4.h qm_pow4.cxx

lennard_jones_1d.o: lennard_jones_1d.h lennard_jones_1d.cxx
	g++ -c lennard_jones_1d.h lennard_jones_1d.cxx
