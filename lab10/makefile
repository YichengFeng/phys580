main.o: percolation_2d.o lab10.cxx
	g++ percolation_2d.o lab10.cxx -o main.o `root-config --cflags --glibs`

percolation_2d.o: percolation_2d.h percolation_2d.cxx
	g++ -c percolation_2d.h percolation_2d.cxx
