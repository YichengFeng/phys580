main.o: molecular_dynamics_2d.o lab12.cxx
	g++ molecular_dynamics_2d.o lab12.cxx -o main.o `root-config --cflags --glibs`

molecular_dynamics_2d.o: molecular_dynamics_2d.h molecular_dynamics_2d.cxx
	g++ -c molecular_dynamics_2d.h molecular_dynamics_2d.cxx
