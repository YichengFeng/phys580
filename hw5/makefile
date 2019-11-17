main.o: ising_model_2d.o molecular_dynamics_2d.o hw5.cxx
	g++ ising_model_2d.o molecular_dynamics_2d.o hw5.cxx -o main.o `root-config --cflags --glibs`

molecular_dynamics_2d.o: molecular_dynamics_2d.h molecular_dynamics_2d.cxx
	g++ -c molecular_dynamics_2d.h molecular_dynamics_2d.cxx

ising_model_2d.o: ising_model_2d.h ising_model_2d.cxx
	g++ -c ising_model_2d.h ising_model_2d.cxx
