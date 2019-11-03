main.o: ising_model_2d.o lab11.cxx
	g++ ising_model_2d.o lab11.cxx -o main.o `root-config --cflags --glibs`

ising_model_2d.o: ising_model_2d.h ising_model_2d.cxx
	g++ -c ising_model_2d.h ising_model_2d.cxx
