main.o: monte_carlo_integrate.o random_walk_2d.o lab9.cxx
	g++ monte_carlo_integrate.o lab9.cxx random_walk_2d.o -o main.o `root-config --cflags --glibs`

random_walk_2d.o: random_walk_2d.h random_walk_2d.cxx
	g++ -c random_walk_2d.h random_walk_2d.cxx

monte_carlo_integrate.o: monte_carlo_integrate.h monte_carlo_integrate.cxx
	g++ -c monte_carlo_integrate.h monte_carlo_integrate.cxx
