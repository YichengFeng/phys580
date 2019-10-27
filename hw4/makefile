main.o: random_walk_2d.o random_walk_3d.o diffusion_2d.o percolation_2d.o hw4.cxx
	g++ random_walk_2d.o random_walk_3d.o diffusion_2d.o percolation_2d.o hw4.cxx -o main.o `root-config --cflags --glibs`

percolation_2d.o: percolation_2d.h percolation_2d.cxx
	g++ -c percolation_2d.h percolation_2d.cxx

diffusion_2d.o: diffusion_2d.h diffusion_2d.cxx
	g++ -c diffusion_2d.h diffusion_2d.cxx

random_walk_3d.o: random_walk_3d.h random_walk_3d.cxx
	g++ -c random_walk_3d.h random_walk_3d.cxx

random_walk_2d.o: random_walk_2d.h random_walk_2d.cxx
	g++ -c random_walk_2d.h random_walk_2d.cxx
