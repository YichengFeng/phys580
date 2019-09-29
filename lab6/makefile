main.o: runge_kutta.o planet_spin_orbit.o lab6.cxx
	g++ runge_kutta.o planet_spin_orbit.o lab6.cxx -o main.o `root-config --cflags --glibs`

planet_spin_orbit.o: planet_spin_orbit.h planet_spin_orbit.cxx
	g++ -c planet_spin_orbit.h planet_spin_orbit.cxx

runge_kutta.o: runge_kutta.h runge_kutta.cxx
	g++ -c runge_kutta.h runge_kutta.cxx

clean:
	rm *.o

