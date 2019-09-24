main.o: runge_kutta.o planet_orbit.o two_planet_orbit.o three_body_2d.o lab5.cxx
	g++ runge_kutta.o planet_orbit.o two_planet_orbit.o three_body_2d.o lab5.cxx -o main.o `root-config --cflags --glibs`

three_body_2d.o: runge_kutta.o three_body_2d.h three_body_2d.cxx
	g++ -c runge_kutta.o three_body_2d.h three_body_2d.cxx

two_planet_orbit.o: runge_kutta.o two_planet_orbit.h two_planet_orbit.cxx
	g++ -c runge_kutta.o two_planet_orbit.h two_planet_orbit.cxx

planet_orbit.o: runge_kutta.o planet_orbit.h planet_orbit.cxx
	g++ -c runge_kutta.o planet_orbit.h planet_orbit.cxx

runge_kutta.o: runge_kutta.h runge_kutta.cxx
	g++ -c runge_kutta.h runge_kutta.cxx


clean:
	rm *.o

