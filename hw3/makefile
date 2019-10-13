main.o: runge_kutta.o relativity_precession.o two_planet_orbit.o capacitor_2d.o charge_dist_3d.o nintegrate_1d.o hw3.cxx
	g++ runge_kutta.o relativity_precession.o two_planet_orbit.o capacitor_2d.o charge_dist_3d.o nintegrate_1d.o hw3.cxx -o main.o `root-config --cflags --glibs`

nintegrate_1d.o: nintegrate_1d.h nintegrate_1d.cxx
	g++ -c nintegrate_1d.h nintegrate_1d.cxx

charge_dist_3d.o: charge_dist_3d.h charge_dist_3d.cxx
	g++ -c charge_dist_3d.h charge_dist_3d.cxx

capacitor_2d.o: capacitor_2d.h capacitor_2d.cxx
	g++ -c capacitor_2d.h capacitor_2d.cxx

relativity_precession.o: relativity_precession.h relativity_precession.cxx
	g++ -c relativity_precession.h relativity_precession.cxx

two_planet_orbit.o: two_planet_orbit.h two_planet_orbit.cxx
	g++ -c two_planet_orbit.h two_planet_orbit.cxx

runge_kutta.o: runge_kutta.h runge_kutta.cxx
	g++ -c runge_kutta.h runge_kutta.cxx

clean:
	rm *.o

