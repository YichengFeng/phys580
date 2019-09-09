main.o: runge_kutta.o sequential_decay.o bicycling.o projectile.o projectile_3d.o hw1.cxx
	g++ runge_kutta.o sequential_decay.o bicycling.o projectile.o projectile_3d.o hw1.cxx -o main.o `root-config --cflags --glibs`

sequential_decay.o: runge_kutta.o sequential_decay.h sequential_decay.cxx
	g++ -c runge_kutta.o sequential_decay.h sequential_decay.cxx

bicycling.o: runge_kutta.o bicycling.h bicycling.cxx
	g++ -c runge_kutta.o bicycling.h bicycling.cxx

projectile.o: runge_kutta.o projectile.h projectile.cxx
	g++ -c runge_kutta.o projectile.h projectile.cxx

projectile_3d.o: runge_kutta.o projectile_3d.h projectile_3d.cxx
	g++ -c runge_kutta.o projectile_3d.h projectile_3d.cxx

runge_kutta.o: runge_kutta.h runge_kutta.cxx
	g++ -c runge_kutta.h runge_kutta.cxx


clean:
	rm *.o


#main.o: runge_kutta.h runge_kutta.cxx sequential_decay.h sequential_decay.cxx bicycling.h bicycling.cxx hw1.cxx
#	g++ runge_kutta.h runge_kutta.cxx sequential_decay.h sequential_decay.cxx bicycling.h bicycling.cxx hw1.cxx -o main.o `root-config --cflags --glibs`
