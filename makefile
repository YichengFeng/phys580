main.o: nuclei_decay.h nuclei_decay.cxx lab1.cxx
	g++ nuclei_decay.h nuclei_decay.cxx lab1.cxx -o main.o `root-config --cflags --glibs`
