main.o: vector_field_2d.o loop_magnetic_field.o lab8.cxx
	g++ vector_field_2d.o loop_magnetic_field.o lab8.cxx -o main.o `root-config --cflags --glibs`

#test.o: vector_field_2d.o test_field_plot.cxx
#	g++ vector_field_2d.o test_field_plot.cxx -o test.o `root-config --cflags --glibs`

loop_magnetic_field.o: loop_magnetic_field.h loop_magnetic_field.cxx
	g++ -c loop_magnetic_field.h loop_magnetic_field.cxx

vector_field_2d.o: vector_field_2d.h vector_field_2d.cxx
	g++ -c vector_field_2d.h vector_field_2d.cxx `root-config --cflags --glibs`
