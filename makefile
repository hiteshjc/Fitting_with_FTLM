CC=g++ -O3 -fopenmp
CODE_SOURCES = ./src/*.C
SOURCES= $(CODE_SOURCES)
#LAPACK= -L /usr/lib -llapack -lblas -lm -lgomp -lgsl -lgslcblas -lpthread 
#LAPACK= -L /cm/shared/apps/lapack/gcc/5.2.0/ -L /usr/lib -L /cm/shared/apps/blas/gcc/5.2.0/lib64/ -llapack -lblas -lm -lgomp -fopenmp -ffast-math -lpthread
LAPACK= -L /cm/shared/apps/lapack/open64/64/3.5.0/ -L /usr/lib -L /cm/shared/apps/blas/gcc/current/lib64/ -llapack -lblas -lm -lgomp -fopenmp -ffast-math -lpthread
BOOST= -I /cm/shared/apps/boost/1.55.0/include/boost/ -I /cm/shared/apps/boost/1.55.0/include/
#GSL_INC= -I /sw/xe/gsl/1.16-2015-04/cnl5.2_gnu4.8.2/include 
#GSL_LIB= -L /sw/xe/gsl/1.16-2015-04/cnl5.2_gnu4.8.2/lib -lgsl  
EXECUTABLE=ed

CODE_OBJECTS=\
	./obj/ed.o	\
	./obj/hamiltonian_spin_functions.o	\
	./obj/main_prog.o \
	./obj/math_utilities.o \
	./obj/matrix_functions.o	\
	./obj/mtrand.o	\
	./obj/number_functions.o \
	./obj/printing_functions.o \
	./obj/search_for.o	\
	./obj/vector_utilities.o \
	./obj/ross.o \
	
	BIT_FORMAT_STR=64	

OBJECTS= $(CODE_OBJECTS)

headers1=./src/*.h 

$(EXECUTABLE): $(OBJECTS) $(headers1)
	$(CC) $(GSL_INC) $(NO_WRITE) $(OBJECTS) $(CFLAGS) $(LAPACK) $(GSL_LIB) -o $@

$(CODE_OBJECTS) : ./obj/%.o : ./src/%.C
	@echo "compiling $<";
	@$(CC) $(GSL_INC) $(BLAS_INC) $(BOOST) -c $< -o 	$@;

clean:
	rm -f $(OBJECTS) $(PROG)
	rm -f $(EXECUTABLE)

depend:
	makedepend $(INCLUDE) $(GSL_INC) $(BOOST) $(BLAS_INC) $(GSL_LIB) $(BLAS_LIB) -- -o $(CFLAGS) -- $(SOURCES) 



# DO NOT DELETE

