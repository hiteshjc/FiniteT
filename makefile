CC=g++ -O3 -fopenmp
CODE_SOURCES = ./src/*.C
#NAUTY_SOURCES = ./nauty/nauty.c ./nauty/nautil.c ./nauty/nausparse.c
SOURCES= $(CODE_SOURCES) $(NAUTY_SOURCES)
#NO_WRITE= -Wno-write-strings -Wno-deprecated -fno-stack-protector
#CFLAGS= -O2 -D
#LAPACK= -L /usr/lib64/ -llapack -lblas -lm -lpthread -lgomp 
#LAPACK= -llapack -lblas -lm -lpthread -lgomp 
LAPACK= -L /opt/lapack/gnu/lib/ -L /opt/modulefiles/applications/.gnu/lapack/.version.3.6.0 -L /cm/shared/apps/lapack/open64/64/3.5.0/ -L /usr/lib -L /cm/shared/apps/blas/gcc/current/lib64/ -llapack -lblas -lm -lgomp -fopenmp -ffast-math -lpthread
BOOST= -I /cm/shared/apps/boost/1.55.0/include/boost/ -I /cm/shared/apps/boost/1.55.0/include/

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
	./obj/spin_functions.o	\
	./obj/spin_model.o	\
	./obj/qddm.o	\
	./obj/vector_utilities.o \

#NAUTY_OBJECTS=\
#	./obj/nauty.o	\
#	./obj/nautil.o	\
#	./obj/nausparse.o	\
	
	BIT_FORMAT_STR=64	

OBJECTS= $(CODE_OBJECTS) $(NAUTY_OBJECTS)

headers1=./src/*.h 
#headers2=./nauty/*.h

$(EXECUTABLE): $(OBJECTS) $(headers1) $(headers2)
	$(CC) $(NO_WRITE) $(OBJECTS) $(CFLAGS) $(LAPACK) -o $@

$(CODE_OBJECTS) : ./obj/%.o : ./src/%.C
	@echo "compiling $<";
	@$(CC) $(BLAS_INC) $(BOOST) -c $< -o 	$@;

#$(NAUTY_OBJECTS) : ./obj/%.o : ./nauty/%.c
#	@echo "compiling $<";
#	@$(CC) $(BLAS_INC) $(BOOST) -c $< -o $@;

clean:
	rm -f $(OBJECTS) $(PROG)
	rm -f $(EXECUTABLE)

depend:
	makedepend $(INCLUDE) $(BOOST) $(BLAS_INC) $(BLAS_LIB) -- -o $(CFLAGS) -- $(SOURCES) 


# DO NOT DELETE

./src/bethe_automorphism.o: ./src/bethe_automorphism.h ./src/global.h
./src/bethe_automorphism.o: /usr/include/stdio.h /usr/include/features.h
./src/bethe_automorphism.o: /usr/include/sys/cdefs.h
./src/bethe_automorphism.o: /usr/include/bits/wordsize.h
./src/bethe_automorphism.o: /usr/include/gnu/stubs.h
./src/bethe_automorphism.o: /usr/include/gnu/stubs-64.h
./src/bethe_automorphism.o: /usr/include/bits/types.h
./src/bethe_automorphism.o: /usr/include/bits/typesizes.h
./src/bethe_automorphism.o: /usr/include/libio.h /usr/include/_G_config.h
./src/bethe_automorphism.o: /usr/include/wchar.h
./src/bethe_automorphism.o: /usr/include/bits/stdio_lim.h
./src/bethe_automorphism.o: /usr/include/bits/sys_errlist.h ./src/matrix.h
./src/bethe_automorphism.o: ./nauty/nausparse.h ./nauty/nauty.h
./src/bethe_automorphism.o: /usr/include/sys/types.h /usr/include/time.h
./src/bethe_automorphism.o: /usr/include/endian.h /usr/include/bits/endian.h
./src/bethe_automorphism.o: /usr/include/bits/byteswap.h
./src/bethe_automorphism.o: /usr/include/sys/select.h
./src/bethe_automorphism.o: /usr/include/bits/select.h
./src/bethe_automorphism.o: /usr/include/bits/sigset.h
./src/bethe_automorphism.o: /usr/include/bits/time.h
./src/bethe_automorphism.o: /usr/include/sys/sysmacros.h
./src/bethe_automorphism.o: /usr/include/bits/pthreadtypes.h
./src/bethe_automorphism.o: /usr/include/unistd.h
./src/bethe_automorphism.o: /usr/include/bits/posix_opt.h
./src/bethe_automorphism.o: /usr/include/bits/confname.h
./src/bethe_automorphism.o: /usr/include/getopt.h /usr/include/stdlib.h
./src/bethe_automorphism.o: /usr/include/alloca.h /usr/include/string.h
./src/bethe_automorphism.o: /usr/include/xlocale.h
./src/bethe_isomorphism.o: ./src/bethe_isomorphism.h ./src/global.h
./src/bethe_isomorphism.o: /usr/include/stdio.h /usr/include/features.h
./src/bethe_isomorphism.o: /usr/include/sys/cdefs.h
./src/bethe_isomorphism.o: /usr/include/bits/wordsize.h
./src/bethe_isomorphism.o: /usr/include/gnu/stubs.h
./src/bethe_isomorphism.o: /usr/include/gnu/stubs-64.h
./src/bethe_isomorphism.o: /usr/include/bits/types.h
./src/bethe_isomorphism.o: /usr/include/bits/typesizes.h /usr/include/libio.h
./src/bethe_isomorphism.o: /usr/include/_G_config.h /usr/include/wchar.h
./src/bethe_isomorphism.o: /usr/include/bits/stdio_lim.h
./src/bethe_isomorphism.o: /usr/include/bits/sys_errlist.h ./src/matrix.h
./src/bethe_isomorphism.o: ./nauty/nausparse.h ./nauty/nauty.h
./src/bethe_isomorphism.o: /usr/include/sys/types.h /usr/include/time.h
./src/bethe_isomorphism.o: /usr/include/endian.h /usr/include/bits/endian.h
./src/bethe_isomorphism.o: /usr/include/bits/byteswap.h
./src/bethe_isomorphism.o: /usr/include/sys/select.h
./src/bethe_isomorphism.o: /usr/include/bits/select.h
./src/bethe_isomorphism.o: /usr/include/bits/sigset.h
./src/bethe_isomorphism.o: /usr/include/bits/time.h
./src/bethe_isomorphism.o: /usr/include/sys/sysmacros.h
./src/bethe_isomorphism.o: /usr/include/bits/pthreadtypes.h
./src/bethe_isomorphism.o: /usr/include/unistd.h
./src/bethe_isomorphism.o: /usr/include/bits/posix_opt.h
./src/bethe_isomorphism.o: /usr/include/bits/confname.h /usr/include/getopt.h
./src/bethe_isomorphism.o: /usr/include/stdlib.h /usr/include/alloca.h
./src/bethe_isomorphism.o: /usr/include/string.h /usr/include/xlocale.h
./src/cayley.o: ./src/number_functions.h ./src/global.h /usr/include/stdio.h
./src/cayley.o: /usr/include/features.h /usr/include/sys/cdefs.h
./src/cayley.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
./src/cayley.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
./src/cayley.o: /usr/include/bits/typesizes.h /usr/include/libio.h
./src/cayley.o: /usr/include/_G_config.h /usr/include/wchar.h
./src/cayley.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
./src/cayley.o: ./src/matrix.h
./src/chris_fm_phi.o: ./src/chris_fm_phi.h ./src/global.h
./src/chris_fm_phi.o: /usr/include/stdio.h /usr/include/features.h
./src/chris_fm_phi.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
./src/chris_fm_phi.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
./src/chris_fm_phi.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
./src/chris_fm_phi.o: /usr/include/libio.h /usr/include/_G_config.h
./src/chris_fm_phi.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
./src/chris_fm_phi.o: /usr/include/bits/sys_errlist.h ./src/matrix.h
./src/chris_fm_phi.o: ./src/cluster.h ./src/number_functions.h
./src/chris_fm_phi.o: ./src/math_utilities.h ./src/vector_utilities.h
./src/chris_fm_phi.o: ./src/bethe_isomorphism.h ./nauty/nausparse.h
./src/chris_fm_phi.o: ./nauty/nauty.h /usr/include/sys/types.h
./src/chris_fm_phi.o: /usr/include/time.h /usr/include/endian.h
./src/chris_fm_phi.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
./src/chris_fm_phi.o: /usr/include/sys/select.h /usr/include/bits/select.h
./src/chris_fm_phi.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
./src/chris_fm_phi.o: /usr/include/sys/sysmacros.h
./src/chris_fm_phi.o: /usr/include/bits/pthreadtypes.h /usr/include/unistd.h
./src/chris_fm_phi.o: /usr/include/bits/posix_opt.h
./src/chris_fm_phi.o: /usr/include/bits/confname.h /usr/include/getopt.h
./src/chris_fm_phi.o: /usr/include/stdlib.h /usr/include/alloca.h
./src/chris_fm_phi.o: /usr/include/string.h /usr/include/xlocale.h
./src/chris_fm_phi.o: ./src/printing_functions.h
./src/cluster.o: ./src/number_functions.h ./src/global.h /usr/include/stdio.h
./src/cluster.o: /usr/include/features.h /usr/include/sys/cdefs.h
./src/cluster.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
./src/cluster.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
./src/cluster.o: /usr/include/bits/typesizes.h /usr/include/libio.h
./src/cluster.o: /usr/include/_G_config.h /usr/include/wchar.h
./src/cluster.o: /usr/include/bits/stdio_lim.h
./src/cluster.o: /usr/include/bits/sys_errlist.h ./src/matrix.h
./src/cluster.o: ./src/cluster.h ./src/math_utilities.h
./src/cluster.o: ./src/vector_utilities.h ./src/bethe_isomorphism.h
./src/cluster.o: ./nauty/nausparse.h ./nauty/nauty.h /usr/include/sys/types.h
./src/cluster.o: /usr/include/time.h /usr/include/endian.h
./src/cluster.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
./src/cluster.o: /usr/include/sys/select.h /usr/include/bits/select.h
./src/cluster.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
./src/cluster.o: /usr/include/sys/sysmacros.h
./src/cluster.o: /usr/include/bits/pthreadtypes.h /usr/include/unistd.h
./src/cluster.o: /usr/include/bits/posix_opt.h /usr/include/bits/confname.h
./src/cluster.o: /usr/include/getopt.h /usr/include/stdlib.h
./src/cluster.o: /usr/include/alloca.h /usr/include/string.h
./src/cluster.o: /usr/include/xlocale.h ./src/printing_functions.h
./src/cluster.o: ./src/search_for.h ./src/chris_fm_phi.h
./src/data_io.o: ./src/data_io.h ./src/global.h /usr/include/stdio.h
./src/data_io.o: /usr/include/features.h /usr/include/sys/cdefs.h
./src/data_io.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
./src/data_io.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
./src/data_io.o: /usr/include/bits/typesizes.h /usr/include/libio.h
./src/data_io.o: /usr/include/_G_config.h /usr/include/wchar.h
./src/data_io.o: /usr/include/bits/stdio_lim.h
./src/data_io.o: /usr/include/bits/sys_errlist.h ./src/matrix.h
./src/data_io.o: ./src/number_functions.h
./src/density_matrix.o: ./src/density_matrix.h ./src/global.h
./src/density_matrix.o: /usr/include/stdio.h /usr/include/features.h
./src/density_matrix.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
./src/density_matrix.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
./src/density_matrix.o: /usr/include/bits/types.h
./src/density_matrix.o: /usr/include/bits/typesizes.h /usr/include/libio.h
./src/density_matrix.o: /usr/include/_G_config.h /usr/include/wchar.h
./src/density_matrix.o: /usr/include/bits/stdio_lim.h
./src/density_matrix.o: /usr/include/bits/sys_errlist.h ./src/matrix.h
./src/density_matrix.o: ./src/matrix_functions.h ./src/printing_functions.h
./src/density_matrix.o: ./src/number_functions.h
./src/dmrg.o: ./src/dmrg.h ./src/cluster.h ./src/global.h
./src/dmrg.o: /usr/include/stdio.h /usr/include/features.h
./src/dmrg.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
./src/dmrg.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
./src/dmrg.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
./src/dmrg.o: /usr/include/libio.h /usr/include/_G_config.h
./src/dmrg.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
./src/dmrg.o: /usr/include/bits/sys_errlist.h ./src/matrix.h
./src/dmrg.o: ./src/number_functions.h ./src/math_utilities.h
./src/dmrg.o: ./src/vector_utilities.h ./src/bethe_isomorphism.h
./src/dmrg.o: ./nauty/nausparse.h ./nauty/nauty.h /usr/include/sys/types.h
./src/dmrg.o: /usr/include/time.h /usr/include/endian.h
./src/dmrg.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
./src/dmrg.o: /usr/include/sys/select.h /usr/include/bits/select.h
./src/dmrg.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
./src/dmrg.o: /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h
./src/dmrg.o: /usr/include/unistd.h /usr/include/bits/posix_opt.h
./src/dmrg.o: /usr/include/bits/confname.h /usr/include/getopt.h
./src/dmrg.o: /usr/include/stdlib.h /usr/include/alloca.h
./src/dmrg.o: /usr/include/string.h /usr/include/xlocale.h
./src/dmrg.o: ./src/printing_functions.h ./src/nrg_refined.h
./src/dmrg.o: ./src/matrix_functions.h ./src/density_matrix.h
./src/dmrg_spin.o: ./src/dmrg_spin.h ./src/cluster.h ./src/global.h
./src/dmrg_spin.o: /usr/include/stdio.h /usr/include/features.h
./src/dmrg_spin.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
./src/dmrg_spin.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
./src/dmrg_spin.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
./src/dmrg_spin.o: /usr/include/libio.h /usr/include/_G_config.h
./src/dmrg_spin.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
./src/dmrg_spin.o: /usr/include/bits/sys_errlist.h ./src/matrix.h
./src/dmrg_spin.o: ./src/number_functions.h ./src/math_utilities.h
./src/dmrg_spin.o: ./src/vector_utilities.h ./src/bethe_isomorphism.h
./src/dmrg_spin.o: ./nauty/nausparse.h ./nauty/nauty.h
./src/dmrg_spin.o: /usr/include/sys/types.h /usr/include/time.h
./src/dmrg_spin.o: /usr/include/endian.h /usr/include/bits/endian.h
./src/dmrg_spin.o: /usr/include/bits/byteswap.h /usr/include/sys/select.h
./src/dmrg_spin.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
./src/dmrg_spin.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
./src/dmrg_spin.o: /usr/include/bits/pthreadtypes.h /usr/include/unistd.h
./src/dmrg_spin.o: /usr/include/bits/posix_opt.h /usr/include/bits/confname.h
./src/dmrg_spin.o: /usr/include/getopt.h /usr/include/stdlib.h
./src/dmrg_spin.o: /usr/include/alloca.h /usr/include/string.h
./src/dmrg_spin.o: /usr/include/xlocale.h ./src/printing_functions.h
./src/dmrg_spin.o: ./src/nrg_spin.h ./src/matrix_functions.h
./src/dmrg_spin.o: ./src/density_matrix.h ./src/bethe_lapack_interface.h
./src/ed.o: ./src/ed.h ./src/hamiltonian.h ./src/global.h
./src/ed.o: /usr/include/stdio.h /usr/include/features.h
./src/ed.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
./src/ed.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
./src/ed.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
./src/ed.o: /usr/include/libio.h /usr/include/_G_config.h
./src/ed.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
./src/ed.o: /usr/include/bits/sys_errlist.h ./src/matrix.h
./src/ed.o: ./src/matrix_functions.h ./src/number_functions.h
./src/ed.o: ./src/printing_functions.h ./src/bethe_lapack_interface.h
./src/ed.o: ./src/math_utilities.h ./src/hamiltonian_spin_functions.h
./src/ed.o: ./src/density_matrix.h
./src/eff_ham.o: ./src/eff_ham.h ./src/hamiltonian.h ./src/global.h
./src/eff_ham.o: /usr/include/stdio.h /usr/include/features.h
./src/eff_ham.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
./src/eff_ham.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
./src/eff_ham.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
./src/eff_ham.o: /usr/include/libio.h /usr/include/_G_config.h
./src/eff_ham.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
./src/eff_ham.o: /usr/include/bits/sys_errlist.h ./src/matrix.h
./src/eff_ham.o: ./src/matrix_functions.h ./src/number_functions.h
./src/eff_ham.o: ./src/printing_functions.h ./src/bethe_lapack_interface.h
./src/eff_ham.o: ./src/math_utilities.h ./src/hamiltonian_spin_functions.h
./src/eff_ham.o: ./src/lanczos_more.h ./src/ed.h ./src/density_matrix.h
./src/exponent.o: ./src/hamiltonian.h ./src/global.h /usr/include/stdio.h
./src/exponent.o: /usr/include/features.h /usr/include/sys/cdefs.h
./src/exponent.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
./src/exponent.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
./src/exponent.o: /usr/include/bits/typesizes.h /usr/include/libio.h
./src/exponent.o: /usr/include/_G_config.h /usr/include/wchar.h
./src/exponent.o: /usr/include/bits/stdio_lim.h
./src/exponent.o: /usr/include/bits/sys_errlist.h ./src/matrix.h ./src/ed.h
./src/exponent.o: ./src/matrix_functions.h ./src/number_functions.h
./src/exponent.o: ./src/printing_functions.h ./src/bethe_lapack_interface.h
./src/exponent.o: ./src/math_utilities.h ./src/hamiltonian_spin_functions.h
./src/exponent.o: ./src/density_matrix.h ./src/xxz_model.h ./src/search_for.h
./src/exponent.o: ./src/dmrg_spin.h ./src/cluster.h ./src/vector_utilities.h
./src/exponent.o: ./src/bethe_isomorphism.h ./nauty/nausparse.h
./src/exponent.o: ./nauty/nauty.h /usr/include/sys/types.h
./src/exponent.o: /usr/include/time.h /usr/include/endian.h
./src/exponent.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
./src/exponent.o: /usr/include/sys/select.h /usr/include/bits/select.h
./src/exponent.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
./src/exponent.o: /usr/include/sys/sysmacros.h
./src/exponent.o: /usr/include/bits/pthreadtypes.h /usr/include/unistd.h
./src/exponent.o: /usr/include/bits/posix_opt.h /usr/include/bits/confname.h
./src/exponent.o: /usr/include/getopt.h /usr/include/stdlib.h
./src/exponent.o: /usr/include/alloca.h /usr/include/string.h
./src/exponent.o: /usr/include/xlocale.h ./src/data_io.h
./src/exponent_old.o: ./src/exponent.h ./src/cluster.h ./src/global.h
./src/exponent_old.o: /usr/include/stdio.h /usr/include/features.h
./src/exponent_old.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
./src/exponent_old.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
./src/exponent_old.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
./src/exponent_old.o: /usr/include/libio.h /usr/include/_G_config.h
./src/exponent_old.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
./src/exponent_old.o: /usr/include/bits/sys_errlist.h ./src/matrix.h
./src/exponent_old.o: ./src/number_functions.h ./src/math_utilities.h
./src/exponent_old.o: ./src/vector_utilities.h ./src/bethe_isomorphism.h
./src/exponent_old.o: ./nauty/nausparse.h ./nauty/nauty.h
./src/exponent_old.o: /usr/include/sys/types.h /usr/include/time.h
./src/exponent_old.o: /usr/include/endian.h /usr/include/bits/endian.h
./src/exponent_old.o: /usr/include/bits/byteswap.h /usr/include/sys/select.h
./src/exponent_old.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
./src/exponent_old.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
./src/exponent_old.o: /usr/include/bits/pthreadtypes.h /usr/include/unistd.h
./src/exponent_old.o: /usr/include/bits/posix_opt.h
./src/exponent_old.o: /usr/include/bits/confname.h /usr/include/getopt.h
./src/exponent_old.o: /usr/include/stdlib.h /usr/include/alloca.h
./src/exponent_old.o: /usr/include/string.h /usr/include/xlocale.h
./src/exponent_old.o: ./src/printing_functions.h ./src/nrg.h
./src/exponent_old.o: ./src/matrix_functions.h ./src/density_matrix.h
./src/exponent_old.o: ./src/nrg_refined.h ./src/data_io.h ./src/hamiltonian.h
./src/exponent_old.o: ./src/ed.h ./src/bethe_lapack_interface.h
./src/exponent_old.o: ./src/hamiltonian_spin_functions.h ./src/xxz_model.h
./src/four_dangling_analysis.o: ./src/four_dangling_analysis.h
./src/four_dangling_analysis.o: ./src/lanczos_more.h ./src/hamiltonian.h
./src/four_dangling_analysis.o: ./src/global.h /usr/include/stdio.h
./src/four_dangling_analysis.o: /usr/include/features.h
./src/four_dangling_analysis.o: /usr/include/sys/cdefs.h
./src/four_dangling_analysis.o: /usr/include/bits/wordsize.h
./src/four_dangling_analysis.o: /usr/include/gnu/stubs.h
./src/four_dangling_analysis.o: /usr/include/gnu/stubs-64.h
./src/four_dangling_analysis.o: /usr/include/bits/types.h
./src/four_dangling_analysis.o: /usr/include/bits/typesizes.h
./src/four_dangling_analysis.o: /usr/include/libio.h /usr/include/_G_config.h
./src/four_dangling_analysis.o: /usr/include/wchar.h
./src/four_dangling_analysis.o: /usr/include/bits/stdio_lim.h
./src/four_dangling_analysis.o: /usr/include/bits/sys_errlist.h
./src/four_dangling_analysis.o: ./src/matrix.h ./src/matrix_functions.h
./src/four_dangling_analysis.o: ./src/number_functions.h
./src/four_dangling_analysis.o: ./src/printing_functions.h
./src/four_dangling_analysis.o: ./src/bethe_lapack_interface.h
./src/four_dangling_analysis.o: ./src/math_utilities.h
./src/four_dangling_analysis.o: ./src/hamiltonian_spin_functions.h ./src/ed.h
./src/four_dangling_analysis.o: ./src/density_matrix.h
./src/generic_dmrg_spin.o: ./src/generic_dmrg_spin.h ./src/global.h
./src/generic_dmrg_spin.o: /usr/include/stdio.h /usr/include/features.h
./src/generic_dmrg_spin.o: /usr/include/sys/cdefs.h
./src/generic_dmrg_spin.o: /usr/include/bits/wordsize.h
./src/generic_dmrg_spin.o: /usr/include/gnu/stubs.h
./src/generic_dmrg_spin.o: /usr/include/gnu/stubs-64.h
./src/generic_dmrg_spin.o: /usr/include/bits/types.h
./src/generic_dmrg_spin.o: /usr/include/bits/typesizes.h /usr/include/libio.h
./src/generic_dmrg_spin.o: /usr/include/_G_config.h /usr/include/wchar.h
./src/generic_dmrg_spin.o: /usr/include/bits/stdio_lim.h
./src/generic_dmrg_spin.o: /usr/include/bits/sys_errlist.h ./src/matrix.h
./src/generic_dmrg_spin.o: ./src/cluster.h ./src/number_functions.h
./src/generic_dmrg_spin.o: ./src/math_utilities.h ./src/vector_utilities.h
./src/generic_dmrg_spin.o: ./src/bethe_isomorphism.h ./nauty/nausparse.h
./src/generic_dmrg_spin.o: ./nauty/nauty.h /usr/include/sys/types.h
./src/generic_dmrg_spin.o: /usr/include/time.h /usr/include/endian.h
./src/generic_dmrg_spin.o: /usr/include/bits/endian.h
./src/generic_dmrg_spin.o: /usr/include/bits/byteswap.h
./src/generic_dmrg_spin.o: /usr/include/sys/select.h
./src/generic_dmrg_spin.o: /usr/include/bits/select.h
./src/generic_dmrg_spin.o: /usr/include/bits/sigset.h
./src/generic_dmrg_spin.o: /usr/include/bits/time.h
./src/generic_dmrg_spin.o: /usr/include/sys/sysmacros.h
./src/generic_dmrg_spin.o: /usr/include/bits/pthreadtypes.h
./src/generic_dmrg_spin.o: /usr/include/unistd.h
./src/generic_dmrg_spin.o: /usr/include/bits/posix_opt.h
./src/generic_dmrg_spin.o: /usr/include/bits/confname.h /usr/include/getopt.h
./src/generic_dmrg_spin.o: /usr/include/stdlib.h /usr/include/alloca.h
./src/generic_dmrg_spin.o: /usr/include/string.h /usr/include/xlocale.h
./src/generic_dmrg_spin.o: ./src/printing_functions.h
./src/generic_dmrg_spin.o: ./src/matrix_functions.h ./src/spin_model.h
./src/generic_dmrg_spin.o: ./src/hamiltonian.h ./src/block.h
./src/generic_dmrg_spin.o: ./src/spin_functions.h ./src/tmp_info.h
./src/generic_dmrg_spin.o: ./src/universe.h ./src/generic_nrg_spin.h
./src/generic_dmrg_spin.o: ./src/density_matrix.h
./src/generic_dmrg_spin.o: ./src/bethe_lapack_interface.h
./src/generic_nrg_spin.o: ./src/generic_nrg_spin.h ./src/global.h
./src/generic_nrg_spin.o: /usr/include/stdio.h /usr/include/features.h
./src/generic_nrg_spin.o: /usr/include/sys/cdefs.h
./src/generic_nrg_spin.o: /usr/include/bits/wordsize.h
./src/generic_nrg_spin.o: /usr/include/gnu/stubs.h
./src/generic_nrg_spin.o: /usr/include/gnu/stubs-64.h
./src/generic_nrg_spin.o: /usr/include/bits/types.h
./src/generic_nrg_spin.o: /usr/include/bits/typesizes.h /usr/include/libio.h
./src/generic_nrg_spin.o: /usr/include/_G_config.h /usr/include/wchar.h
./src/generic_nrg_spin.o: /usr/include/bits/stdio_lim.h
./src/generic_nrg_spin.o: /usr/include/bits/sys_errlist.h ./src/matrix.h
./src/generic_nrg_spin.o: ./src/cluster.h ./src/number_functions.h
./src/generic_nrg_spin.o: ./src/math_utilities.h ./src/vector_utilities.h
./src/generic_nrg_spin.o: ./src/bethe_isomorphism.h ./nauty/nausparse.h
./src/generic_nrg_spin.o: ./nauty/nauty.h /usr/include/sys/types.h
./src/generic_nrg_spin.o: /usr/include/time.h /usr/include/endian.h
./src/generic_nrg_spin.o: /usr/include/bits/endian.h
./src/generic_nrg_spin.o: /usr/include/bits/byteswap.h
./src/generic_nrg_spin.o: /usr/include/sys/select.h
./src/generic_nrg_spin.o: /usr/include/bits/select.h
./src/generic_nrg_spin.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
./src/generic_nrg_spin.o: /usr/include/sys/sysmacros.h
./src/generic_nrg_spin.o: /usr/include/bits/pthreadtypes.h
./src/generic_nrg_spin.o: /usr/include/unistd.h /usr/include/bits/posix_opt.h
./src/generic_nrg_spin.o: /usr/include/bits/confname.h /usr/include/getopt.h
./src/generic_nrg_spin.o: /usr/include/stdlib.h /usr/include/alloca.h
./src/generic_nrg_spin.o: /usr/include/string.h /usr/include/xlocale.h
./src/generic_nrg_spin.o: ./src/printing_functions.h ./src/matrix_functions.h
./src/generic_nrg_spin.o: ./src/spin_model.h ./src/hamiltonian.h
./src/generic_nrg_spin.o: ./src/block.h ./src/spin_functions.h
./src/generic_nrg_spin.o: ./src/tmp_info.h ./src/universe.h ./src/ed.h
./src/generic_nrg_spin.o: ./src/bethe_lapack_interface.h
./src/generic_nrg_spin.o: ./src/hamiltonian_spin_functions.h
./src/generic_nrg_spin.o: ./src/density_matrix.h
./src/get_ai.o: ./src/get_ai.h ./src/bethe_lapack_interface.h ./src/global.h
./src/get_ai.o: /usr/include/stdio.h /usr/include/features.h
./src/get_ai.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
./src/get_ai.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
./src/get_ai.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
./src/get_ai.o: /usr/include/libio.h /usr/include/_G_config.h
./src/get_ai.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
./src/get_ai.o: /usr/include/bits/sys_errlist.h ./src/matrix.h
./src/get_ai.o: ./src/printing_functions.h ./src/matrix_functions.h
./src/hamiltonian_spin_functions.o: ./src/hamiltonian_spin_functions.h
./src/hamiltonian_spin_functions.o: ./src/global.h /usr/include/stdio.h
./src/hamiltonian_spin_functions.o: /usr/include/features.h
./src/hamiltonian_spin_functions.o: /usr/include/sys/cdefs.h
./src/hamiltonian_spin_functions.o: /usr/include/bits/wordsize.h
./src/hamiltonian_spin_functions.o: /usr/include/gnu/stubs.h
./src/hamiltonian_spin_functions.o: /usr/include/gnu/stubs-64.h
./src/hamiltonian_spin_functions.o: /usr/include/bits/types.h
./src/hamiltonian_spin_functions.o: /usr/include/bits/typesizes.h
./src/hamiltonian_spin_functions.o: /usr/include/libio.h
./src/hamiltonian_spin_functions.o: /usr/include/_G_config.h
./src/hamiltonian_spin_functions.o: /usr/include/wchar.h
./src/hamiltonian_spin_functions.o: /usr/include/bits/stdio_lim.h
./src/hamiltonian_spin_functions.o: /usr/include/bits/sys_errlist.h
./src/hamiltonian_spin_functions.o: ./src/matrix.h ./src/number_functions.h
./src/hamiltonian_spin_functions.o: ./src/spin_functions.h
./src/henley_1d.o: ./src/global.h /usr/include/stdio.h
./src/henley_1d.o: /usr/include/features.h /usr/include/sys/cdefs.h
./src/henley_1d.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
./src/henley_1d.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
./src/henley_1d.o: /usr/include/bits/typesizes.h /usr/include/libio.h
./src/henley_1d.o: /usr/include/_G_config.h /usr/include/wchar.h
./src/henley_1d.o: /usr/include/bits/stdio_lim.h
./src/henley_1d.o: /usr/include/bits/sys_errlist.h ./src/matrix.h
./src/henley_1d.o: ./src/henley_1d.h ./src/hamiltonian.h
./src/henley_1d.o: ./src/hamiltonian_spin_functions.h
./src/henley_1d.o: ./src/number_functions.h ./src/search_for.h
./src/henley_1d.o: ./src/printing_functions.h
./src/hida.o: ./src/hida.h ./src/global.h /usr/include/stdio.h
./src/hida.o: /usr/include/features.h /usr/include/sys/cdefs.h
./src/hida.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
./src/hida.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
./src/hida.o: /usr/include/bits/typesizes.h /usr/include/libio.h
./src/hida.o: /usr/include/_G_config.h /usr/include/wchar.h
./src/hida.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
./src/hida.o: ./src/matrix.h ./src/printing_functions.h
./src/hida.o: ./src/bethe_lapack_interface.h ./src/math_utilities.h
./src/hida.o: ./src/spin_model.h ./src/hamiltonian.h ./src/ed.h
./src/hida.o: ./src/matrix_functions.h ./src/number_functions.h
./src/hida.o: ./src/hamiltonian_spin_functions.h ./src/density_matrix.h
./src/lanczos_more.bk.Feb02.o: ./src/lanczos_more.h ./src/hamiltonian.h
./src/lanczos_more.bk.Feb02.o: ./src/global.h /usr/include/stdio.h
./src/lanczos_more.bk.Feb02.o: /usr/include/features.h
./src/lanczos_more.bk.Feb02.o: /usr/include/sys/cdefs.h
./src/lanczos_more.bk.Feb02.o: /usr/include/bits/wordsize.h
./src/lanczos_more.bk.Feb02.o: /usr/include/gnu/stubs.h
./src/lanczos_more.bk.Feb02.o: /usr/include/gnu/stubs-64.h
./src/lanczos_more.bk.Feb02.o: /usr/include/bits/types.h
./src/lanczos_more.bk.Feb02.o: /usr/include/bits/typesizes.h
./src/lanczos_more.bk.Feb02.o: /usr/include/libio.h /usr/include/_G_config.h
./src/lanczos_more.bk.Feb02.o: /usr/include/wchar.h
./src/lanczos_more.bk.Feb02.o: /usr/include/bits/stdio_lim.h
./src/lanczos_more.bk.Feb02.o: /usr/include/bits/sys_errlist.h ./src/matrix.h
./src/lanczos_more.bk.Feb02.o: ./src/matrix_functions.h
./src/lanczos_more.bk.Feb02.o: ./src/number_functions.h
./src/lanczos_more.bk.Feb02.o: ./src/printing_functions.h
./src/lanczos_more.bk.Feb02.o: ./src/bethe_lapack_interface.h
./src/lanczos_more.bk.Feb02.o: ./src/math_utilities.h
./src/lanczos_more.bk.Feb02.o: ./src/hamiltonian_spin_functions.h ./src/ed.h
./src/lanczos_more.bk.Feb02.o: ./src/density_matrix.h
./src/lanczos_more.o: ./src/lanczos_more.h ./src/hamiltonian.h ./src/global.h
./src/lanczos_more.o: /usr/include/stdio.h /usr/include/features.h
./src/lanczos_more.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
./src/lanczos_more.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
./src/lanczos_more.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
./src/lanczos_more.o: /usr/include/libio.h /usr/include/_G_config.h
./src/lanczos_more.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
./src/lanczos_more.o: /usr/include/bits/sys_errlist.h ./src/matrix.h
./src/lanczos_more.o: ./src/matrix_functions.h ./src/number_functions.h
./src/lanczos_more.o: ./src/printing_functions.h
./src/lanczos_more.o: ./src/bethe_lapack_interface.h ./src/math_utilities.h
./src/lanczos_more.o: ./src/hamiltonian_spin_functions.h ./src/ed.h
./src/lanczos_more.o: ./src/density_matrix.h
./src/laplacian.o: ./src/laplacian.h ./src/global.h /usr/include/stdio.h
./src/laplacian.o: /usr/include/features.h /usr/include/sys/cdefs.h
./src/laplacian.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
./src/laplacian.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
./src/laplacian.o: /usr/include/bits/typesizes.h /usr/include/libio.h
./src/laplacian.o: /usr/include/_G_config.h /usr/include/wchar.h
./src/laplacian.o: /usr/include/bits/stdio_lim.h
./src/laplacian.o: /usr/include/bits/sys_errlist.h ./src/matrix.h
./src/laplacian.o: ./src/matrix_functions.h ./src/number_functions.h
./src/laplacian.o: ./src/bethe_lapack_interface.h ./src/printing_functions.h
./src/main_prog.o: ./src/global.h /usr/include/stdio.h
./src/main_prog.o: /usr/include/features.h /usr/include/sys/cdefs.h
./src/main_prog.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
./src/main_prog.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
./src/main_prog.o: /usr/include/bits/typesizes.h /usr/include/libio.h
./src/main_prog.o: /usr/include/_G_config.h /usr/include/wchar.h
./src/main_prog.o: /usr/include/bits/stdio_lim.h
./src/main_prog.o: /usr/include/bits/sys_errlist.h ./src/matrix.h
./src/main_prog.o: ./src/data_io.h ./src/block.h ./src/spin_functions.h
./src/main_prog.o: ./src/tmp_info.h ./src/tests.h ./src/number_functions.h
./src/main_prog.o: ./src/math_utilities.h ./src/vector_utilities.h
./src/main_prog.o: ./src/bethe_isomorphism.h ./nauty/nausparse.h
./src/main_prog.o: ./nauty/nauty.h /usr/include/sys/types.h
./src/main_prog.o: /usr/include/time.h /usr/include/endian.h
./src/main_prog.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
./src/main_prog.o: /usr/include/sys/select.h /usr/include/bits/select.h
./src/main_prog.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
./src/main_prog.o: /usr/include/sys/sysmacros.h
./src/main_prog.o: /usr/include/bits/pthreadtypes.h /usr/include/unistd.h
./src/main_prog.o: /usr/include/bits/posix_opt.h /usr/include/bits/confname.h
./src/main_prog.o: /usr/include/getopt.h /usr/include/stdlib.h
./src/main_prog.o: /usr/include/alloca.h /usr/include/string.h
./src/main_prog.o: /usr/include/xlocale.h ./src/printing_functions.h
./src/main_prog.o: ./src/cluster.h ./src/search_for.h ./src/mtrand.h
./src/main_prog.o: ./src/hamiltonian.h ./src/hida.h
./src/main_prog.o: ./src/bethe_lapack_interface.h ./src/tfi_model.h
./src/main_prog.o: ./src/xxz_model.h ./src/spin_model.h ./src/henley_1d.h
./src/main_prog.o: ./src/spinless_model.h ./src/spinless_boson_model.h
./src/main_prog.o: ./src/spin_waves.h ./src/matrix_functions.h
./src/main_prog.o: ./src/laplacian.h ./src/ed.h
./src/main_prog.o: ./src/hamiltonian_spin_functions.h ./src/density_matrix.h
./src/main_prog.o: ./src/lanczos_more.h ./src/nrg.h ./src/nrg_refined.h
./src/main_prog.o: ./src/nrg_spin.h ./src/generic_nrg_spin.h ./src/universe.h
./src/main_prog.o: ./src/tree_nrg_spin.h ./src/tensor.h
./src/main_prog.o: ./src/tree_dmrg_spin.h ./src/generic_dmrg_spin.h
./src/main_prog.o: ./src/dmrg.h ./src/dmrg_spin.h ./src/exponent.h
./src/main_prog.o: ./src/exponent_old.h ./src/singlet_triplet.h ./src/qddm.h
./src/main_prog.o: ./src/two_dangling_analysis.h
./src/main_prog.o: ./src/four_dangling_analysis.h ./src/sumi_ct.h
./src/main_prog.o: ./src/get_ai.h
./src/math_utilities.o: ./src/math_utilities.h ./src/global.h
./src/math_utilities.o: /usr/include/stdio.h /usr/include/features.h
./src/math_utilities.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
./src/math_utilities.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
./src/math_utilities.o: /usr/include/bits/types.h
./src/math_utilities.o: /usr/include/bits/typesizes.h /usr/include/libio.h
./src/math_utilities.o: /usr/include/_G_config.h /usr/include/wchar.h
./src/math_utilities.o: /usr/include/bits/stdio_lim.h
./src/math_utilities.o: /usr/include/bits/sys_errlist.h ./src/matrix.h
./src/math_utilities.o: ./src/mtrand.h ./src/number_functions.h
./src/math_utilities.o: ./src/search_for.h
./src/matrix_functions.o: ./src/bethe_lapack_interface.h
./src/matrix_functions.o: ./src/matrix_functions.h ./src/global.h
./src/matrix_functions.o: /usr/include/stdio.h /usr/include/features.h
./src/matrix_functions.o: /usr/include/sys/cdefs.h
./src/matrix_functions.o: /usr/include/bits/wordsize.h
./src/matrix_functions.o: /usr/include/gnu/stubs.h
./src/matrix_functions.o: /usr/include/gnu/stubs-64.h
./src/matrix_functions.o: /usr/include/bits/types.h
./src/matrix_functions.o: /usr/include/bits/typesizes.h /usr/include/libio.h
./src/matrix_functions.o: /usr/include/_G_config.h /usr/include/wchar.h
./src/matrix_functions.o: /usr/include/bits/stdio_lim.h
./src/matrix_functions.o: /usr/include/bits/sys_errlist.h ./src/matrix.h
./src/matrix_functions.o: ./src/printing_functions.h ./src/math_utilities.h
./src/mtrand.o: ./src/mtrand.h
./src/nrg.o: ./src/nrg.h ./src/global.h /usr/include/stdio.h
./src/nrg.o: /usr/include/features.h /usr/include/sys/cdefs.h
./src/nrg.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
./src/nrg.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
./src/nrg.o: /usr/include/bits/typesizes.h /usr/include/libio.h
./src/nrg.o: /usr/include/_G_config.h /usr/include/wchar.h
./src/nrg.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
./src/nrg.o: ./src/matrix.h ./src/cluster.h ./src/number_functions.h
./src/nrg.o: ./src/math_utilities.h ./src/vector_utilities.h
./src/nrg.o: ./src/bethe_isomorphism.h ./nauty/nausparse.h ./nauty/nauty.h
./src/nrg.o: /usr/include/sys/types.h /usr/include/time.h
./src/nrg.o: /usr/include/endian.h /usr/include/bits/endian.h
./src/nrg.o: /usr/include/bits/byteswap.h /usr/include/sys/select.h
./src/nrg.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
./src/nrg.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
./src/nrg.o: /usr/include/bits/pthreadtypes.h /usr/include/unistd.h
./src/nrg.o: /usr/include/bits/posix_opt.h /usr/include/bits/confname.h
./src/nrg.o: /usr/include/getopt.h /usr/include/stdlib.h
./src/nrg.o: /usr/include/alloca.h /usr/include/string.h
./src/nrg.o: /usr/include/xlocale.h ./src/printing_functions.h
./src/nrg.o: ./src/matrix_functions.h ./src/density_matrix.h ./src/ed.h
./src/nrg.o: ./src/hamiltonian.h ./src/bethe_lapack_interface.h
./src/nrg.o: ./src/hamiltonian_spin_functions.h
./src/nrg_refined.o: ./src/nrg_refined.h ./src/global.h /usr/include/stdio.h
./src/nrg_refined.o: /usr/include/features.h /usr/include/sys/cdefs.h
./src/nrg_refined.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
./src/nrg_refined.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
./src/nrg_refined.o: /usr/include/bits/typesizes.h /usr/include/libio.h
./src/nrg_refined.o: /usr/include/_G_config.h /usr/include/wchar.h
./src/nrg_refined.o: /usr/include/bits/stdio_lim.h
./src/nrg_refined.o: /usr/include/bits/sys_errlist.h ./src/matrix.h
./src/nrg_refined.o: ./src/cluster.h ./src/number_functions.h
./src/nrg_refined.o: ./src/math_utilities.h ./src/vector_utilities.h
./src/nrg_refined.o: ./src/bethe_isomorphism.h ./nauty/nausparse.h
./src/nrg_refined.o: ./nauty/nauty.h /usr/include/sys/types.h
./src/nrg_refined.o: /usr/include/time.h /usr/include/endian.h
./src/nrg_refined.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
./src/nrg_refined.o: /usr/include/sys/select.h /usr/include/bits/select.h
./src/nrg_refined.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
./src/nrg_refined.o: /usr/include/sys/sysmacros.h
./src/nrg_refined.o: /usr/include/bits/pthreadtypes.h /usr/include/unistd.h
./src/nrg_refined.o: /usr/include/bits/posix_opt.h
./src/nrg_refined.o: /usr/include/bits/confname.h /usr/include/getopt.h
./src/nrg_refined.o: /usr/include/stdlib.h /usr/include/alloca.h
./src/nrg_refined.o: /usr/include/string.h /usr/include/xlocale.h
./src/nrg_refined.o: ./src/printing_functions.h ./src/matrix_functions.h
./src/nrg_refined.o: ./src/ed.h ./src/hamiltonian.h
./src/nrg_refined.o: ./src/bethe_lapack_interface.h
./src/nrg_refined.o: ./src/hamiltonian_spin_functions.h
./src/nrg_refined.o: ./src/density_matrix.h
./src/nrg_spin.o: ./src/nrg_spin.h ./src/global.h /usr/include/stdio.h
./src/nrg_spin.o: /usr/include/features.h /usr/include/sys/cdefs.h
./src/nrg_spin.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
./src/nrg_spin.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
./src/nrg_spin.o: /usr/include/bits/typesizes.h /usr/include/libio.h
./src/nrg_spin.o: /usr/include/_G_config.h /usr/include/wchar.h
./src/nrg_spin.o: /usr/include/bits/stdio_lim.h
./src/nrg_spin.o: /usr/include/bits/sys_errlist.h ./src/matrix.h
./src/nrg_spin.o: ./src/cluster.h ./src/number_functions.h
./src/nrg_spin.o: ./src/math_utilities.h ./src/vector_utilities.h
./src/nrg_spin.o: ./src/bethe_isomorphism.h ./nauty/nausparse.h
./src/nrg_spin.o: ./nauty/nauty.h /usr/include/sys/types.h
./src/nrg_spin.o: /usr/include/time.h /usr/include/endian.h
./src/nrg_spin.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
./src/nrg_spin.o: /usr/include/sys/select.h /usr/include/bits/select.h
./src/nrg_spin.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
./src/nrg_spin.o: /usr/include/sys/sysmacros.h
./src/nrg_spin.o: /usr/include/bits/pthreadtypes.h /usr/include/unistd.h
./src/nrg_spin.o: /usr/include/bits/posix_opt.h /usr/include/bits/confname.h
./src/nrg_spin.o: /usr/include/getopt.h /usr/include/stdlib.h
./src/nrg_spin.o: /usr/include/alloca.h /usr/include/string.h
./src/nrg_spin.o: /usr/include/xlocale.h ./src/printing_functions.h
./src/nrg_spin.o: ./src/matrix_functions.h ./src/ed.h ./src/hamiltonian.h
./src/nrg_spin.o: ./src/bethe_lapack_interface.h
./src/nrg_spin.o: ./src/hamiltonian_spin_functions.h ./src/density_matrix.h
./src/number_functions.o: ./src/number_functions.h ./src/global.h
./src/number_functions.o: /usr/include/stdio.h /usr/include/features.h
./src/number_functions.o: /usr/include/sys/cdefs.h
./src/number_functions.o: /usr/include/bits/wordsize.h
./src/number_functions.o: /usr/include/gnu/stubs.h
./src/number_functions.o: /usr/include/gnu/stubs-64.h
./src/number_functions.o: /usr/include/bits/types.h
./src/number_functions.o: /usr/include/bits/typesizes.h /usr/include/libio.h
./src/number_functions.o: /usr/include/_G_config.h /usr/include/wchar.h
./src/number_functions.o: /usr/include/bits/stdio_lim.h
./src/number_functions.o: /usr/include/bits/sys_errlist.h ./src/matrix.h
./src/optimizer.o: ./src/hamiltonian_spin_functions.h ./src/global.h
./src/optimizer.o: /usr/include/stdio.h /usr/include/features.h
./src/optimizer.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
./src/optimizer.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
./src/optimizer.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
./src/optimizer.o: /usr/include/libio.h /usr/include/_G_config.h
./src/optimizer.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
./src/optimizer.o: /usr/include/bits/sys_errlist.h ./src/matrix.h
./src/printing_functions.o: ./src/printing_functions.h ./src/global.h
./src/printing_functions.o: /usr/include/stdio.h /usr/include/features.h
./src/printing_functions.o: /usr/include/sys/cdefs.h
./src/printing_functions.o: /usr/include/bits/wordsize.h
./src/printing_functions.o: /usr/include/gnu/stubs.h
./src/printing_functions.o: /usr/include/gnu/stubs-64.h
./src/printing_functions.o: /usr/include/bits/types.h
./src/printing_functions.o: /usr/include/bits/typesizes.h
./src/printing_functions.o: /usr/include/libio.h /usr/include/_G_config.h
./src/printing_functions.o: /usr/include/wchar.h
./src/printing_functions.o: /usr/include/bits/stdio_lim.h
./src/printing_functions.o: /usr/include/bits/sys_errlist.h ./src/matrix.h
./src/qddm.o: ./src/qddm.h ./src/lanczos_more.h ./src/hamiltonian.h
./src/qddm.o: ./src/global.h /usr/include/stdio.h /usr/include/features.h
./src/qddm.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
./src/qddm.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
./src/qddm.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
./src/qddm.o: /usr/include/libio.h /usr/include/_G_config.h
./src/qddm.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
./src/qddm.o: /usr/include/bits/sys_errlist.h ./src/matrix.h
./src/qddm.o: ./src/matrix_functions.h ./src/number_functions.h
./src/qddm.o: ./src/printing_functions.h ./src/bethe_lapack_interface.h
./src/qddm.o: ./src/math_utilities.h ./src/hamiltonian_spin_functions.h
./src/qddm.o: ./src/ed.h ./src/density_matrix.h
./src/search_for.o: ./src/global.h /usr/include/stdio.h
./src/search_for.o: /usr/include/features.h /usr/include/sys/cdefs.h
./src/search_for.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
./src/search_for.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
./src/search_for.o: /usr/include/bits/typesizes.h /usr/include/libio.h
./src/search_for.o: /usr/include/_G_config.h /usr/include/wchar.h
./src/search_for.o: /usr/include/bits/stdio_lim.h
./src/search_for.o: /usr/include/bits/sys_errlist.h ./src/matrix.h
./src/search_for.o: ./src/search_for.h
./src/singlet_triplet.o: ./src/singlet_triplet.h ./src/lanczos_more.h
./src/singlet_triplet.o: ./src/hamiltonian.h ./src/global.h
./src/singlet_triplet.o: /usr/include/stdio.h /usr/include/features.h
./src/singlet_triplet.o: /usr/include/sys/cdefs.h
./src/singlet_triplet.o: /usr/include/bits/wordsize.h
./src/singlet_triplet.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
./src/singlet_triplet.o: /usr/include/bits/types.h
./src/singlet_triplet.o: /usr/include/bits/typesizes.h /usr/include/libio.h
./src/singlet_triplet.o: /usr/include/_G_config.h /usr/include/wchar.h
./src/singlet_triplet.o: /usr/include/bits/stdio_lim.h
./src/singlet_triplet.o: /usr/include/bits/sys_errlist.h ./src/matrix.h
./src/singlet_triplet.o: ./src/matrix_functions.h ./src/number_functions.h
./src/singlet_triplet.o: ./src/printing_functions.h
./src/singlet_triplet.o: ./src/bethe_lapack_interface.h
./src/singlet_triplet.o: ./src/math_utilities.h
./src/singlet_triplet.o: ./src/hamiltonian_spin_functions.h ./src/ed.h
./src/singlet_triplet.o: ./src/density_matrix.h
./src/spin_functions.o: ./src/spin_functions.h ./src/global.h
./src/spin_functions.o: /usr/include/stdio.h /usr/include/features.h
./src/spin_functions.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
./src/spin_functions.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
./src/spin_functions.o: /usr/include/bits/types.h
./src/spin_functions.o: /usr/include/bits/typesizes.h /usr/include/libio.h
./src/spin_functions.o: /usr/include/_G_config.h /usr/include/wchar.h
./src/spin_functions.o: /usr/include/bits/stdio_lim.h
./src/spin_functions.o: /usr/include/bits/sys_errlist.h ./src/matrix.h
./src/spin_hole_model.o: ./src/spinless_model.h ./src/hamiltonian.h
./src/spin_hole_model.o: ./src/global.h /usr/include/stdio.h
./src/spin_hole_model.o: /usr/include/features.h /usr/include/sys/cdefs.h
./src/spin_hole_model.o: /usr/include/bits/wordsize.h
./src/spin_hole_model.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
./src/spin_hole_model.o: /usr/include/bits/types.h
./src/spin_hole_model.o: /usr/include/bits/typesizes.h /usr/include/libio.h
./src/spin_hole_model.o: /usr/include/_G_config.h /usr/include/wchar.h
./src/spin_hole_model.o: /usr/include/bits/stdio_lim.h
./src/spin_hole_model.o: /usr/include/bits/sys_errlist.h ./src/matrix.h
./src/spinless_boson_model.o: ./src/global.h /usr/include/stdio.h
./src/spinless_boson_model.o: /usr/include/features.h
./src/spinless_boson_model.o: /usr/include/sys/cdefs.h
./src/spinless_boson_model.o: /usr/include/bits/wordsize.h
./src/spinless_boson_model.o: /usr/include/gnu/stubs.h
./src/spinless_boson_model.o: /usr/include/gnu/stubs-64.h
./src/spinless_boson_model.o: /usr/include/bits/types.h
./src/spinless_boson_model.o: /usr/include/bits/typesizes.h
./src/spinless_boson_model.o: /usr/include/libio.h /usr/include/_G_config.h
./src/spinless_boson_model.o: /usr/include/wchar.h
./src/spinless_boson_model.o: /usr/include/bits/stdio_lim.h
./src/spinless_boson_model.o: /usr/include/bits/sys_errlist.h ./src/matrix.h
./src/spinless_boson_model.o: ./src/spinless_boson_model.h
./src/spinless_boson_model.o: ./src/hamiltonian.h
./src/spinless_boson_model.o: ./src/hamiltonian_spin_functions.h
./src/spinless_boson_model.o: ./src/number_functions.h ./src/search_for.h
./src/spinless_boson_model.o: ./src/printing_functions.h
./src/spinless_model.o: ./src/global.h /usr/include/stdio.h
./src/spinless_model.o: /usr/include/features.h /usr/include/sys/cdefs.h
./src/spinless_model.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
./src/spinless_model.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
./src/spinless_model.o: /usr/include/bits/typesizes.h /usr/include/libio.h
./src/spinless_model.o: /usr/include/_G_config.h /usr/include/wchar.h
./src/spinless_model.o: /usr/include/bits/stdio_lim.h
./src/spinless_model.o: /usr/include/bits/sys_errlist.h ./src/matrix.h
./src/spinless_model.o: ./src/spinless_model.h ./src/hamiltonian.h
./src/spinless_model.o: ./src/hamiltonian_spin_functions.h
./src/spinless_model.o: ./src/number_functions.h ./src/search_for.h
./src/spinless_model.o: ./src/printing_functions.h
./src/spin_model.o: ./src/global.h /usr/include/stdio.h
./src/spin_model.o: /usr/include/features.h /usr/include/sys/cdefs.h
./src/spin_model.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
./src/spin_model.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
./src/spin_model.o: /usr/include/bits/typesizes.h /usr/include/libio.h
./src/spin_model.o: /usr/include/_G_config.h /usr/include/wchar.h
./src/spin_model.o: /usr/include/bits/stdio_lim.h
./src/spin_model.o: /usr/include/bits/sys_errlist.h ./src/matrix.h
./src/spin_model.o: ./src/spin_model.h ./src/hamiltonian.h
./src/spin_model.o: ./src/hamiltonian_spin_functions.h
./src/spin_model.o: ./src/number_functions.h ./src/search_for.h
./src/spin_model.o: ./src/printing_functions.h ./src/spin_functions.h
./src/spin_waves.o: ./src/spin_waves.h ./src/global.h /usr/include/stdio.h
./src/spin_waves.o: /usr/include/features.h /usr/include/sys/cdefs.h
./src/spin_waves.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
./src/spin_waves.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
./src/spin_waves.o: /usr/include/bits/typesizes.h /usr/include/libio.h
./src/spin_waves.o: /usr/include/_G_config.h /usr/include/wchar.h
./src/spin_waves.o: /usr/include/bits/stdio_lim.h
./src/spin_waves.o: /usr/include/bits/sys_errlist.h ./src/matrix.h
./src/spin_waves.o: ./src/matrix_functions.h ./src/number_functions.h
./src/spin_waves.o: ./src/bethe_lapack_interface.h ./src/printing_functions.h
./src/spin_waves.o: ./src/cluster.h ./src/math_utilities.h
./src/spin_waves.o: ./src/vector_utilities.h ./src/bethe_isomorphism.h
./src/spin_waves.o: ./nauty/nausparse.h ./nauty/nauty.h
./src/spin_waves.o: /usr/include/sys/types.h /usr/include/time.h
./src/spin_waves.o: /usr/include/endian.h /usr/include/bits/endian.h
./src/spin_waves.o: /usr/include/bits/byteswap.h /usr/include/sys/select.h
./src/spin_waves.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
./src/spin_waves.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
./src/spin_waves.o: /usr/include/bits/pthreadtypes.h /usr/include/unistd.h
./src/spin_waves.o: /usr/include/bits/posix_opt.h
./src/spin_waves.o: /usr/include/bits/confname.h /usr/include/getopt.h
./src/spin_waves.o: /usr/include/stdlib.h /usr/include/alloca.h
./src/spin_waves.o: /usr/include/string.h /usr/include/xlocale.h
./src/sumi_ct.o: ./src/sumi_ct.h ./src/matrix.h ./src/number_functions.h
./src/sumi_ct.o: ./src/global.h /usr/include/stdio.h /usr/include/features.h
./src/sumi_ct.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
./src/sumi_ct.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
./src/sumi_ct.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
./src/sumi_ct.o: /usr/include/libio.h /usr/include/_G_config.h
./src/sumi_ct.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
./src/sumi_ct.o: /usr/include/bits/sys_errlist.h ./src/lanczos_more.h
./src/sumi_ct.o: ./src/hamiltonian.h ./src/matrix_functions.h
./src/sumi_ct.o: ./src/printing_functions.h ./src/bethe_lapack_interface.h
./src/sumi_ct.o: ./src/math_utilities.h ./src/hamiltonian_spin_functions.h
./src/sumi_ct.o: ./src/ed.h ./src/density_matrix.h
./src/tensor.o: ./src/tensor.h ./src/global.h /usr/include/stdio.h
./src/tensor.o: /usr/include/features.h /usr/include/sys/cdefs.h
./src/tensor.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
./src/tensor.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
./src/tensor.o: /usr/include/bits/typesizes.h /usr/include/libio.h
./src/tensor.o: /usr/include/_G_config.h /usr/include/wchar.h
./src/tensor.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
./src/tensor.o: ./src/matrix.h ./src/bethe_lapack_interface.h
./src/tests.o: ./src/tests.h ./src/global.h /usr/include/stdio.h
./src/tests.o: /usr/include/features.h /usr/include/sys/cdefs.h
./src/tests.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
./src/tests.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
./src/tests.o: /usr/include/bits/typesizes.h /usr/include/libio.h
./src/tests.o: /usr/include/_G_config.h /usr/include/wchar.h
./src/tests.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
./src/tests.o: ./src/matrix.h ./src/number_functions.h ./src/math_utilities.h
./src/tests.o: ./src/vector_utilities.h ./src/bethe_isomorphism.h
./src/tests.o: ./nauty/nausparse.h ./nauty/nauty.h /usr/include/sys/types.h
./src/tests.o: /usr/include/time.h /usr/include/endian.h
./src/tests.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
./src/tests.o: /usr/include/sys/select.h /usr/include/bits/select.h
./src/tests.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
./src/tests.o: /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h
./src/tests.o: /usr/include/unistd.h /usr/include/bits/posix_opt.h
./src/tests.o: /usr/include/bits/confname.h /usr/include/getopt.h
./src/tests.o: /usr/include/stdlib.h /usr/include/alloca.h
./src/tests.o: /usr/include/string.h /usr/include/xlocale.h
./src/tests.o: ./src/printing_functions.h ./src/cluster.h ./src/search_for.h
./src/tests.o: ./src/data_io.h ./src/spin_waves.h ./src/matrix_functions.h
./src/tests.o: ./src/laplacian.h ./src/chris_fm_phi.h
./src/tfi_model.o: ./src/global.h /usr/include/stdio.h
./src/tfi_model.o: /usr/include/features.h /usr/include/sys/cdefs.h
./src/tfi_model.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
./src/tfi_model.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
./src/tfi_model.o: /usr/include/bits/typesizes.h /usr/include/libio.h
./src/tfi_model.o: /usr/include/_G_config.h /usr/include/wchar.h
./src/tfi_model.o: /usr/include/bits/stdio_lim.h
./src/tfi_model.o: /usr/include/bits/sys_errlist.h ./src/matrix.h
./src/tfi_model.o: ./src/tfi_model.h ./src/hamiltonian.h
./src/tfi_model.o: ./src/hamiltonian_spin_functions.h
./src/tfi_model.o: ./src/number_functions.h ./src/search_for.h
./src/tree_dmrg_spin.o: ./src/tree_dmrg_spin.h ./src/global.h
./src/tree_dmrg_spin.o: /usr/include/stdio.h /usr/include/features.h
./src/tree_dmrg_spin.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
./src/tree_dmrg_spin.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
./src/tree_dmrg_spin.o: /usr/include/bits/types.h
./src/tree_dmrg_spin.o: /usr/include/bits/typesizes.h /usr/include/libio.h
./src/tree_dmrg_spin.o: /usr/include/_G_config.h /usr/include/wchar.h
./src/tree_dmrg_spin.o: /usr/include/bits/stdio_lim.h
./src/tree_dmrg_spin.o: /usr/include/bits/sys_errlist.h ./src/matrix.h
./src/tree_dmrg_spin.o: ./src/cluster.h ./src/number_functions.h
./src/tree_dmrg_spin.o: ./src/math_utilities.h ./src/vector_utilities.h
./src/tree_dmrg_spin.o: ./src/bethe_isomorphism.h ./nauty/nausparse.h
./src/tree_dmrg_spin.o: ./nauty/nauty.h /usr/include/sys/types.h
./src/tree_dmrg_spin.o: /usr/include/time.h /usr/include/endian.h
./src/tree_dmrg_spin.o: /usr/include/bits/endian.h
./src/tree_dmrg_spin.o: /usr/include/bits/byteswap.h
./src/tree_dmrg_spin.o: /usr/include/sys/select.h /usr/include/bits/select.h
./src/tree_dmrg_spin.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
./src/tree_dmrg_spin.o: /usr/include/sys/sysmacros.h
./src/tree_dmrg_spin.o: /usr/include/bits/pthreadtypes.h
./src/tree_dmrg_spin.o: /usr/include/unistd.h /usr/include/bits/posix_opt.h
./src/tree_dmrg_spin.o: /usr/include/bits/confname.h /usr/include/getopt.h
./src/tree_dmrg_spin.o: /usr/include/stdlib.h /usr/include/alloca.h
./src/tree_dmrg_spin.o: /usr/include/string.h /usr/include/xlocale.h
./src/tree_dmrg_spin.o: ./src/printing_functions.h ./src/matrix_functions.h
./src/tree_dmrg_spin.o: ./src/spin_model.h ./src/hamiltonian.h ./src/block.h
./src/tree_dmrg_spin.o: ./src/spin_functions.h ./src/tmp_info.h
./src/tree_dmrg_spin.o: ./src/universe.h ./src/tree_nrg_spin.h ./src/tensor.h
./src/tree_dmrg_spin.o: ./src/density_matrix.h ./src/bethe_lapack_interface.h
./src/tree_nrg_spin.o: ./src/tree_nrg_spin.h ./src/global.h
./src/tree_nrg_spin.o: /usr/include/stdio.h /usr/include/features.h
./src/tree_nrg_spin.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
./src/tree_nrg_spin.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
./src/tree_nrg_spin.o: /usr/include/bits/types.h
./src/tree_nrg_spin.o: /usr/include/bits/typesizes.h /usr/include/libio.h
./src/tree_nrg_spin.o: /usr/include/_G_config.h /usr/include/wchar.h
./src/tree_nrg_spin.o: /usr/include/bits/stdio_lim.h
./src/tree_nrg_spin.o: /usr/include/bits/sys_errlist.h ./src/matrix.h
./src/tree_nrg_spin.o: ./src/cluster.h ./src/number_functions.h
./src/tree_nrg_spin.o: ./src/math_utilities.h ./src/vector_utilities.h
./src/tree_nrg_spin.o: ./src/bethe_isomorphism.h ./nauty/nausparse.h
./src/tree_nrg_spin.o: ./nauty/nauty.h /usr/include/sys/types.h
./src/tree_nrg_spin.o: /usr/include/time.h /usr/include/endian.h
./src/tree_nrg_spin.o: /usr/include/bits/endian.h
./src/tree_nrg_spin.o: /usr/include/bits/byteswap.h /usr/include/sys/select.h
./src/tree_nrg_spin.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
./src/tree_nrg_spin.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
./src/tree_nrg_spin.o: /usr/include/bits/pthreadtypes.h /usr/include/unistd.h
./src/tree_nrg_spin.o: /usr/include/bits/posix_opt.h
./src/tree_nrg_spin.o: /usr/include/bits/confname.h /usr/include/getopt.h
./src/tree_nrg_spin.o: /usr/include/stdlib.h /usr/include/alloca.h
./src/tree_nrg_spin.o: /usr/include/string.h /usr/include/xlocale.h
./src/tree_nrg_spin.o: ./src/printing_functions.h ./src/matrix_functions.h
./src/tree_nrg_spin.o: ./src/spin_model.h ./src/hamiltonian.h ./src/block.h
./src/tree_nrg_spin.o: ./src/spin_functions.h ./src/tmp_info.h
./src/tree_nrg_spin.o: ./src/universe.h ./src/tensor.h ./src/ed.h
./src/tree_nrg_spin.o: ./src/bethe_lapack_interface.h
./src/tree_nrg_spin.o: ./src/hamiltonian_spin_functions.h
./src/tree_nrg_spin.o: ./src/density_matrix.h
./src/two_dangling_analysis.o: ./src/two_dangling_analysis.h
./src/two_dangling_analysis.o: ./src/lanczos_more.h ./src/hamiltonian.h
./src/two_dangling_analysis.o: ./src/global.h /usr/include/stdio.h
./src/two_dangling_analysis.o: /usr/include/features.h
./src/two_dangling_analysis.o: /usr/include/sys/cdefs.h
./src/two_dangling_analysis.o: /usr/include/bits/wordsize.h
./src/two_dangling_analysis.o: /usr/include/gnu/stubs.h
./src/two_dangling_analysis.o: /usr/include/gnu/stubs-64.h
./src/two_dangling_analysis.o: /usr/include/bits/types.h
./src/two_dangling_analysis.o: /usr/include/bits/typesizes.h
./src/two_dangling_analysis.o: /usr/include/libio.h /usr/include/_G_config.h
./src/two_dangling_analysis.o: /usr/include/wchar.h
./src/two_dangling_analysis.o: /usr/include/bits/stdio_lim.h
./src/two_dangling_analysis.o: /usr/include/bits/sys_errlist.h ./src/matrix.h
./src/two_dangling_analysis.o: ./src/matrix_functions.h
./src/two_dangling_analysis.o: ./src/number_functions.h
./src/two_dangling_analysis.o: ./src/printing_functions.h
./src/two_dangling_analysis.o: ./src/bethe_lapack_interface.h
./src/two_dangling_analysis.o: ./src/math_utilities.h
./src/two_dangling_analysis.o: ./src/hamiltonian_spin_functions.h ./src/ed.h
./src/two_dangling_analysis.o: ./src/density_matrix.h
./src/vector_utilities.o: ./src/vector_utilities.h ./src/global.h
./src/vector_utilities.o: /usr/include/stdio.h /usr/include/features.h
./src/vector_utilities.o: /usr/include/sys/cdefs.h
./src/vector_utilities.o: /usr/include/bits/wordsize.h
./src/vector_utilities.o: /usr/include/gnu/stubs.h
./src/vector_utilities.o: /usr/include/gnu/stubs-64.h
./src/vector_utilities.o: /usr/include/bits/types.h
./src/vector_utilities.o: /usr/include/bits/typesizes.h /usr/include/libio.h
./src/vector_utilities.o: /usr/include/_G_config.h /usr/include/wchar.h
./src/vector_utilities.o: /usr/include/bits/stdio_lim.h
./src/vector_utilities.o: /usr/include/bits/sys_errlist.h ./src/matrix.h
./src/xxz_model.o: ./src/global.h /usr/include/stdio.h
./src/xxz_model.o: /usr/include/features.h /usr/include/sys/cdefs.h
./src/xxz_model.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
./src/xxz_model.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
./src/xxz_model.o: /usr/include/bits/typesizes.h /usr/include/libio.h
./src/xxz_model.o: /usr/include/_G_config.h /usr/include/wchar.h
./src/xxz_model.o: /usr/include/bits/stdio_lim.h
./src/xxz_model.o: /usr/include/bits/sys_errlist.h ./src/matrix.h
./src/xxz_model.o: ./src/xxz_model.h ./src/hamiltonian.h
./src/xxz_model.o: ./src/hamiltonian_spin_functions.h
./src/xxz_model.o: ./src/number_functions.h ./src/search_for.h
./src/xxz_model.o: ./src/printing_functions.h
./nauty/nauty.o: ./nauty/nauty.h /usr/include/stdio.h /usr/include/features.h
./nauty/nauty.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
./nauty/nauty.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
./nauty/nauty.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
./nauty/nauty.o: /usr/include/libio.h /usr/include/_G_config.h
./nauty/nauty.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
./nauty/nauty.o: /usr/include/bits/sys_errlist.h /usr/include/sys/types.h
./nauty/nauty.o: /usr/include/time.h /usr/include/endian.h
./nauty/nauty.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
./nauty/nauty.o: /usr/include/sys/select.h /usr/include/bits/select.h
./nauty/nauty.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
./nauty/nauty.o: /usr/include/sys/sysmacros.h
./nauty/nauty.o: /usr/include/bits/pthreadtypes.h /usr/include/unistd.h
./nauty/nauty.o: /usr/include/bits/posix_opt.h /usr/include/bits/confname.h
./nauty/nauty.o: /usr/include/getopt.h /usr/include/stdlib.h
./nauty/nauty.o: /usr/include/alloca.h /usr/include/string.h
./nauty/nauty.o: /usr/include/xlocale.h
./nauty/nautil.o: ./nauty/nauty.h /usr/include/stdio.h
./nauty/nautil.o: /usr/include/features.h /usr/include/sys/cdefs.h
./nauty/nautil.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
./nauty/nautil.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
./nauty/nautil.o: /usr/include/bits/typesizes.h /usr/include/libio.h
./nauty/nautil.o: /usr/include/_G_config.h /usr/include/wchar.h
./nauty/nautil.o: /usr/include/bits/stdio_lim.h
./nauty/nautil.o: /usr/include/bits/sys_errlist.h /usr/include/sys/types.h
./nauty/nautil.o: /usr/include/time.h /usr/include/endian.h
./nauty/nautil.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
./nauty/nautil.o: /usr/include/sys/select.h /usr/include/bits/select.h
./nauty/nautil.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
./nauty/nautil.o: /usr/include/sys/sysmacros.h
./nauty/nautil.o: /usr/include/bits/pthreadtypes.h /usr/include/unistd.h
./nauty/nautil.o: /usr/include/bits/posix_opt.h /usr/include/bits/confname.h
./nauty/nautil.o: /usr/include/getopt.h /usr/include/stdlib.h
./nauty/nautil.o: /usr/include/alloca.h /usr/include/string.h
./nauty/nautil.o: /usr/include/xlocale.h
./nauty/nausparse.o: ./nauty/nausparse.h ./nauty/nauty.h /usr/include/stdio.h
./nauty/nausparse.o: /usr/include/features.h /usr/include/sys/cdefs.h
./nauty/nausparse.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
./nauty/nausparse.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
./nauty/nausparse.o: /usr/include/bits/typesizes.h /usr/include/libio.h
./nauty/nausparse.o: /usr/include/_G_config.h /usr/include/wchar.h
./nauty/nausparse.o: /usr/include/bits/stdio_lim.h
./nauty/nausparse.o: /usr/include/bits/sys_errlist.h /usr/include/sys/types.h
./nauty/nausparse.o: /usr/include/time.h /usr/include/endian.h
./nauty/nausparse.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
./nauty/nausparse.o: /usr/include/sys/select.h /usr/include/bits/select.h
./nauty/nausparse.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
./nauty/nausparse.o: /usr/include/sys/sysmacros.h
./nauty/nausparse.o: /usr/include/bits/pthreadtypes.h /usr/include/unistd.h
./nauty/nausparse.o: /usr/include/bits/posix_opt.h
./nauty/nausparse.o: /usr/include/bits/confname.h /usr/include/getopt.h
./nauty/nausparse.o: /usr/include/stdlib.h /usr/include/alloca.h
./nauty/nausparse.o: /usr/include/string.h /usr/include/xlocale.h
