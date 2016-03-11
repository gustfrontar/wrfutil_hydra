MPICC=mpifccpx

all:
	${MPICC} -o spawn spawn.c
	${MPICC} -o dummy-mpi dummy-mpi.c

clean:
	@if test -e spawn; then\
	  rm spawn *~;\
	fi
