cvdf   := $(realpath ../../../src)
target := vdf.x

ifndef MPI_ENABLE
MPI_ENABLE := 1
endif

ifeq (${MPI_ENABLE},1)
cc := $(shell which mpicxx)
else
cc := $(shell which g++)
endif

compflags :=
compflags += -DMPI_ENABLE=${MPI_ENABLE}

flags :=
flags += -fconcepts-diagnostics-depth=3

main:
	${cc} -std=c++20 -g -O3 ${flags} ${compflags} -I${cvdf} -c diffkernel.cc -o diffkernel.o
	${cc} -std=c++20 -g -O3 ${flags} ${compflags} -I${cvdf} -c main.cc -o main.o
	${cc} *.o -o ${target}

run: main
	./${target}

clean:
	rm -f ${target}
