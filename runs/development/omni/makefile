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
	${cc} -std=c++20 -g -O3 -fmax-errors=4 ${flags} ${compflags} -I${SPADE} -I${PTL}/include main.cc -o ${target} -L${PTL}/lib -lPTL

run: main
	./${target}

clean:
	rm -f ${target}
