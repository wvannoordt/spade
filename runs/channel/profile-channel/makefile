spade  := $(realpath ../../../src)
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
flags += -fsanitize=undefined,address

iflags := 
iflags += -I${spade}
iflags += -I${PTL}/include

lflags := 
lflags += -L${PTL}/lib -lPTL

main:
	${cc} -std=c++20 -g -O0 ${flags} ${compflags} ${iflags} main.cc -o ${target} ${lflags}

run: main
	./${target}

clean:
	rm -f prf-ct.x prf-gc.x
