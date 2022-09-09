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
ifeq (${sanny},1)
flags += -fsanitize=undefined,address -fstack-protector-all
endif

main:
	${cc} -std=c++20 -g3 -O3 ${flags} ${compflags} -I${cvdf} main.cc -o ${target}

run: main
	./${target}

clean:
	rm -f ${target}
