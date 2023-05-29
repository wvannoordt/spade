spade   := ${SPADE}
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


iflags :=
iflags += -I${spade}/src
iflags += -I${SCIDF}/src
iflags += -I${HYWALL_PATH}/include

main:
	${cc} --version
	${cc} -std=c++20 -g -O3 ${flags} ${compflags} ${iflags} main.cc -o ${target} ${lflags}

run: main
	./${target}

raw:
	${cc} -E -std=c++20 -g -O3 ${flags} ${compflags} ${iflags} main.cc > raw.cc

clean:
	${MAKE} -C ${HYWALL_PATH} -f makefile clean
	rm -f ${target}
