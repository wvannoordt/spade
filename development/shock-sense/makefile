target := vdf.x

ifndef MPI_ENABLE
MPI_ENABLE := 1
endif

ifeq (${MPI_ENABLE},1)
cc := $(shell which mpicxx)
else
cc := $(shell which g++)
endif

ifndef use_cuda
use_cuda := 0
endif

compflags :=
compflags += -DMPI_ENABLE=${MPI_ENABLE}
compflags += -DUSE_SOURCE_LOCATION=0

sanny := 0

#Note that including -g dramatically increases compile time
flags :=
flags += -std=c++20 -O3 -fconcepts-diagnostics-depth=3 -fmax-errors=4
ifeq (${sanny},1)
flags += -fsanitize=undefined,address -fstack-protector-all
endif

incl :=
incl += -I${SPADE}/src
incl += -I${SYMD}/src

link :=
main:
	${cc} ${flags} ${compflags} ${incl} main.cc -o ${target} ${link}

run: main
	./${target}

clean:
	rm -f ${target}
