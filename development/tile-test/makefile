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

ifeq (${use_cuda},1)
cc_host = ${cc}
cc := $(shell which nvcc) -ccbin=${cc_host}
endif

compflags :=
compflags += -DMPI_ENABLE=${MPI_ENABLE}
compflags += -DUSE_SOURCE_LOCATION=0

sanny := 1

flags :=
flags += -std=c++20 -g -O3 -fconcepts-diagnostics-depth=3 -fmax-errors=4
ifeq (${sanny},1)
flags += -fsanitize=undefined,address -fstack-protector-all
endif

ifeq (${use_cuda},1)
flags_old := ${flags}
flags = -x cu -std=c++20 -Xcompiler "${flags_old}"
endif

main:
	${cc} ${flags} ${compflags} -I${SPADE} -I${PTL}/include main.cc -o ${target} -L${PTL}/lib -lPTL

run: main
	./${target}

clean:
	rm -f ${target}
