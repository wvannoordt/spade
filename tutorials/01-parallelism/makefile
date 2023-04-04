target := vdf.x
cc := $(shell which mpicxx)


sanny := 0

flags :=
flags += -std=c++20 -O3 -fconcepts-diagnostics-depth=3 -fmax-errors=4

ifeq (${symbols},1)
flags += -g
endif

ifeq (${sanny},1)
flags += -fsanitize=undefined,address -fstack-protector-all
endif

main:
	${cc} ${flags} -I${SPADE}/src main.cc -o ${target}

run: main
	./${target}

clean:
	rm -f ${target}