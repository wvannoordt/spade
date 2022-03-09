cvdf   := $(realpath ../../src)
target := vdf.x
main:
	g++ -std=c++20 -g -O0 -I${cvdf} main.cc -o ${target}

run: main
	./${target}

clean:
	rm -f ${target}
