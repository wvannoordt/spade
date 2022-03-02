cvdf   := $(realpath ../../src)
target := vdf.x
main:
	g++ -std=c++20 -O3 -I${cvdf} main.cc -o ${target}

run: main
	./${target}

clean:
	rm -f ${target}