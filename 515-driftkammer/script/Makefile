default:
	g++ -std=c++21 -Wall -Wpedantic -O2 `root-config --libs --ldflags` `root-config --cflags` analysis.cpp -o a


run:
	make default
	./a ../data/langzeit.root
