all: nonograms nonograms-v1

# -std=c++11
%.o: %.cc
	g++ -Wall -g -c -O6 -o $@ $<

nonograms:	main.o
	g++ -Wall -g -o $@ $^
	
nonograms-v1:	nonograms-v1.o
	g++ -Wall -g -o $@ $^

clean:
	rm -f main.o nonograms-v1.o nonograms nonograms-v1

