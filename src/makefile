CC = g++
CXXFLAGS = -O3 -std=c++11

OBJS = visit_writer.o IOdict.o main.o flDriver.o lung.o duct.o lobule.o gas.o

flPROG: $(OBJS)
	$(CC) $(CXXFLAGS) -o flPROG $(OBJS)

moveup:
	cp flPROG ../flPROG

clean:
	rm flPROG
	rm -f $(OBJS)