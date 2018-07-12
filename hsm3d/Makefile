
OBJECTS = common.o coordinates.o hough3d.o io.o scan3d.o


# Uncomment one of these two depending on your platform (Linux/Os X)

# uncomment these lines in Linux
LIBRARY = libh3d.so
LIBFLAGS = -shared

# uncomment these lines in Os X
#LIBRARY = libh3d.dylib
#LIBFLAGS = -dynamiclib -ljson

EXECUTABLES = testcoordinates  testcomplete testvisibility testhough3d

C++ = g++ -std=c++17


LDFLAGS =  -L. -L/usr/local/lib  -L/home/aditya/ANN/ann_1.1.2/lib -lm -lh3d -ljson-c -lANN

all: $(LIBRARY) $(EXECUTABLES)

$(LIBRARY): $(OBJECTS)
	$(C++) $(LIBFLAGS) $(OBJECTS) -o $(LIBRARY)

testcoordinates: $(LIBRARY) testcoordinates.o
	$(C++)   testcoordinates.o -o testcoordinates $(LDFLAGS)

testvisibility: $(LIBRARY) testvisibility.o
	$(C++)   testvisibility.o -o testvisibility $(LDFLAGS)

testhough3d: $(LIBRARY) testhough3d.o
	$(C++)   $(LDFLAGS) testhough3d.o  -o testhough3d $(LDFLAGS)

testcomplete: $(LIBRARY) testcomplete.o
	$(C++)  testcomplete.o  -o testcomplete $(LDFLAGS)



FLAGS= -Wall   -I/usr/local/include -I/home/aditya/ANN/ann_1.1.2/include -fPIC  -g -O2 -DDOUBLE_PRECISION 


%.o : %.cpp
	$(C++) $(CFLAGS)  $(FLAGS) -c $<

clean:
	rm -rf *.o *~ $(EXECUTABLES) $(LIBRARY)

doc:
	doxygen Doxyfile

cleandoc:
	rm -rf html
