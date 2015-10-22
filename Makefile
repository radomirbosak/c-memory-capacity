#WALL=-Wall
#STD=-std=c99

#LIBS  = -lkernel32 -luser32 -lgdi32 -lopengl32
#MATHLIBS = -lgsl -lgslcblas -lm
MATHLIBS = -L/usr/lib/atlas-base -lgsl -lcblas -latlas -lm
#PIC = -fPIC
#WALL = -Wall
CFLAGS = -c $(STD) $(WALL) $(PIC)

all: Foo

%.o: %.c
	gcc -O3 -march=native $(CFLAGS) $< -o $@

Foo: foo.o matrix.o randomvar.o
	gcc -O3 -march=native -o $@ $^ $(MATHLIBS)