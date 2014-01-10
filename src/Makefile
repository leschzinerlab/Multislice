CC=gcc
CFLAGS=-I. -m32 
LIBS=-lm
AUTOSLICE= autoslic_Q.o tiffsubs.o slicelib.o fft2dc.o memory.o 
POISSON= add_poisson.o tiffsubs.o slicelib.o fft2dc.o memory.o
IMAGE= image_phase.o tiffsubs.o slicelib.o fft2dc.o memory.o

%.o: %.c  
	$(CC) -c -o $@ $^ $(CFLAGS) $(LIBS)

all: autoslice_Q add_poisson image

autoslice_Q: $(AUTOSLICE)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

add_poisson: $(POISSON)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

image: $(IMAGE)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -f *.o core

