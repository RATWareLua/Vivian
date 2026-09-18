CC      := clang
CFLAGS  := -std=c23 -O2 -Wall -Wextra -I../../include
LIBSRC  := ../../src/vivi.c ../../src/dna.c ../../src/genome.c \
           ../../src/chromosome.c ../../src/channel.c ../../src/pool.c \
           ../../src/parity.c ../../src/rs.c ../../src/cell.c ../../src/organism.c \
           ../../src/model.c

.PHONY: all clean

all: vivi_calib

vivi_calib: vivi_calib.c $(LIBSRC) $(wildcard ../../include/vivi/*.h)
	$(CC) $(CFLAGS) -o $@ vivi_calib.c $(LIBSRC) -lm

clean:
	rm -f vivi_calib vivi_calib.exe
