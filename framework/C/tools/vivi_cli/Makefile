CC      := clang
CFLAGS  := -std=c23 -O2 -Wall -Wextra -I../../include
LIBSRC  := ../../src/vivi.c ../../src/dna.c ../../src/genome.c \
           ../../src/chromosome.c ../../src/channel.c ../../src/pool.c \
           ../../src/parity.c ../../src/cell.c ../../src/organism.c

.PHONY: all clean

all: vivi_cli

vivi_cli: vivi_cli.c $(LIBSRC) $(wildcard ../../include/vivi/*.h)
	$(CC) $(CFLAGS) -o $@ vivi_cli.c $(LIBSRC)

clean:
	rm -f vivi_cli vivi_cli.exe
