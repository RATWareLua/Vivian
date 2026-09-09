# Genetic storage framework (C23) -- clang-only build
CC      := clang
AR      := llvm-ar
CFLAGS  := -std=c23 -O2 -Wall -Wextra -Iinclude
LIBSRC  := src/vivi.c src/dna.c src/genome.c src/chromosome.c src/cell.c src/organism.c
HDRS    := $(wildcard include/gsy/*.h)
LIBOBJ  := $(LIBSRC:.c=.o)

.PHONY: all lib test demo clean

all: libgsy.a vivi_test.exe vivi_demo.exe

# framework only -- link vivi.lib into your own entrypoint
lib: libgsy.a

libgsy.a: $(LIBOBJ) $(HDRS)
	$(AR) rcs $@ $(LIBOBJ)

src/%.o: src/%.c $(HDRS)
	$(CC) $(CFLAGS) -c $< -o $@

vivi_test.exe: test/test.c $(LIBSRC) $(HDRS)
	$(CC) $(CFLAGS) -o $@ test/test.c $(LIBSRC)

vivi_demo.exe: demo/demo.c $(LIBSRC) $(HDRS)
	$(CC) $(CFLAGS) -o $@ demo/demo.c $(LIBSRC)

test: vivi_test.exe
	./vivi_test.exe

demo: vivi_demo.exe
	./vivi_demo.exe

clean:
	rm -f libgsy.a vivi_test.exe vivi_demo.exe src/*.o

