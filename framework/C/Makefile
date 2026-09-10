# Vivian framework (C23) -- clang-only build
CC      := clang
AR      := llvm-ar
CFLAGS  := -std=c23 -O2 -Wall -Wextra -Iinclude
LIBSRC  := src/vivi.c src/dna.c src/genome.c src/chromosome.c src/cell.c src/organism.c
HDRS    := $(wildcard include/vivi/*.h)
LIBOBJ  := $(LIBSRC:.c=.o)

.PHONY: all lib test demo density asan clean

all: vivi.a vivi_test.exe vivi_demo.exe

# framework only -- link vivi.a into your own entrypoint
lib: vivi.a

vivi.a: $(LIBOBJ) $(HDRS)
	$(AR) rcs $@ $(LIBOBJ)

src/%.o: src/%.c $(HDRS)
	$(CC) $(CFLAGS) -c $< -o $@

vivi_test.exe: test/test.c $(LIBSRC) $(HDRS)
	$(CC) $(CFLAGS) -o $@ test/test.c $(LIBSRC)

vivi_demo.exe: demo/demo.c $(LIBSRC) $(HDRS)
	$(CC) $(CFLAGS) -o $@ demo/demo.c $(LIBSRC)

density.exe: test/density.c $(LIBSRC) $(HDRS)
	$(CC) $(CFLAGS) -o $@ test/density.c $(LIBSRC)

test: vivi_test.exe
	./vivi_test.exe

demo: vivi_demo.exe
	./vivi_demo.exe

density: density.exe
	./density.exe

asan: 
	$(CC) -std=c23 -O1 -g -fsanitize=address,undefined -Iinclude -o vivi_test_asan test/test.c $(LIBSRC)
	./vivi_test_asan

clean:
	rm -f vivi.a vivi_test.exe vivi_demo.exe density.exe vivi_test_asan src/*.o
