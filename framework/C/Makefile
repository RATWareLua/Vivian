# Vivian framework (C23) -- clang-only build
CC      := clang
AR      := llvm-ar
CFLAGS  := -std=c23 -O2 -Wall -Wextra -Iinclude
LIBSRC  := src/vivi.c src/dna.c src/genome.c src/chromosome.c src/channel.c src/pool.c src/parity.c src/cell.c src/organism.c
HDRS    := $(wildcard include/vivi/*.h)
LIBOBJ  := $(LIBSRC:.c=.o)

.PHONY: all lib test demo density sim asan clean

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

sim.exe: test/sim.c $(LIBSRC) $(HDRS)
	$(CC) $(CFLAGS) -o $@ test/sim.c $(LIBSRC)

sim: sim.exe
	./sim.exe --size 256 --trials 1000 --p-sub 0.0005 --header
	./sim.exe --size 256 --trials 1000 --p-sub 0.001
	./sim.exe --size 256 --trials 1000 --p-sub 0.002
	./sim.exe --size 256 --trials 1000 --p-sub 0.004 --p-drop 0.01
	./sim.exe --size 1024 --gene-raw 64 --trials 1000 --p-sub 0.001 --header
	./sim.exe --size 1024 --gene-raw 64 --access --trials 1000 --p-sub 0.001 --p-access 0.05 --p-cross 0.01
	./sim.exe --size 1024 --gene-raw 64 --library --trials 1000 --p-sub 0.001 --header
	./sim.exe --size 1024 --gene-raw 64 --library --trials 1000 --p-sub 0.001 --replicas 3
	./sim.exe --size 1024 --gene-raw 64 --library --trials 1000 --p-sub 0.001 --parity 4
	./sim.exe --size 1024 --gene-raw 64 --library --trials 1000 --p-sub 0.001 --replicas 2 --parity 4

asan: 
	$(CC) -std=c23 -O1 -g -fsanitize=address,undefined -DVIVI_TEST_TRACK -Iinclude -o vivi_test_asan test/test.c $(LIBSRC)
	./vivi_test_asan

clean:
	rm -f vivi.a vivi_test.exe vivi_demo.exe density.exe sim.exe vivi_test_asan src/*.o
