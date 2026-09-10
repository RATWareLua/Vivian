# Vivian framework (C23) -- clang-only build
CC      := clang
AR      := llvm-ar
CFLAGS  := -std=c23 -O2 -Wall -Wextra -Iinclude
LIBSRC  := src/vivi.c src/dna.c src/genome.c src/chromosome.c src/channel.c src/pool.c src/parity.c src/cell.c src/organism.c
HDRS    := $(wildcard include/vivi/*.h)
LIBOBJ  := $(LIBSRC:.c=.o)

FUZZ_BIN := $(patsubst test/fuzz/%.c,vivi_fuzz_%,$(wildcard test/fuzz/*.c))

.PHONY: all lib test demo density sim research xcheck fuzz fuzz-smoke asan clean

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

# deterministic byte-parity scenarios shared with framework/luau/xcheck.lua;
# CI diffs the two outputs (portable -- runs anywhere clang runs)
xcheck.exe: test/xcheck.c $(LIBSRC) $(HDRS)
	$(CC) $(CFLAGS) -o $@ test/xcheck.c $(LIBSRC)

xcheck: xcheck.exe
	@./xcheck.exe

# full parameter sweeps -> research_whole.csv, research_access.csv, research_library.csv
research: sim.exe
	./sim.exe --size 1024 --gene-raw 64 --trials 1000 --p-sub 0 --header > research_whole.csv
	@for p in 0.0002 0.0005 0.001 0.002 0.004; do \
		./sim.exe --size 1024 --gene-raw 64 --trials 1000 --p-sub $$p >> research_whole.csv; done
	./sim.exe --size 1024 --gene-raw 64 --access --trials 1000 --p-sub 0.001 --p-access 0 --header > research_access.csv
	@for a in 0.01 0.05 0.1; do \
		./sim.exe --size 1024 --gene-raw 64 --access --trials 1000 --p-sub 0.001 --p-access $$a --p-cross 0.01 >> research_access.csv; done
	./sim.exe --size 1024 --gene-raw 64 --access --trials 1000 --p-sub 0.001 --p-access 0 --p-cross 0.01 --p-primer 0.05 >> research_access.csv
	./sim.exe --size 1024 --gene-raw 64 --access --trials 1000 --p-sub 0.001 --p-access 0.05 --p-cross 0.01 --p-primer 0.05 >> research_access.csv
	./sim.exe --size 1024 --gene-raw 64 --library --trials 1000 --p-sub 0.001 --header > research_library.csv
	@for k in 2 3 4; do \
		./sim.exe --size 1024 --gene-raw 64 --library --trials 1000 --p-sub 0.001 --replicas $$k >> research_library.csv; done
	@for m in 2 4 8; do \
		./sim.exe --size 1024 --gene-raw 64 --library --trials 1000 --p-sub 0.001 --parity $$m >> research_library.csv; done
	./sim.exe --size 1024 --gene-raw 64 --library --trials 1000 --p-sub 0.001 --replicas 2 --parity 4 >> research_library.csv
	@echo "wrote research_whole.csv research_access.csv research_library.csv"

# libFuzzer harnesses (POSIX only; needs clang -fsanitize=fuzzer)
fuzz: $(FUZZ_BIN)

vivi_fuzz_%: test/fuzz/%.c $(LIBSRC) $(HDRS)
	$(CC) -std=c23 -O1 -g -fsanitize=fuzzer,address,undefined -Iinclude -o $@ $< $(LIBSRC)

fuzz-smoke: $(FUZZ_BIN)
	@for f in $(FUZZ_BIN); do \
		echo "== $$f"; \
		./$$f -max_total_time=10 -max_len=2048 -print_final_stats=1 || exit 1; \
	done

clean:
	rm -f vivi.a vivi_test.exe vivi_demo.exe density.exe sim.exe xcheck.exe vivi_test_asan $(FUZZ_BIN) research_*.csv src/*.o
