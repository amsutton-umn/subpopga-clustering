CC=gcc
CFLAGS=-ggdb -Wall
#CFLAGS=-O3 -DNDEBUG
# use -DNDEBUG to disable assertions
LDFLAGS=-llzma -lm -lglpk
BIN=subpopga
SRCS := $(wildcard src/*.c)
OBJS := $(patsubst %.c,%.o,$(SRCS))

$(BIN): $(OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

%.o : %.c
	$(CC) -c $(CFLAGS) $< -o $@

.PHONY: clean
clean:
	rm -f $(BIN) $(OBJS)
