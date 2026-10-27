.PHONY: all clean

NAME = tp

CC = gcc
CFLAGS += -std=c89 -O2 -s
CFLAGS += -Wall -Wextra -pedantic

ifeq ($(OS),Windows_NT)
EXE = $(NAME).exe
else
EXE = $(NAME)
endif

all: $(EXE)

clean:
	rm -f $(EXE)

$(EXE): *.c
	$(CC) -o $@ $^ $(CFLAGS)
