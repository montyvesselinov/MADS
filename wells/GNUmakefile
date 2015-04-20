PROG = wells
CC = gcc 
CFLAGS = -Wall -lm
LDFLAGS = -lm
OBJ = wells.o
$(PROG): $(OBJ)
	$(CC) $(OBJ) $(LDFLAGS) -o $(PROG)

wells.o: wells.c wells.h

clean:
	rm -f $(PROG) $(OBJ)

tar:
	tar -cvzf wells.tgz `hg st -c | awk '{print $$2}'` .hg


