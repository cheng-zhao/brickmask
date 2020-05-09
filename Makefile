CC = gcc
LIBS = -lm -lcfitsio
INCL = 
OPTS = -O3 -Wall $(INCL) $(LIBS) #-DOMP -fopenmp
SRCS = baselib.c brickmask.c find_brick.c load_conf.c read_data.c save_res.c
EXEC = brickmask

all:
	$(CC) $(SRCS) -o $(EXEC) $(OPTS)

clean:
	rm -f $(EXEC)

