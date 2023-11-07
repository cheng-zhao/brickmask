CC = gcc
MPICC = mpicc
CFLAGS = -std=c99 -O3 -Wall
LIBS = -lm -lcfitsio
INCL = -Isrc -Ilib -Iio
# Set "USE_MPI = T" to enable MPI parallelisation
USE_MPI = T
# Uncomment the following line for eBOSS ELG masks
#CFLAGS += -DEBOSS -DFAST_FITS_IMG

# Settings for CFITSIO
CFITSIO_DIR = 
ifneq ($(CFITSIO_DIR),)
  LIBS += -L$(CFITSIO_DIR)/lib
  INCL += -I$(CFITSIO_DIR)/include
endif

# Settings for MPI
ifeq ($(USE_MPI), T)
  TARGET=BRICKMASK_MPI
  CFLAGS += -DMPI -Wno-stringop-overflow
else
  TARGET=BRICKMASK_NOMPI
endif

SRCS = $(wildcard src/*.c lib/*.c io/*.c)
EXEC = BRICKMASK

all: $(TARGET)

BRICKMASK_NOMPI:
	$(CC) $(CFLAGS) -o $(EXEC) $(SRCS) $(LIBS) $(INCL)

BRICKMASK_MPI:
	$(MPICC) $(CFLAGS) -o $(EXEC) $(SRCS) $(LIBS) $(INCL)

clean:
	rm $(EXEC)
