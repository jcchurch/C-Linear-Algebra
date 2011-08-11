CC = /usr/bin/cc
RM  = /bin/rm
CFLAGS = -O2 -lm

LIBRARY = matrix.o L2_distance.o matrixadv.o qr.o eigen.o qsort.c svd.o

TEST_APS = qr_decomposition_test invtest eigen_test quicksort_test svd_test

all: $(TEST_APS)

%.o: $*.c $*.h
	$(CC) $(CFLAGS) -c $*.c

$(TEST_APS): $(LIBRARY)
	$(CC) $(CFLAGS) $(LIBRARY) $@.c -o $@

clean:
	-$(RM) *.o
	-$(RM) $(TEST_APS)
