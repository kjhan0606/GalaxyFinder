# In-memory stellar linker. It does not build opfof.exe.
CC ?= mpicc
CFLAGS = -std=gnu11 -O3 -fPIC -shared -DINTEL -DINDEX -DVarPM -DXYZDBL \
	-DWGROUPSIZE=8 -DNMEG=17000L -DINCLUDE_TREE_FORCE -D_LARGE_FILES \
	-DSAVESLICE -DPMSEEDFORCE -DQUADHILBERT -DNENER=0 -DNPRE=8 \
	-DREAD_SINK -DNCHEM=9 -DNDUST=4
CPPFLAGS = -I. -I..

.PHONY: all clean test

all: libopfof_sb.so

libopfof_sb.so: sb_adapter.c Treewalk.fof.ordered.c fof.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ sb_adapter.c Treewalk.fof.ordered.c -lm

clean:
	rm -f libopfof_sb.so
