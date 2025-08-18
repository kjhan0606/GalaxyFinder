#
# Makefile for PM/Tree code pyul.exe - April 23, 2001 - JJD
#

#--- Executable

PROGS =  opfof.exe 

#--- Objects to build

OBJS		= opfof.o Time.o Treewalk.fof.ordered.o  read.o  header.o ramses2read.o

OS = LINUX_GCC



#--- Targets

SUBDIRS = 

all: 
	for subdir in $(SUBDIRS); do \
		(cd $$subdir && $(MAKE) all);\
	done 
	$(MAKE) this

this: $(OBJS)
	rm -f $(PROGS) 
	$(MAKE) $(PROGS)
new: 
	$(MAKE) clean
	$(MAKE) all

clean: 
	rm -f $(OBJS)
	for subdir in $(SUBDIRS); do \
		(cd $$subdir && $(MAKE) clean);\
	done 


#mkpower:
#	rm -f $(MKPOW) 
#	rm -f $(OBJPOW) 
#	$(MAKE) $(MKPOW)

default: all

#$(MKPOW): $(OBJPOW)



include Rules.make
