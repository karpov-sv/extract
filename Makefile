# Makefile for standalone tools and FAVOR library (for external usage, static and shared)

CC       = gcc $(PROFILE) -std=gnu89
PROFILE  = -pg
DEBUG    = -DDEBUG
GDEBUG   = -g
OPTIMIZE = -O2 -msse4 -ffast-math
INCLUDES = -I.
THREADED = -D_REENTRANT -D_GNU_SOURCE -DCORO_UCONTEXT -D_XOPEN_SOURCE
FLAGS    = -Wall -fPIC -fms-extensions -fno-strict-aliasing -fopenmp
CFLAGS   = $(FLAGS) $(THREADED) $(GDEBUG) $(DEBUG) $(OPTIMIZE) $(INCLUDES)
LDLIBS  = -L. -lm -lcfitsio -lpthread -lgomp

ifeq ($(shell uname), Linux)
LDLIBS  += -ltcmalloc -lcrypt
endif

# GSL part
LDLIBS += -lgsl -lgslcblas

# GLib part
CFLAGS += `pkg-config --cflags glib-2.0`
LDLIBS += `pkg-config --libs glib-2.0`

TARGETS = \
	extract \

LIBFAVOR_OBJS = \
	utils.o \
	time_str.o \
	coords.o \
	image.o \
		image_keywords.o \
		image_fits.o \
		image_clean.o \
	kdtree.o \
	psf.o \
	mpfit.o \

EXTRACT_OBJS = \
	extract.o \
		extract_peaks.o \

all: depend $(TARGETS)

$(TARGETS): $(LIBFAVOR_OBJS)

extract: $(EXTRACT_OBJS)

clean:
	rm -f *~ *.o */*.o $(TARGETS) $(SCHEDULER_OBJS)

depend:
	@echo "# DO NOT DELETE THIS LINE -- make depend depends on it." >Makefile.depend
	@makedepend *.c -fMakefile.depend -I/usr/local/include 2>/dev/null

-include Makefile.depend # DO NOT DELETE
# DO NOT DELETE
