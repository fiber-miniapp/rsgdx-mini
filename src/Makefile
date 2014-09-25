
FC = mpif90
FFLAGS = -O3 -openmp -parallel -xHost
LDFLAGS = $(BLAS)

BLAS = -L/usr/local/intel/composer_xe_2013/mkl/lib/intel64 \
-Wl,--start-group \
-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core \
-Wl,--end-group

OBJS = \
	m_precisions.o \
	m_comm.o \
	m_constant.o \
	m_conv.o \
	m_gmt.o \
	m_hmatrices.o \
	m_hmatrices_rkmatrixxvector_index.o \
	m_outputval.o \
	mpi_comm_timer.o \
	timer.o \
	visualizer.o \
	convolution.o \
	coseismic_slip.o \
	RSGDX_MO_multi.o

.SUFFIXES: .o .f90 .c

.f90.o: precision_pri.inc
	$(FC)  $(FFLAGS) -o $@ -c $<

%.o: %.mod

all: rsgdx.ex

rsgdx.ex: $(OBJS)
	$(FC)  $(FFLAGS) $(LDFLAGS) -o $@ $(OBJS)

clean:
	rm -f *.o *.mod rsgdx.ex

