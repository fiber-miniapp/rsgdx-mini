
FC = mpifrtpx
#	FFLAGS = -Kfast,openmp
FFLAGS = -Kfast,parallel,openmp
LDFLAGS = -SCALAPACK -SSL2BLAMP

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

