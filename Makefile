FCMPI=mpiifort
FC=ifort

#FFLAGS += -fno-inline -no-ip -fp-stack-check -fpe0 -ftz -O3 
FFLAGS += -ftrapuv -fno-inline -no-ip -check all -check noarg_temp_created -fp-stack-check -fpe0 -ftz -O0 
#FFLAGS += -qopenmp

LIB=

BETA=../exe/beta.exe
TARGETS= $(BETA)

all: $(TARGETS)

MOD= \
      kinds.mod \
      jp_pkind.mod \
      jp_pkind2.mod \
      jp_pietc.mod \
      jp_pietc_s.mod \
      jp_pmat.mod \
      jp_pmat4.mod \
      jp_pbfil.mod \
      jp_pbfil2.mod \
      jp_pbfil3.mod \
      mg_input.mod \
      mg_output.mod \
      mg_mppstuff.mod \
      mg_domain.mod \
      mg_parameter.mod \
      mg_bocos.mod \
      mg_transfer.mod \
      mg_generations.mod \
      mg_interpolate.mod \
      mg_filtering.mod \
      mg_timers.mod \
      mg_entrymod.mod \
      mg_intstate.mod 

OBJ= \
      kinds.o \
      jp_pkind.o \
      jp_pkind2.o \
      jp_pietc.o \
      jp_pietc_s.o \
      jp_pmat.o \
      jp_pmat4.o \
      jp_pbfil.o \
      jp_pbfil2.o \
      jp_pbfil3.o \
      mg_input.o \
      mg_output.o \
      mg_mppstuff.o \
      mg_domain.o \
      mg_parameter.o \
      mg_bocos.o \
      mg_transfer.o \
      mg_generations.o \
      mg_interpolate.o \
      mg_filtering.o \
      mg_timers.o \
      mg_entrymod.o \
      mg_intstate.o \
      RBETA_TEST.o


RM= /bin/rm -f

$(TARGETS): $(OBJ)
	$(FCMPI) $(FFLAGS) $(OBJ) -I ./ -o $(BETA)

%.o:%.f90
	$(FCMPI) $(FFLAGS) -c $<

clean:
	$(RM) $(OBJ) $(MOD) $(BETA)

# dependencies

RBETA_TEST.o: kinds.o mg_entrymod.o mg_mppstuff.o mg_filtering.o mg_timers.o \
              mg_transfer.o mg_parameter.o
jp_pbfil.o: kinds.o jp_pietc.o jp_pmat.o mg_parameter.o
jp_pbfil2.o: jp_pkind.o
jp_pbfil3.o: jp_pkind.o jp_pkind2.o jp_pietc.o jp_pmat.o jp_pmat4.o jp_pbfil2.o
jp_pietc.o: jp_pkind.o
jp_pietc_s.o: jp_pkind.o
jp_pmat.o: jp_pkind.o jp_pietc.o
jp_pmat4.o: jp_pkind.o jp_pietc.o jp_pietc_s.o
mg_domain.o: kinds.o mg_parameter.o mg_mppstuff.o
mg_mppstuff.o: kinds.o mg_parameter.o
mg_parameter.o: kinds.o jp_pietc.o
mg_bocos.o: kinds.o mg_mppstuff.o mg_parameter.o mg_domain.o
mg_transfer.o: kinds.o mg_parameter.o mg_intstate.o  \
               mg_interpolate.o mg_bocos.o mg_mppstuff.o mg_timers.o
mg_generations.o: kinds.o mg_parameter.o mg_mppstuff.o mg_bocos.o mg_intstate.o \
                   mg_timers.o mg_output.o
mg_intstate.o: kinds.o jp_pkind2.o mg_parameter.o mg_mppstuff.o jp_pbfil.o
mg_input.o: kinds.o mg_parameter.o mg_mppstuff.o
mg_output.o: kinds.o mg_parameter.o mg_mppstuff.o
mg_interpolate.o: kinds.o mg_parameter.o mg_intstate.o mg_mppstuff.o
mg_filtering.o: kinds.o mg_timers.o mg_parameter.o mg_intstate.o mg_generations.o \
                mg_bocos.o jp_pbfil.o mg_mppstuff.o mg_output.o
mg_timers.o: kinds.o mg_mppstuff.o
mg_entrymod.o: kinds.o mg_parameter.o mg_mppstuff.o mg_domain.o mg_intstate.o \
               mg_intstate.o mg_interpolate.o mg_input.o mg_output.o


