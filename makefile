# Makefile written by Arijit on 24/05/17

FORT=gfortran
CC=g++

MPIF90=mpif90
MPICC=mpicxx

FORTFLAGS=-cpp -ffree-line-length-0
LIBS= -lfftw3 #-lm -llapack -lblas  

HEADERS_PATH=./src -I/usr/include
HEADERS=

SRC_PATH=./src
EXE_FILES_SRC_PATH=./src
_SOURCES= const.f90 parm.f90  fft.f90 fft_module.f90 integral.f90 \
	 imaginary_time_frequency_FFTs.f90 green.f90 dmft.f90 ipt.f90  eq-dmft_HM.f90 # all the modules which are arranged in such a way the most dependent goes right and less to left

OBJ=$(_SOURCES:.f90=.o) # create object file list

EXE_SOURCES= main_HM_eq.f90 # fft_test4.f90 fft_test.f90 initialise.f90 \
	#fft_test2.f90  fft_test3.f90  fft_test4.f90  fft_test.f90  # main program file name list name list


EXES=$(EXE_SOURCES:.f90=.x) # making exectable file name list

all: $(OBJ) $(EXES) clean_obj

# it will recursively do the .x file form string EXE wich is made from EXE_SOURCES:.f90
$(EXES):%.x: $(EXE_FILES_SRC_PATH)/%.f90 $(OBJ) 
	$(FORT) $(FORTFLAGS) -I$(HEADERS_PATH) $< -o $@ $(OBJ) $(LIBS)

$(OBJ):%.o: $(SRC_PATH)/%.f90 #$(CHEADERS_PATH)/$(CHEADERS)
	$(FORT) $(FORTFLAGS) -I$(HEADERS_PATH) -c $< -o $@

clean_obj: 
	rm *.o *.mod
clean:
	rm *.x 

