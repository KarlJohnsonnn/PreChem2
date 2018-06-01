IN1=M_DEF_Mac

include $(IN1) 

OBJo = -L$(LIB_O) -lchemieo
OBJg = -L$(LIB_D) -lchemieg

all: PreChem2 PreChem2_dbg
       
PreChem2:  optimized  
	$(F90) $(OPT_O) $(KINCL_O) -o PreChem2.exe \
        PreChem2.f90 $(OP) $(OBJo) $(CL);

PreChem2_dbg:  debug
	$(F90) $(OPT_D) $(KINCL_D) -o PreChem2_dbg.exe \
        PreChem2.f90 $(OP) $(OBJg) $(CL);

optimized: 
	@make -f Make_src "IN2=$(IN1)" "LIB2=$(LIB_O)" "OPT2=$(OPT_O)" "CHEM=chemieo" "KINCL=$(KINCL_O)"

debug: 
	@make -f Make_src "IN2=$(IN1)" "LIB2=$(LIB_D)" "OPT2=$(OPT_D)" "CHEM=chemieg" "KINCL=$(KINCL_D)"

clean:
	rm -f *.o *.mod LIB*/*

distclean:
	rm -f *.exe

