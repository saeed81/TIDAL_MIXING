#
include make.macro
#
LIB_SRC = src/par_tide.f90 \
          src/arraydef.f90 \
          src/in_manager.f90 \
	  src/phycst.f90 \
	  src/bessel_jy0.f90 \
	  src/bessel_jy1.f90 \
	  src/grn_tide.f90\
          src/lib_otps.f90\
	  src/tmx_array.f90\
	  src/io_ezcdf.f90\
	  src/arayini.f90\
	  src/domini.f90 \
	  src/wrtout.f90\
          src/soltmx.f90\
          src/tmx.f90
#
LIB_OBJ = $(LIB_SRC:.f90=.o)
#
LIBTMX = lib/libtmx.a
#
LIB = -L./lib -ltmx

#
EXEC_TMX = ./tmx
#

.SUFFIXES:
.SUFFIXES: .f90 .o
#
all:    dirmod $(EXEC_TMX) #rmobj rmmod 
	@echo "                 "
	@echo "                 "
	@echo "                 "
	#@echo \(.o .mod removed safely\) 
	@echo "                 "
	@echo "                 " 
	@echo   TMX model is OK
	@echo "                 "
	@echo "                 " 
#
dirmod:
	@if [ ! -d ./mod ]; then mkdir ./mod ; fi
# 

#$(EXEC_TMX) :   src/model.f90 $(LIBTMX) 
#	$(FC) $(FF) -o $(EXEC_TMX) src/model.f90 $(LIB) $(LIB_CDF)
#
$(EXEC_TMX) : $(LIBTMX) model.o
	@echo "                          "
	@echo "                          " 
	$(FC) $(FF) -o $(EXEC_TMX) model.o $(LIBTMX) $(LIB_CDF)

#
# main program
model.o : src/model.f90
	$(FC) $(FF) -c src/model.f90
#
$(LIBTMX): $(LIB_OBJ)
	@echo ""
	@echo $(LIB_OBJ)
	@echo ""
	@mkdir -p lib
	ar -rv lib/libtmx.a  $(LIB_OBJ)
	ranlib lib/libtmx.a
	@echo ""
#
#
.f90.o: $(LIB_SRC)
	$(FC) -c $(FF) $< -o $*.o
#

#
rmobj :
	@echo "                "
	@echo "                "
	@\rm -f src/*.o model.o 

rmmod  : 
	@echo "                "
	@echo "                "
	@\rm -f mod/*.mod

.PHONY: clean

clean:
	rm -f *.x *.out *.err *.log *~ out src/*.o mod/* lib/* *.nc src/*.nc src/*.x src/datanc/*.nc src/*.err src/*.out output/* 
#
install:
	@cp -f *.x src/ && cd src && ./*.x
#
unistall:
	rm -rf bin
#