#-----------------------------------------------------------------------
# GNUmakefile
# May-2005	Harrison B. Prosper
# 19-Oct-2005   HBP Add rules to make Root shared library
# 02-Mar-2011   HBP Convert to Root
#-----------------------------------------------------------------------
program:= blimit
#-----------------------------------------------------------------------
ifndef ROOTSYS
$(error Environment variable ROOTSYS undefined)
endif
ifndef verbose
AT	:=@
endif
#-----------------------------------------------------------------------
C++		:= g++
LD		:= g++
#-----------------------------------------------------------------------
SYSLIBS		:= -lm -lc
ROOTFLAGS	:= $(shell root-config --cflags)
ROOTLIBS	:= $(shell root-config --libs) $(SYSLIBS)
CPPFLAGS	:= -I. -D__MAIN__ $(ROOTFLAGS)
CXXFLAGS	:= -c -ggdb -O2
LDFLAGS		:= -ggdb -O2
LDSHARE		:= -rdynamic -shared -Wl,-E -Wl,-soname,$(program).
#----------------------------------------------------------------------
source	:= blimit.cc
object	:= blimit.o
#-----------------------------------------------------------------------
bin:	$(program)

$(program)	: $(object)
	@echo "---> Linking $@"
	$(AT)$(LD) $(LDFLAGS) $(object) $(ROOTLIBS) $(LIBS) -o $@
	@echo ""

$(object)	: $(source)
	@echo "---> Compiling $<"
	$(AT)$(C++) $(CPPFLAGS) $(CXXFLAGS) -o $@ $< 

clean	:
		@rm -rf	*.o $(program)


