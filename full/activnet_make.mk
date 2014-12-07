###########################################################################
## Makefile generated for MATLAB file/project 'activnet_lib'. 
## 
## Makefile     : activnet_make.mk
## Generated on : Tue Oct 14 17:24:54 2014
## MATLAB Coder version: 2.7 (R2014b)
## 
## Build Info:
## 
## Final product: $(RELATIVE_PATH_TO_ANCHOR)/activnet
## Product type : executable
## 
###########################################################################

###########################################################################
## MACROS
###########################################################################

# Macro Descriptions:
# PRODUCT_NAME            Name of the system to build
# MAKEFILE                Name of this makefile
# COMPUTER                Computer type. See the MATLAB "computer" command.

PRODUCT_NAME              = activnet_gen
MAKEFILE                  = activnet_make.mk
START_DIR                 = ./activnet_gen
RELATIVE_PATH_TO_ANCHOR   = .


#------------------------
# BUILD TOOL COMMANDS
#------------------------

# C Compiler: 
CC = gcc

# Linker: 
LD = gcc


#----------------------------------------
# code optimization in gcc Build Configuration
#----------------------------------------

CFLAGS               = -c -O3


###########################################################################
## OUTPUT INFO
###########################################################################

PRODUCT = $(RELATIVE_PATH_TO_ANCHOR)/activnet
PRODUCT_TYPE = "executable"
BUILD_TYPE = "Executable"

###########################################################################
## INCLUDE PATHS
###########################################################################


INCLUDES = -I$(START_DIR)


###########################################################################
## OBJECTS
###########################################################################

OBJS = activnet_gen.o rt_nonfinite.o rtGetNaN.o rtGetInf.o main.o

###########################################################################
## SYSTEM LIBRARIES
###########################################################################

SYSTEM_LIBS = -lm


#---------------
# C Compiler
#---------------


CFLAGS += $(INCLUDES)




all : $(PRODUCT)
	@echo "### Successfully generated all binary outputs."



###########################################################################
## FINAL TARGET
###########################################################################

#-------------------------------------------
# Create a standalone executable            
#-------------------------------------------

$(PRODUCT) : $(OBJS) 
	$(LD) $(LDFLAGS) -o $(PRODUCT) $(OBJS) $(SYSTEM_LIBS)
	@rm -f $(OBJS)

###########################################################################
## INTERMEDIATE TARGETS
###########################################################################

#---------------------
# SOURCE-TO-OBJECT
#---------------------

%.o : %.c
	$(CC) $(CFLAGS) -o "$@" "$<"


%.o : $(START_DIR)/%.c
	$(CC) $(CFLAGS) -o "$@" "$<"



###########################################################################
## MISCELLANEOUS TARGETS
###########################################################################

info : 
	@echo "### PRODUCT = $(PRODUCT)"
	@echo "### PRODUCT_TYPE = $(PRODUCT_TYPE)"
	@echo "### BUILD_TYPE = $(BUILD_TYPE)"
	@echo "### INCLUDES = $(INCLUDES)"
	@echo "### OBJS = $(OBJS)"
	@echo "### LIBS = $(LIBS)"
	@echo "### SYSTEM_LIBS = $(SYSTEM_LIBS)"
	@echo "### CFLAGS = $(CFLAGS)"
	

clean : 
	@rm -f $(PRODUCT)
	@rm -f $(OBJS)
	@echo "### Deleted all derived files."


