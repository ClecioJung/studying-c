# Based on http://make.mad-scientist.net/papers/advanced-auto-dependency-generation/
# https://gist.github.com/maxtruxa/4b3929e118914ccef057f8a05c614b0f
# https://spin.atomicobject.com/2016/08/26/makefile-c-projects/

# Recursive wildcard
rwildcard=$(foreach d,$(wildcard $(1:=/*)),$(call rwildcard,$d,$2) $(filter $(subst *,%,$2),$d))

# ----------------------------------------
# Project definitions
# ----------------------------------------

# Folders
BDIR = bin
DDIR = .deps
SDIR = src
LIB_ODIR = .obj
LIB_DDIR = .deps/lib
LIB_DIR = lib

# Statically linked library
LIB_SLL = $(LIB_ODIR)/libmath.a

# Source files
SRCS = $(call rwildcard,$(SDIR),*.c)
# Dependency files
DEPS = $(addprefix $(DDIR)/,$(patsubst %,%.d,$(patsubst %,%,$(basename $(notdir $(SRCS))))))
# Executable target files
BINS = $(addprefix $(BDIR)/,$(patsubst %,%,$(basename $(notdir $(SRCS)))))
# Library source files
LIB_SRCS = $(wildcard $(LIB_DIR)/*.c)
# Library dependency files
LIB_DEPS = $(addprefix $(LIB_DDIR)/,$(patsubst %,%.d,$(patsubst %,%,$(basename $(notdir $(LIB_SRCS))))))
# Library object files
LIB_OBJS = $(addprefix $(LIB_ODIR)/,$(patsubst %,%.o,$(patsubst %,%,$(basename $(notdir $(LIB_SRCS))))))

# ----------------------------------------
# Compiler and linker definitions
# ----------------------------------------

# Compiler and linker
CC = gcc
AR = ar

# Libraries
INCLUDES = $(LIB_SLL)

# Flags for compiler
CFLAGS = -W -Wall -Wextra -pedantic -Werror -O2 -std=c11

# ----------------------------------------
# Fomating macros
# ----------------------------------------

BOLD = \033[1m
NORMAL = \033[0m
RED = \033[0;31m
GREEN = \033[0;32m

# ----------------------------------------
# Compilation and linking rules
# ----------------------------------------

all: $(LIB_SLL) $(BINS)

$(LIB_SLL): $(LIB_OBJS) | $(LIB_ODIR)
	@ echo "${GREEN}Building statically linked library: ${BOLD}$@${GREEN}, using dependencies: ${BOLD}$^${NORMAL}"
	$(AR) rcs $(LIB_SLL) $(LIB_OBJS)

$(LIB_ODIR)/%.o : $(LIB_DIR)/%.c
$(LIB_ODIR)/%.o : $(LIB_DIR)/%.c $(LIB_DDIR)/%.d | $(LIB_DDIR) $(LIB_ODIR)
	@ echo "${GREEN}Building target: ${BOLD}$@${GREEN}, using dependencies: ${BOLD}$^${NORMAL}"
	$(CC) $(CFLAGS) -MT $@ -MMD -MP -MF $(LIB_DDIR)/$*.Td -c $(filter %.c %.s %.o,$^) -o $@
	mv -f $(LIB_DDIR)/$*.Td $(LIB_DDIR)/$*.d && touch $@

$(BDIR)/%: $(SDIR)/**/%.c
$(BDIR)/%: $(SDIR)/**/%.c $(LIB_SLL) $(DDIR)/%.d | $(DDIR) $(BDIR)
	@ echo "${GREEN}Building target: ${BOLD}$@${GREEN}, using dependencies: ${BOLD}$^${NORMAL}"
	$(CC) $(CFLAGS) -MT $@ -MMD -MP -MF $(DDIR)/$*.Td $(filter %.c %.s %.o,$^) -o $@ $(INCLUDES)
	mv -f $(DDIR)/$*.Td $(DDIR)/$*.d && touch $@

$(DDIR)/%.d: ;
.PRECIOUS: $(DDIR)/%.d

$(LIB_DDIR)/%.d: ;
.PRECIOUS: $(LIB_DDIR)/%.d

-include $(DEPS)
-include $(LIB_DEPS)

# ----------------------------------------
# Script rules
# ----------------------------------------

test: $(BINS)
	@ for t in $(BINS); do \
		echo "${GREEN}Running test: ${BOLD}$$t${NORMAL}" ; \
		./$$t 1>/dev/null ; \
		if [ $$? -ne 0 ] ; then \
			echo "${RED}Test ${BOLD}$$t${RED} failed${NORMAL}" ; \
			exit 1 ; \
		fi ; \
	done
	@echo "${GREEN}Success, all tests passed.${NORMAL}"

$(BDIR) $(DDIR) $(LIB_ODIR) $(LIB_DDIR):
	@ echo "${GREEN}Creating directory: ${BOLD}$@${NORMAL}"
	mkdir -p $@

clean:
	@ echo "${GREEN}Cleaning up${NORMAL}"
	rm -rf $(BDIR) $(DDIR) $(LIB_ODIR) $(LIB_DDIR) *.d *.o *.a *.so

remade: clean all

.PHONY: all test clean remade

# ----------------------------------------