# Based on http://make.mad-scientist.net/papers/advanced-auto-dependency-generation/
# https://gist.github.com/maxtruxa/4b3929e118914ccef057f8a05c614b0f
# https://spin.atomicobject.com/2016/08/26/makefile-c-projects/

# ----------------------------------------
# Project definitions
# ----------------------------------------

# Folders
DDIR = .deps
SDIR = src

# Recursive wildcard
rwildcard=$(foreach d,$(wildcard $(1:=/*)),$(call rwildcard,$d,$2) $(filter $(subst *,%,$2),$d))

# Source files
SRCS = $(call rwildcard,$(SDIR),*.c)
#SRCS = $(wildcard $(SDIR)/**/*.c $(SDIR)/**/*.cpp)

# Target files
TRGS = $(patsubst %,%,$(basename $(notdir $(SRCS))))

# Dependency files
DEPS = $(addprefix $(DDIR)/,$(patsubst %,%.d,$(TRGS)))

# Gitignore file
GITIGNORE = .gitignore

# ----------------------------------------
# Compiler and linker definitions
# ----------------------------------------

# Compiler and linker
CC = gcc
CXX = g++

# Libraries
INCLUDES = 

# Flags for compiler
CFLAGS = -W -Wall -Wextra -pedantic -O2 -std=c11
CXXFLAGS = -W -Wall -Wextra -pedantic -O2 -std=c++11
DEPFLAGS = -MT $@ -MMD -MP -MF $(DDIR)/$*.Td

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

all: $(TRGS)

%: $(SDIR)/**/%.c
%: $(SDIR)/**/%.c $(DDIR)/%.d | $(DDIR)
	@ echo "${GREEN}Building target: ${BOLD}$@${GREEN}, using dependencies: ${BOLD}$^${NORMAL}"
	$(CC) $(CFLAGS) $(DEPFLAGS) $(filter %.c %.s %.o,$^) -o $@ $(INCLUDES)
	mv -f $(DDIR)/$*.Td $(DDIR)/$*.d && touch $@
	@if [ $(shell grep -Fxq "$@" ${GITIGNORE}; echo $$?) = 1 ] ; then \
		echo "${GREEN}Adding target ${BOLD}$@${GREEN} to ${GITIGNORE}${NORMAL}" ; \
		echo "$@" >> ${GITIGNORE} ; \
	fi

%: $(SDIR)/**/%.cpp
%: $(SDIR)/**/%.cpp $(DDIR)/%.d | $(DDIR)
	@ echo "${GREEN}Building target: ${BOLD}$@${GREEN}, using dependencies: ${BOLD}$^${NORMAL}"
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) $(filter %.cpp %.s %.o,$^) -o $@ $(INCLUDES)
	mv -f $(DDIR)/$*.Td $(DDIR)/$*.d && touch $@
	@if [ $(shell grep -Fxq "$@" ${GITIGNORE}; echo $$?) = 1 ] ; then \
		echo "${GREEN}Adding target ${BOLD}$@${GREEN} to ${GITIGNORE}${NORMAL}" ; \
		echo "$@" >> ${GITIGNORE} ; \
	fi

$(DDIR)/%.d: ;
.PRECIOUS: $(DDIR)/%.d

-include $(DEPS)

# ----------------------------------------
# Script rules
# ----------------------------------------

$(DDIR):
	@ echo "${GREEN}Creating directory: ${BOLD}$@${NORMAL}"
	mkdir -p $@

clean:
	@ echo "${GREEN}Cleaning up${NORMAL}"
	rm -rf $(DDIR)/ $(TRGS)

remade: clean all

.PHONY: all clean remade

# ----------------------------------------