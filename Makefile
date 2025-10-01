# Makefile compiles C/C++ code with several optimisations
#
# 1. Configuration option allows separate debug, release and profile builds 
#  
#  DEBUGGING mode
# 
#	 -g: Produces debugging information in the OS's native format (stabs, COFF, XCOFF, or DWARF 2).
#	 -ggdb: Produces debugging information specifically intended for gdb (defaults to level2).
#	 -ggdb3: Produces extra debugging information, for example: including macro definitions.
#	 -Wall: Enables a base set of warnings, but not all.
#	 -Wextra: Enables additional warnings, such as missing function prototypes.
#	 -Wpedantic: Enforces strict standard compliance.
#	 -Wshadow: Warns when a variable declaration shadows another (prevents subtle bugs).
#	 -Wconversion / -Wsign-conversion: Warns about implicit type conversions that may change values.
#	 -Wnull-dereference: Detects null pointer dereferences.
#	 -Wformat=2: Stricter checks for printf/scanf format strings.
#	 -Wduplicated-cond / -Wduplicated-branches: Detects duplicate conditions or branches in if/switch statements.
#	 -Wlogical-op: Warns about suspicious use of logical operators.
#	 -Wuseless-cast: Warns about redundant type casts.
#
# 	RELEASE mode 
#
#	- Optimization & Performance Flags:
#
#		-O3: Enables the highest level of optimization.
#		    Includes function inlining, loop unrolling, and vectorization.
#		    Can lead to significant performance improvements but might increase binary size.
#		    If binary size is a concern, consider -O2 instead.
#			 Is possibly volatile - possible to break working code if undefined behaviour invoked
#		-ffast-math: Enables aggressive floating-point optimizations.
#		    Assumes no special cases (e.g., NaNs, infinities).
#		    Speeds up floating-point operations but can lead to non-standard behavior.
#		    Use cautiously in scientific computing where exact IEEE compliance is required.
#		-march=native: Optimizes the generated code for the host CPU.
#		    Allows the compiler to use CPU-specific instructions and optimizations.
#		    Improves performance but makes binaries less portable.
#		    If distributing the binary, consider -mtune=generic instead.
#		-flto: Enables Link-Time Optimization (LTO).
#		    Optimizes across compilation units, reducing binary size and improving execution speed.
#		    Can inline functions across files for better performance.
#		    Also added to LDFLAGS for full effect.
#
#	- Debugging & Error Handling:
#
#		-DNDEBUG: Disables assert() checks. In release builds, assertions are typically not needed,
# 			 improving performance. If you still need debugging capabilities, consider defining NDEBUG selectively.
#
#	- Warning Flags (catch code quality if minor changes to release code):
#
#		-Wall: Enables a basic set of useful warnings.
#		-Wextra: Enables additional warnings not included in -Wall, such as missing function prototypes.
#		-Wpedantic: Enforces strict compliance with the C/C++ standard.
#		    Helps catch non-standard extensions and potential portability issues.
#
#	- Math-Specific Optimization:
#
#		-fno-math-errno: Disables setting errno after math functions like sqrt() and log().
#		    Saves CPU cycles but removes error reporting via errno.
#		    Safe for most high-performance applications that don’t rely on errno.
#
#	- Linker Flags:
#
#	    LDFLAGS += -flto: Ensures Link-Time Optimization (LTO) is also applied during linking.
#
#	PROFILE mode
#	
#	-O2: Agressive optimisation not allowed so reduced to level 2
#	-pg: Add profiling flags to functions for gprof profiling
#	-fno-omit-frame-pointer: Improve frame tracking
#	
# 2. Depending on code structure makefile will find and compile for multiple source files
# 	 Will find all files in source dir with same extension
#	 To extend to both c and cpp files:
#
#	 - Find both source extensions and generate obj files:
#
#		SRC_EXTS=c cpp
#		SOURCES=$(foreach ext, $(SRC_EXTS), $(wildcard $(SRC_DIR)/*.$(ext)))
#		OBJECTS=$(patsubst $(SRC_DIR)/%,$(OBJ_DIR)/%, $(SOURCES:.c=.o) $(SOURCES:.cpp=.o))
#
#	 - Add compilation rules for each type - note linking same for both
#
#		(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
#			(CC) -MMD -MP -MF $(DEPS_DIR)/$*.d -c $(CFLAGS) $< -o $@
#
#		(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
#			(CXX) -MMD -MP -MF $(DEPS_DIR)/$*.d -c $(CFLAGS) $< -o $@
#
# 3. Project structure
#	- Auto generates object dir (obj) and dependencies dir (dep)
#		-p suppresses error message for existing dirs
# 	- Option to also generate Doxygen docs and doxygen config files
#	
# 4. Auto generate a dependencies (.d) file which is a makefile it will tell make whether 
#    headers have been edited and which object files need to be recompiled
# 
#	 -MMD: Produces dependencies file - not including system req (-MD to include)
#	 -MP: Makes phony targets for all dependencies
#	 -MF: Overides default dep file - so can place in directory

# Automatically detect all source files in SRC_DIR
SRC_DIR=.
SRC_EXT=c
SOURCES=$(wildcard $(SRC_DIR)/*.$(SRC_EXT))
EXECUTABLE=MoT-Voellmy
DATE=2025-05-20 # $(shell date +%Y-%m-%d)

# Compiler and Flags
CXX=gcc

# Directories
OBJ_DIR=obj
DEPS_DIR=deps
DOCS_DIR=docs
DOC_SRC=doc_src
DOXYFILE=$(DOC_SRC)/Doxyfile

# Configuration Handling
LDFLAGS += -lm

CONF?=release
ifeq ($(CONF),release)
#    CFLAGS += -O3 -ffast-math -DNDEBUG -Wall -Wextra -Wpedantic -march=native -flto=auto -fno-math-errno
    CFLAGS += -O2 -DNDEBUG -Wall -Wextra -Wpedantic
else ifeq ($(CONF),debug)
    CFLAGS += -ggdb3 -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -Wnull-dereference \
           -Wformat=2 -Wduplicated-cond -Wduplicated-branches -Wlogical-op -Wuseless-cast
else ifeq ($(CONF),profile)
    CFLAGS += -pg -O2 -Wall -Wextra -Wpedantic -fno-omit-frame-pointer
    LDFLAGS += -pg
else
    $(error "CONF should be 'release','debug' or 'profile'")
endif

# Build for linux or windows
COMP?=linux
ifeq ($(COMP),linux)
	EXECUTABLE := $(EXECUTABLE)-linux.$(DATE)
else ifeq ($(COMP),static)
	EXECUTABLE := $(EXECUTABLE)-linux-static.$(DATE).exe
	CFLAGS += -static
else ifeq ($(COMP),windows)
	CXX=x86_64-w64-mingw32-gcc
	EXECUTABLE := $(EXECUTABLE)-win64.$(DATE).exe
else ifeq ($(COMP),macos)
	CXX=clang
	EXECUTABLE := $(EXECUTABLE)-macOS.$(DATE)
	CFLAGS += -arch arm64 -arch x86_64
endif

# Generate object files list dynamically
OBJECTS=$(patsubst $(SRC_DIR)/%.$(SRC_EXT), $(OBJ_DIR)/%.o, $(SOURCES))

# Default target
all: gendirs $(EXECUTABLE)

# Ensure necessary directories exist
gendirs:
	@mkdir -p $(OBJ_DIR) $(DEPS_DIR)

# Linking
$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(OBJECTS) $(LDFLAGS) -o $@

# Compilation rule
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.$(SRC_EXT)
	$(CXX) -MMD -MP -MF $(DEPS_DIR)/$*.d -c $(CFLAGS) $< -o $@

# Phony targets
.PHONY: run clean rebuild print debug help verbose docs doxyconfig cleandocs
debug:
	@echo "Sources: $(SOURCES)"
	@echo "Objects: $(OBJECTS)"

run: all
	./$(EXECUTABLE)

clean:
	rm -f $(EXECUTABLE) $(OBJECTS) $(wildcard $(DEPS_DIR)/*.d)

rebuild: clean all

print:
	@echo $(SOURCES)
	
help:
	@echo "Available targets:"
	@echo "  make all       - Compile the project"
	@echo "  make run       - Build and run the executable"
	@echo "  make debug     - Print source and object file lists"
	@echo "  make rebuild   - Clean and rebuild"
	@echo "  make docs		- Build documentation with doxygen"
	@echo "  make clean     - Remove compiled files"
	@echo "  make cleandocs - Remove documentation files (not configuration)"
	@echo "Compilation modes:"
	@echo "  make CONF=release	- Default build"
	@echo "  make CONF=debug	- Build with debug flags for gdb"
	@echo "  make CONF=release	- Build with profile flags for gprof"
	@echo "Build modes:"
	@echo "  make COMP=linux	- Default build"
	@echo "  make COMP=windows	- Build for Windows using mingw32"
	
verbose: 
	$(MAKE) V=1

# Construct and edit doxygen source files 
# To view html: docs/html/index.html
# To view pdf: docs/latex/refman.pdf
docs: doxyconfig
	@doxygen $(DOXYFILE)
	@if [ -f $(DOCS_DIR)/latex/Makefile ]; then $(MAKE) -C $(DOCS_DIR)/latex; fi

doxyconfig:
	@mkdir -p $(DOCS_DIR) $(DOC_SRC)
	@if [ ! -f $(DOXYFILE) ]; then doxygen -g $(DOXYFILE); fi
	@case "$$(uname)" in \
	  Darwin*) SED_INPLACE="sed -i ''" ;; \
	  *)       SED_INPLACE="sed -i" ;; \
	esac; \
	$$SED_INPLACE 's|^OUTPUT_DIRECTORY .*|OUTPUT_DIRECTORY = $(DOCS_DIR)|' $(DOXYFILE); \
	$$SED_INPLACE 's|^INPUT .*|INPUT = $(SRC_DIR) $(DOC_SRC)|' $(DOXYFILE); \
	$$SED_INPLACE 's|^RECURSIVE .*|RECURSIVE = YES|' $(DOXYFILE); \
	$$SED_INPLACE 's|^PROJECT_NAME .*|PROJECT_NAME = \"$(EXECUTABLE)\"|' $(DOXYFILE); \
	$$SED_INPLACE 's|^GENERATE_HTML .*|GENERATE_HTML = YES|' $(DOXYFILE); \
	$$SED_INPLACE 's|^GENERATE_LATEX .*|GENERATE_LATEX = YES|' $(DOXYFILE); \
	$$SED_INPLACE 's|^GENERATE_MAN .*|GENERATE_MAN = NO|' $(DOXYFILE); \
	$$SED_INPLACE 's|^GENERATE_XML .*|GENERATE_XML = NO|' $(DOXYFILE); \
	$$SED_INPLACE 's|^GENERATE_RTF .*|GENERATE_RTF = NO|' $(DOXYFILE); \
	$$SED_INPLACE 's|^EXTRACT_ALL .*|EXTRACT_ALL = YES|' $(DOXYFILE); \
	$$SED_INPLACE 's|^HAVE_DOT .*|HAVE_DOT = YES|' $(DOXYFILE); \
	$$SED_INPLACE 's|^GENERATE_TREEVIEW .*|GENERATE_TREEVIEW = YES|' $(DOXYFILE); \
	$$SED_INPLACE 's|^CALL_GRAPH .*|CALL_GRAPH = YES|' $(DOXYFILE); \
	$$SED_INPLACE 's|^CALLER_GRAPH .*|CALLER_GRAPH = YES|' $(DOXYFILE)

cleandocs:
	rm -rf $(DOCS_DIR)/html/* $(DOCS_DIR)/latex/*
	
	
# Dependency handling
-include $(wildcard $(DEPS_DIR)/*.d)

