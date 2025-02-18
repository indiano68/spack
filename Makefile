# Compiler settings
CXX        := g++
CXXFLAGS   := -Wall -Wextra -std=c++17 -Iinclude

# Archiver settings for static library
AR         := ar
ARFLAGS    := rcs

# Directories
SRCDIR     := src
TESTSRCDIR := test
BUILDDIR   := build
OBJDIR     := $(BUILDDIR)/obj
LIBDIR     := $(BUILDDIR)/lib
TESTBINDIR := $(BUILDDIR)/test

# Library output
LIBNAME    := mylib
LIBARCH    := lib$(LIBNAME).a
LIB        := $(LIBDIR)/$(LIBARCH)

# Library source files and corresponding object files
SRCFILES   := $(wildcard $(SRCDIR)/*.cpp)
OBJFILES   := $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SRCFILES))

# Test source files: every .cpp in TESTSRCDIR
TEST_SRCS  := $(wildcard $(TESTSRCDIR)/*.cpp)
# Test executables: output them to TESTBINDIR, with the same base name as the source file.
TEST_EXES  := $(patsubst $(TESTSRCDIR)/%.cpp, $(TESTBINDIR)/%, $(TEST_SRCS))

# Phony targets
.PHONY: all lib tests clean

# Default: build everything.
all: lib tests

# --------------------------
# Library build rules
# --------------------------

lib: $(LIB)

# Archive object files into the static library.
$(LIB): $(OBJFILES)
	@mkdir -p $(LIBDIR)
	$(AR) $(ARFLAGS) $@ $^

# Pattern rule to compile each source file into an object file.
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# --------------------------
# Test build rules
# --------------------------

tests: $(TEST_EXES)

# Pattern rule to compile each test source file into an executable placed in TESTBINDIR.
$(TESTBINDIR)/%: $(TESTSRCDIR)/%.cpp $(LIB)
	@mkdir -p $(TESTBINDIR)
	$(CXX) $(CXXFLAGS) $< -L$(LIBDIR) -l$(LIBNAME) -o $@

# --------------------------
# Clean up build artifacts
# --------------------------
clean:
	rm -rf $(BUILDDIR)
