# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canoncical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /org/centers/cvc/software/share/usr.linux.x86_64/bin/cmake

# The command to remove a file.
RM = /org/centers/cvc/software/share/usr.linux.x86_64/bin/cmake -E remove -f

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /org/centers/cvc/software/share/usr.linux.x86_64/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /h1/apereira/F2Dock

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /h1/apereira/F2Dock

# Include any dependencies generated for this target.
include src/fast-LJ/CMakeFiles/fast-LJ.dir/depend.make

# Include the progress variables for this target.
include src/fast-LJ/CMakeFiles/fast-LJ.dir/progress.make

# Include the compile flags for this target's objects.
include src/fast-LJ/CMakeFiles/fast-LJ.dir/flags.make

src/fast-LJ/CMakeFiles/fast-LJ.dir/LJTest.cpp.o: src/fast-LJ/CMakeFiles/fast-LJ.dir/flags.make
src/fast-LJ/CMakeFiles/fast-LJ.dir/LJTest.cpp.o: src/fast-LJ/LJTest.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /h1/apereira/F2Dock/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/fast-LJ/CMakeFiles/fast-LJ.dir/LJTest.cpp.o"
	cd /h1/apereira/F2Dock/src/fast-LJ && /h1/apereira/OpenMPI/bin/mpic++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/fast-LJ.dir/LJTest.cpp.o -c /h1/apereira/F2Dock/src/fast-LJ/LJTest.cpp

src/fast-LJ/CMakeFiles/fast-LJ.dir/LJTest.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fast-LJ.dir/LJTest.cpp.i"
	cd /h1/apereira/F2Dock/src/fast-LJ && /h1/apereira/OpenMPI/bin/mpic++  $(CXX_DEFINES) $(CXX_FLAGS) -E /h1/apereira/F2Dock/src/fast-LJ/LJTest.cpp > CMakeFiles/fast-LJ.dir/LJTest.cpp.i

src/fast-LJ/CMakeFiles/fast-LJ.dir/LJTest.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fast-LJ.dir/LJTest.cpp.s"
	cd /h1/apereira/F2Dock/src/fast-LJ && /h1/apereira/OpenMPI/bin/mpic++  $(CXX_DEFINES) $(CXX_FLAGS) -S /h1/apereira/F2Dock/src/fast-LJ/LJTest.cpp -o CMakeFiles/fast-LJ.dir/LJTest.cpp.s

src/fast-LJ/CMakeFiles/fast-LJ.dir/LJTest.cpp.o.requires:
.PHONY : src/fast-LJ/CMakeFiles/fast-LJ.dir/LJTest.cpp.o.requires

src/fast-LJ/CMakeFiles/fast-LJ.dir/LJTest.cpp.o.provides: src/fast-LJ/CMakeFiles/fast-LJ.dir/LJTest.cpp.o.requires
	$(MAKE) -f src/fast-LJ/CMakeFiles/fast-LJ.dir/build.make src/fast-LJ/CMakeFiles/fast-LJ.dir/LJTest.cpp.o.provides.build
.PHONY : src/fast-LJ/CMakeFiles/fast-LJ.dir/LJTest.cpp.o.provides

src/fast-LJ/CMakeFiles/fast-LJ.dir/LJTest.cpp.o.provides.build: src/fast-LJ/CMakeFiles/fast-LJ.dir/LJTest.cpp.o
.PHONY : src/fast-LJ/CMakeFiles/fast-LJ.dir/LJTest.cpp.o.provides.build

src/fast-LJ/CMakeFiles/fast-LJ.dir/fastLJ.cpp.o: src/fast-LJ/CMakeFiles/fast-LJ.dir/flags.make
src/fast-LJ/CMakeFiles/fast-LJ.dir/fastLJ.cpp.o: src/fast-LJ/fastLJ.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /h1/apereira/F2Dock/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/fast-LJ/CMakeFiles/fast-LJ.dir/fastLJ.cpp.o"
	cd /h1/apereira/F2Dock/src/fast-LJ && /h1/apereira/OpenMPI/bin/mpic++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/fast-LJ.dir/fastLJ.cpp.o -c /h1/apereira/F2Dock/src/fast-LJ/fastLJ.cpp

src/fast-LJ/CMakeFiles/fast-LJ.dir/fastLJ.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fast-LJ.dir/fastLJ.cpp.i"
	cd /h1/apereira/F2Dock/src/fast-LJ && /h1/apereira/OpenMPI/bin/mpic++  $(CXX_DEFINES) $(CXX_FLAGS) -E /h1/apereira/F2Dock/src/fast-LJ/fastLJ.cpp > CMakeFiles/fast-LJ.dir/fastLJ.cpp.i

src/fast-LJ/CMakeFiles/fast-LJ.dir/fastLJ.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fast-LJ.dir/fastLJ.cpp.s"
	cd /h1/apereira/F2Dock/src/fast-LJ && /h1/apereira/OpenMPI/bin/mpic++  $(CXX_DEFINES) $(CXX_FLAGS) -S /h1/apereira/F2Dock/src/fast-LJ/fastLJ.cpp -o CMakeFiles/fast-LJ.dir/fastLJ.cpp.s

src/fast-LJ/CMakeFiles/fast-LJ.dir/fastLJ.cpp.o.requires:
.PHONY : src/fast-LJ/CMakeFiles/fast-LJ.dir/fastLJ.cpp.o.requires

src/fast-LJ/CMakeFiles/fast-LJ.dir/fastLJ.cpp.o.provides: src/fast-LJ/CMakeFiles/fast-LJ.dir/fastLJ.cpp.o.requires
	$(MAKE) -f src/fast-LJ/CMakeFiles/fast-LJ.dir/build.make src/fast-LJ/CMakeFiles/fast-LJ.dir/fastLJ.cpp.o.provides.build
.PHONY : src/fast-LJ/CMakeFiles/fast-LJ.dir/fastLJ.cpp.o.provides

src/fast-LJ/CMakeFiles/fast-LJ.dir/fastLJ.cpp.o.provides.build: src/fast-LJ/CMakeFiles/fast-LJ.dir/fastLJ.cpp.o
.PHONY : src/fast-LJ/CMakeFiles/fast-LJ.dir/fastLJ.cpp.o.provides.build

# Object files for target fast-LJ
fast__LJ_OBJECTS = \
"CMakeFiles/fast-LJ.dir/LJTest.cpp.o" \
"CMakeFiles/fast-LJ.dir/fastLJ.cpp.o"

# External object files for target fast-LJ
fast__LJ_EXTERNAL_OBJECTS =

lib/libfast-LJ.a: src/fast-LJ/CMakeFiles/fast-LJ.dir/LJTest.cpp.o
lib/libfast-LJ.a: src/fast-LJ/CMakeFiles/fast-LJ.dir/fastLJ.cpp.o
lib/libfast-LJ.a: src/fast-LJ/CMakeFiles/fast-LJ.dir/build.make
lib/libfast-LJ.a: src/fast-LJ/CMakeFiles/fast-LJ.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX static library ../../lib/libfast-LJ.a"
	cd /h1/apereira/F2Dock/src/fast-LJ && $(CMAKE_COMMAND) -P CMakeFiles/fast-LJ.dir/cmake_clean_target.cmake
	cd /h1/apereira/F2Dock/src/fast-LJ && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/fast-LJ.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/fast-LJ/CMakeFiles/fast-LJ.dir/build: lib/libfast-LJ.a
.PHONY : src/fast-LJ/CMakeFiles/fast-LJ.dir/build

src/fast-LJ/CMakeFiles/fast-LJ.dir/requires: src/fast-LJ/CMakeFiles/fast-LJ.dir/LJTest.cpp.o.requires
src/fast-LJ/CMakeFiles/fast-LJ.dir/requires: src/fast-LJ/CMakeFiles/fast-LJ.dir/fastLJ.cpp.o.requires
.PHONY : src/fast-LJ/CMakeFiles/fast-LJ.dir/requires

src/fast-LJ/CMakeFiles/fast-LJ.dir/clean:
	cd /h1/apereira/F2Dock/src/fast-LJ && $(CMAKE_COMMAND) -P CMakeFiles/fast-LJ.dir/cmake_clean.cmake
.PHONY : src/fast-LJ/CMakeFiles/fast-LJ.dir/clean

src/fast-LJ/CMakeFiles/fast-LJ.dir/depend:
	cd /h1/apereira/F2Dock && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /h1/apereira/F2Dock /h1/apereira/F2Dock/src/fast-LJ /h1/apereira/F2Dock /h1/apereira/F2Dock/src/fast-LJ /h1/apereira/F2Dock/src/fast-LJ/CMakeFiles/fast-LJ.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/fast-LJ/CMakeFiles/fast-LJ.dir/depend

