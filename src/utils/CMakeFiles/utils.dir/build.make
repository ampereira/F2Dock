# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
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
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E remove -f

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/local/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/andre/F2Dock-refactored

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/andre/F2Dock-refactored

# Include any dependencies generated for this target.
include src/utils/CMakeFiles/utils.dir/depend.make

# Include the progress variables for this target.
include src/utils/CMakeFiles/utils.dir/progress.make

# Include the compile flags for this target's objects.
include src/utils/CMakeFiles/utils.dir/flags.make

src/utils/CMakeFiles/utils.dir/utils.cpp.o: src/utils/CMakeFiles/utils.dir/flags.make
src/utils/CMakeFiles/utils.dir/utils.cpp.o: src/utils/utils.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/andre/F2Dock-refactored/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/utils/CMakeFiles/utils.dir/utils.cpp.o"
	cd /Users/andre/F2Dock-refactored/src/utils && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/utils.dir/utils.cpp.o -c /Users/andre/F2Dock-refactored/src/utils/utils.cpp

src/utils/CMakeFiles/utils.dir/utils.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/utils.dir/utils.cpp.i"
	cd /Users/andre/F2Dock-refactored/src/utils && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /Users/andre/F2Dock-refactored/src/utils/utils.cpp > CMakeFiles/utils.dir/utils.cpp.i

src/utils/CMakeFiles/utils.dir/utils.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/utils.dir/utils.cpp.s"
	cd /Users/andre/F2Dock-refactored/src/utils && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /Users/andre/F2Dock-refactored/src/utils/utils.cpp -o CMakeFiles/utils.dir/utils.cpp.s

src/utils/CMakeFiles/utils.dir/utils.cpp.o.requires:
.PHONY : src/utils/CMakeFiles/utils.dir/utils.cpp.o.requires

src/utils/CMakeFiles/utils.dir/utils.cpp.o.provides: src/utils/CMakeFiles/utils.dir/utils.cpp.o.requires
	$(MAKE) -f src/utils/CMakeFiles/utils.dir/build.make src/utils/CMakeFiles/utils.dir/utils.cpp.o.provides.build
.PHONY : src/utils/CMakeFiles/utils.dir/utils.cpp.o.provides

src/utils/CMakeFiles/utils.dir/utils.cpp.o.provides.build: src/utils/CMakeFiles/utils.dir/utils.cpp.o

# Object files for target utils
utils_OBJECTS = \
"CMakeFiles/utils.dir/utils.cpp.o"

# External object files for target utils
utils_EXTERNAL_OBJECTS =

lib/libutils.a: src/utils/CMakeFiles/utils.dir/utils.cpp.o
lib/libutils.a: src/utils/CMakeFiles/utils.dir/build.make
lib/libutils.a: src/utils/CMakeFiles/utils.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX static library ../../lib/libutils.a"
	cd /Users/andre/F2Dock-refactored/src/utils && $(CMAKE_COMMAND) -P CMakeFiles/utils.dir/cmake_clean_target.cmake
	cd /Users/andre/F2Dock-refactored/src/utils && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/utils.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/utils/CMakeFiles/utils.dir/build: lib/libutils.a
.PHONY : src/utils/CMakeFiles/utils.dir/build

src/utils/CMakeFiles/utils.dir/requires: src/utils/CMakeFiles/utils.dir/utils.cpp.o.requires
.PHONY : src/utils/CMakeFiles/utils.dir/requires

src/utils/CMakeFiles/utils.dir/clean:
	cd /Users/andre/F2Dock-refactored/src/utils && $(CMAKE_COMMAND) -P CMakeFiles/utils.dir/cmake_clean.cmake
.PHONY : src/utils/CMakeFiles/utils.dir/clean

src/utils/CMakeFiles/utils.dir/depend:
	cd /Users/andre/F2Dock-refactored && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/andre/F2Dock-refactored /Users/andre/F2Dock-refactored/src/utils /Users/andre/F2Dock-refactored /Users/andre/F2Dock-refactored/src/utils /Users/andre/F2Dock-refactored/src/utils/CMakeFiles/utils.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/utils/CMakeFiles/utils.dir/depend

