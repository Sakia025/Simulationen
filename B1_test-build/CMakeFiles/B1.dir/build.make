# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.16.4/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.16.4/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/smller/Simulationen/B1_test

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/smller/Simulationen/B1_test-build

# Utility rule file for B1.

# Include the progress variables for this target.
include CMakeFiles/B1.dir/progress.make

CMakeFiles/B1: exampleB1


B1: CMakeFiles/B1
B1: CMakeFiles/B1.dir/build.make

.PHONY : B1

# Rule to build all files generated by this target.
CMakeFiles/B1.dir/build: B1

.PHONY : CMakeFiles/B1.dir/build

CMakeFiles/B1.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/B1.dir/cmake_clean.cmake
.PHONY : CMakeFiles/B1.dir/clean

CMakeFiles/B1.dir/depend:
	cd /Users/smller/Simulationen/B1_test-build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/smller/Simulationen/B1_test /Users/smller/Simulationen/B1_test /Users/smller/Simulationen/B1_test-build /Users/smller/Simulationen/B1_test-build /Users/smller/Simulationen/B1_test-build/CMakeFiles/B1.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/B1.dir/depend

