# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.8

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
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/home/ding/桌面/tfhe-master原版 (三次修改FFT)-3密钥合并/src"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/home/ding/桌面/tfhe-master原版 (三次修改FFT)-3密钥合并/build"

# Include any dependencies generated for this target.
include test/CMakeFiles/test-lwe-fftw.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/test-lwe-fftw.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/test-lwe-fftw.dir/flags.make

test/CMakeFiles/test-lwe-fftw.dir/test-lwe.cpp.o: test/CMakeFiles/test-lwe-fftw.dir/flags.make
test/CMakeFiles/test-lwe-fftw.dir/test-lwe.cpp.o: /home/ding/桌面/tfhe-master原版\ (三次修改FFT)-3密钥合并/src/test/test-lwe.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/ding/桌面/tfhe-master原版 (三次修改FFT)-3密钥合并/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/test-lwe-fftw.dir/test-lwe.cpp.o"
	cd "/home/ding/桌面/tfhe-master原版 (三次修改FFT)-3密钥合并/build/test" && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/test-lwe-fftw.dir/test-lwe.cpp.o -c "/home/ding/桌面/tfhe-master原版 (三次修改FFT)-3密钥合并/src/test/test-lwe.cpp"

test/CMakeFiles/test-lwe-fftw.dir/test-lwe.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test-lwe-fftw.dir/test-lwe.cpp.i"
	cd "/home/ding/桌面/tfhe-master原版 (三次修改FFT)-3密钥合并/build/test" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/ding/桌面/tfhe-master原版 (三次修改FFT)-3密钥合并/src/test/test-lwe.cpp" > CMakeFiles/test-lwe-fftw.dir/test-lwe.cpp.i

test/CMakeFiles/test-lwe-fftw.dir/test-lwe.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test-lwe-fftw.dir/test-lwe.cpp.s"
	cd "/home/ding/桌面/tfhe-master原版 (三次修改FFT)-3密钥合并/build/test" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/ding/桌面/tfhe-master原版 (三次修改FFT)-3密钥合并/src/test/test-lwe.cpp" -o CMakeFiles/test-lwe-fftw.dir/test-lwe.cpp.s

test/CMakeFiles/test-lwe-fftw.dir/test-lwe.cpp.o.requires:

.PHONY : test/CMakeFiles/test-lwe-fftw.dir/test-lwe.cpp.o.requires

test/CMakeFiles/test-lwe-fftw.dir/test-lwe.cpp.o.provides: test/CMakeFiles/test-lwe-fftw.dir/test-lwe.cpp.o.requires
	$(MAKE) -f test/CMakeFiles/test-lwe-fftw.dir/build.make test/CMakeFiles/test-lwe-fftw.dir/test-lwe.cpp.o.provides.build
.PHONY : test/CMakeFiles/test-lwe-fftw.dir/test-lwe.cpp.o.provides

test/CMakeFiles/test-lwe-fftw.dir/test-lwe.cpp.o.provides.build: test/CMakeFiles/test-lwe-fftw.dir/test-lwe.cpp.o


# Object files for target test-lwe-fftw
test__lwe__fftw_OBJECTS = \
"CMakeFiles/test-lwe-fftw.dir/test-lwe.cpp.o"

# External object files for target test-lwe-fftw
test__lwe__fftw_EXTERNAL_OBJECTS =

test/test-lwe-fftw: test/CMakeFiles/test-lwe-fftw.dir/test-lwe.cpp.o
test/test-lwe-fftw: test/CMakeFiles/test-lwe-fftw.dir/build.make
test/test-lwe-fftw: libtfhe/libtfhe-fftw.a
test/test-lwe-fftw: /usr/local/lib/libfftw3.a
test/test-lwe-fftw: test/CMakeFiles/test-lwe-fftw.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/home/ding/桌面/tfhe-master原版 (三次修改FFT)-3密钥合并/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable test-lwe-fftw"
	cd "/home/ding/桌面/tfhe-master原版 (三次修改FFT)-3密钥合并/build/test" && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test-lwe-fftw.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/test-lwe-fftw.dir/build: test/test-lwe-fftw

.PHONY : test/CMakeFiles/test-lwe-fftw.dir/build

test/CMakeFiles/test-lwe-fftw.dir/requires: test/CMakeFiles/test-lwe-fftw.dir/test-lwe.cpp.o.requires

.PHONY : test/CMakeFiles/test-lwe-fftw.dir/requires

test/CMakeFiles/test-lwe-fftw.dir/clean:
	cd "/home/ding/桌面/tfhe-master原版 (三次修改FFT)-3密钥合并/build/test" && $(CMAKE_COMMAND) -P CMakeFiles/test-lwe-fftw.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/test-lwe-fftw.dir/clean

test/CMakeFiles/test-lwe-fftw.dir/depend:
	cd "/home/ding/桌面/tfhe-master原版 (三次修改FFT)-3密钥合并/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/home/ding/桌面/tfhe-master原版 (三次修改FFT)-3密钥合并/src" "/home/ding/桌面/tfhe-master原版 (三次修改FFT)-3密钥合并/src/test" "/home/ding/桌面/tfhe-master原版 (三次修改FFT)-3密钥合并/build" "/home/ding/桌面/tfhe-master原版 (三次修改FFT)-3密钥合并/build/test" "/home/ding/桌面/tfhe-master原版 (三次修改FFT)-3密钥合并/build/test/CMakeFiles/test-lwe-fftw.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : test/CMakeFiles/test-lwe-fftw.dir/depend

