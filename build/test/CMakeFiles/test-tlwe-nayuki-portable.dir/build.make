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
include test/CMakeFiles/test-tlwe-nayuki-portable.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/test-tlwe-nayuki-portable.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/test-tlwe-nayuki-portable.dir/flags.make

test/CMakeFiles/test-tlwe-nayuki-portable.dir/test-tlwe.cpp.o: test/CMakeFiles/test-tlwe-nayuki-portable.dir/flags.make
test/CMakeFiles/test-tlwe-nayuki-portable.dir/test-tlwe.cpp.o: /home/ding/桌面/tfhe-master原版\ (三次修改FFT)-3密钥合并/src/test/test-tlwe.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/ding/桌面/tfhe-master原版 (三次修改FFT)-3密钥合并/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/test-tlwe-nayuki-portable.dir/test-tlwe.cpp.o"
	cd "/home/ding/桌面/tfhe-master原版 (三次修改FFT)-3密钥合并/build/test" && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/test-tlwe-nayuki-portable.dir/test-tlwe.cpp.o -c "/home/ding/桌面/tfhe-master原版 (三次修改FFT)-3密钥合并/src/test/test-tlwe.cpp"

test/CMakeFiles/test-tlwe-nayuki-portable.dir/test-tlwe.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test-tlwe-nayuki-portable.dir/test-tlwe.cpp.i"
	cd "/home/ding/桌面/tfhe-master原版 (三次修改FFT)-3密钥合并/build/test" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/ding/桌面/tfhe-master原版 (三次修改FFT)-3密钥合并/src/test/test-tlwe.cpp" > CMakeFiles/test-tlwe-nayuki-portable.dir/test-tlwe.cpp.i

test/CMakeFiles/test-tlwe-nayuki-portable.dir/test-tlwe.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test-tlwe-nayuki-portable.dir/test-tlwe.cpp.s"
	cd "/home/ding/桌面/tfhe-master原版 (三次修改FFT)-3密钥合并/build/test" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/ding/桌面/tfhe-master原版 (三次修改FFT)-3密钥合并/src/test/test-tlwe.cpp" -o CMakeFiles/test-tlwe-nayuki-portable.dir/test-tlwe.cpp.s

test/CMakeFiles/test-tlwe-nayuki-portable.dir/test-tlwe.cpp.o.requires:

.PHONY : test/CMakeFiles/test-tlwe-nayuki-portable.dir/test-tlwe.cpp.o.requires

test/CMakeFiles/test-tlwe-nayuki-portable.dir/test-tlwe.cpp.o.provides: test/CMakeFiles/test-tlwe-nayuki-portable.dir/test-tlwe.cpp.o.requires
	$(MAKE) -f test/CMakeFiles/test-tlwe-nayuki-portable.dir/build.make test/CMakeFiles/test-tlwe-nayuki-portable.dir/test-tlwe.cpp.o.provides.build
.PHONY : test/CMakeFiles/test-tlwe-nayuki-portable.dir/test-tlwe.cpp.o.provides

test/CMakeFiles/test-tlwe-nayuki-portable.dir/test-tlwe.cpp.o.provides.build: test/CMakeFiles/test-tlwe-nayuki-portable.dir/test-tlwe.cpp.o


# Object files for target test-tlwe-nayuki-portable
test__tlwe__nayuki__portable_OBJECTS = \
"CMakeFiles/test-tlwe-nayuki-portable.dir/test-tlwe.cpp.o"

# External object files for target test-tlwe-nayuki-portable
test__tlwe__nayuki__portable_EXTERNAL_OBJECTS =

test/test-tlwe-nayuki-portable: test/CMakeFiles/test-tlwe-nayuki-portable.dir/test-tlwe.cpp.o
test/test-tlwe-nayuki-portable: test/CMakeFiles/test-tlwe-nayuki-portable.dir/build.make
test/test-tlwe-nayuki-portable: libtfhe/libtfhe-nayuki-portable.a
test/test-tlwe-nayuki-portable: test/CMakeFiles/test-tlwe-nayuki-portable.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/home/ding/桌面/tfhe-master原版 (三次修改FFT)-3密钥合并/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable test-tlwe-nayuki-portable"
	cd "/home/ding/桌面/tfhe-master原版 (三次修改FFT)-3密钥合并/build/test" && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test-tlwe-nayuki-portable.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/test-tlwe-nayuki-portable.dir/build: test/test-tlwe-nayuki-portable

.PHONY : test/CMakeFiles/test-tlwe-nayuki-portable.dir/build

test/CMakeFiles/test-tlwe-nayuki-portable.dir/requires: test/CMakeFiles/test-tlwe-nayuki-portable.dir/test-tlwe.cpp.o.requires

.PHONY : test/CMakeFiles/test-tlwe-nayuki-portable.dir/requires

test/CMakeFiles/test-tlwe-nayuki-portable.dir/clean:
	cd "/home/ding/桌面/tfhe-master原版 (三次修改FFT)-3密钥合并/build/test" && $(CMAKE_COMMAND) -P CMakeFiles/test-tlwe-nayuki-portable.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/test-tlwe-nayuki-portable.dir/clean

test/CMakeFiles/test-tlwe-nayuki-portable.dir/depend:
	cd "/home/ding/桌面/tfhe-master原版 (三次修改FFT)-3密钥合并/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/home/ding/桌面/tfhe-master原版 (三次修改FFT)-3密钥合并/src" "/home/ding/桌面/tfhe-master原版 (三次修改FFT)-3密钥合并/src/test" "/home/ding/桌面/tfhe-master原版 (三次修改FFT)-3密钥合并/build" "/home/ding/桌面/tfhe-master原版 (三次修改FFT)-3密钥合并/build/test" "/home/ding/桌面/tfhe-master原版 (三次修改FFT)-3密钥合并/build/test/CMakeFiles/test-tlwe-nayuki-portable.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : test/CMakeFiles/test-tlwe-nayuki-portable.dir/depend

