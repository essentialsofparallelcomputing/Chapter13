cmake_minimum_required(VERSION 3.10)
project (ShallowWater)

set (CMAKE_C_STANDARD 99)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
#if (NOT TARGET distclean)
#   INCLUDE( cmake/Modules/TargetDistclean.cmake OPTIONAL)
#endif (NOT TARGET distclean)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O3")

if (NOT CMAKE_OPENACC_VERBOSE)
    set(CMAKE_OPENACC_VERBOSE true)
endif (NOT CMAKE_OPENACC_VERBOSE)

find_package(Vector)
if (CMAKE_VECTOR_VERBOSE)
    set(VECTOR_C_FLAGS "${VECTOR_C_FLAGS} ${VECTOR_C_VERBOSE}")
endif (CMAKE_VECTOR_VERBOSE)
if (CMAKE_VECTOR_NOVEC)
    set(VECTOR_C_FLAGS "${VECTOR_NOVEC_C_FLAGS}")
endif (CMAKE_VECTOR_NOVEC)

find_package(OpenACC)
if ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU" AND "${CMAKE_C_COMPILER_VERSION}" VERSION_LESS_EQUAL "7.4.0")
   set(OpenACC_C_FOUND FALSE)
endif ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU" AND "${CMAKE_C_COMPILER_VERSION}" VERSION_LESS_EQUAL "7.4.0")

if (CMAKE_OPENACC_VERBOSE)
    set(OpenACC_C_FLAGS "${OpenACC_C_FLAGS} ${OpenACC_C_VERBOSE}")
endif (CMAKE_OPENACC_VERBOSE)

add_subdirectory(graphics)

if (ENABLE_GRAPHICS)
   add_definitions(-DHAVE_GRAPHICS=1)
   include_directories(graphics)
endif (ENABLE_GRAPHICS)

# Adds build target of ShallowWater with source code files
add_executable(ShallowWater ShallowWater.c memory.c timer.c)
target_link_libraries(ShallowWater -lm)
if (ENABLE_GRAPHICS)
   target_link_libraries(ShallowWater graphics ${FREEGLUT_LIBRARIES})
endif (ENABLE_GRAPHICS)

if (ImageMagick_MagickWand_FOUND)
   add_definitions(-DHAVE_MAGICKWAND)
endif (ImageMagick_MagickWand_FOUND)

# Adds build target of ShallowWater_par1 with source code files
if (OpenACC_C_FOUND)
   add_executable(ShallowWater_par1 ShallowWater_par1.c memory.c timer.c)
   set_source_files_properties(ShallowWater_par1.c PROPERTIES COMPILE_FLAGS "${OpenACC_C_FLAGS} ${VECTOR_C_FLAGS}")
   set_target_properties(ShallowWater_par1 PROPERTIES LINK_FLAGS "${OpenACC_C_FLAGS}")
   target_link_libraries(ShallowWater_par1 -lm)
endif (OpenACC_C_FOUND)

# Adds build target of ShallowWater_par2 with source code files
if (OpenACC_C_FOUND)
   add_executable(ShallowWater_par2 ShallowWater_par2.c memory.c timer.c)
   set_source_files_properties(ShallowWater_par2.c PROPERTIES COMPILE_FLAGS "${OpenACC_C_FLAGS} ${VECTOR_C_FLAGS}")
   set_target_properties(ShallowWater_par2 PROPERTIES LINK_FLAGS "${OpenACC_C_FLAGS}")
   target_link_libraries(ShallowWater_par2 -lm)
endif (OpenACC_C_FOUND)

# Adds build target of ShallowWater_par3 with source code files
if (OpenACC_C_FOUND)
   add_executable(ShallowWater_par3 ShallowWater_par3.c memory.c timer.c)
   set_source_files_properties(ShallowWater_par3.c PROPERTIES COMPILE_FLAGS "${OpenACC_C_FLAGS} ${VECTOR_C_FLAGS}")
   set_target_properties(ShallowWater_par3 PROPERTIES LINK_FLAGS "${OpenACC_C_FLAGS}")
   target_link_libraries(ShallowWater_par3 -lm)
endif (OpenACC_C_FOUND)

# Adds build target of ShallowWater_par4 with source code files
if (OpenACC_C_FOUND)
   add_executable(ShallowWater_par4 ShallowWater_par4.c memory.c timer.c)
   set_source_files_properties(ShallowWater_par4.c PROPERTIES COMPILE_FLAGS "${OpenACC_C_FLAGS} ${VECTOR_C_FLAGS}")
   set_target_properties(ShallowWater_par4 PROPERTIES LINK_FLAGS "${OpenACC_C_FLAGS}")
   target_link_libraries(ShallowWater_par4 -lm)
endif (OpenACC_C_FOUND)

# Cleanup
add_custom_target(distclean COMMAND rm -rf CMakeCache.txt CMakeFiles
                  Makefile cmake_install.cmake ShallowWater.dSYM ipo_out.optrpt)
