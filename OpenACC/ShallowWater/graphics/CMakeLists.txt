cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
# --- custom targets: ---
#if (NOT TARGET distclean)
#   INCLUDE( cmake/Modules/TargetDistclean.cmake OPTIONAL)
#endif (NOT TARGET distclean)

# enable @rpath in the install name for any shared library being built
set(CMAKE_MACOSX_RPATH 1)

project (GRAPHICS)
# The version number.
set (GRAPHICS_VERSION_MAJOR 2)
set (GRAPHICS_VERSION_MINOR 0.7)
set (GRAPHICS_VERSION "${GRAPHICS_VERSION_MAJOR}.${GRAPHICS_VERSION_MINOR}")

if (NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif(NOT CMAKE_BUILD_TYPE)

# uncomment to build debug
#set(CMAKE_BUILD_TYPE Debug)

if ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
   # using Clang
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ftree-vectorize")
   #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ftree-vectorize -ftree-vectorizer-verbose=5")
   set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -g -O0 -Wextra")
   set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g -O3")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
   # using GCC
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ftree-vectorize")
   #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ftree-vectorize -ftree-vectorizer-verbose=5")
   set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -g -O0 -Wextra")
   set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g -O3")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
   # using Intel C
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "MSVC")
   # using Visual Studio C++
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "XL")
   # using IBM XL Compiler
   set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g -O0")
   set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g -O3")
endif()

if (APPLE AND "${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native -mtune=native -Wa,-q")
else (APPLE AND "${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native -mtune=native")
endif (APPLE AND "${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")


#set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -pedantic")

# FFLAGS depend on the compiler
#if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
#  # gfortran
#  set (CMAKE_Fortran_FLAGS "-cpp -fdefault-real-8 -fdefault-double-8 -ffixed-line-length-132")
#  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -g")
#  set (CMAKE_Fortran_FLAGS_DEBUG   "-g -fbounds-check -Wall")
#  set (CMAKE_Fortran_FLAGS_GCOV   "-g -fprofile-arcs -ftest-coverage -fbounds-check -Wall")
#elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
#  #set (CMAKE_Fortran_FLAGS "-std90 -cpp -r8 -132 -fp-model source -fp-speculation strict -vec-report3 -align all -fPIC")
#  #set (CMAKE_Fortran_FLAGS_RELEASE "-g -O3 -xHost")
#  set (CMAKE_Fortran_FLAGS "-cpp -r8 -extend_source")
#  set (CMAKE_Fortran_FLAGS_RELEASE "-g -O3 -xHost")
#  set (CMAKE_Fortran_FLAGS_DEBUG   "-g -O0 -traceback -fp-stack-check -ftrapuv -check bounds -check uninit -check noarg_temp_created -fpe0")
#elseif (CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
#  set (CMAKE_Fortran_FLAGS "-r8 -Mpreprocess -Mextend ")
#  set (CMAKE_Fortran_FLAGS_RELEASE "-g -O3")
#  set (CMAKE_Fortran_FLAGS_DEBUG   "-g -Mbounds -traceback")
#elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
#  #set (CMAKE_Fortran_FLAGS   "-s default32")
#  set (CMAKE_Fortran_FLAGS_DEBUG   "-g -O0")
#  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
#elseif (CMAKE_Fortran_COMPILER_ID MATCHES "XL")
#  #set (CMAKE_Fortran_FLAGS   "-s default32")
#  set (CMAKE_Fortran_FLAGS_DEBUG   "-g -O0")
#  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
#else ()
#  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
#  message ("Fortran compiler name: " ${Fortran_COMPILER_NAME})
#  message ("Fortran compiler id: " ${CMAKE_Fortran_COMPILER_ID})
#  message ("TopLevel: No optimized Fortran compiler flags are known, we just try -O2...")
#  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
#  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
#endif ()

if (NOT DEFINED OPENCL_FOUND)
   find_package(OpenCL)
   #  OPENCL_FOUND        : TRUE if found, FALSE otherwise
   #  OPENCL_INCLUDE_DIRS : Include directories for OpenCL
   #  OPENCL_LIBRARIES    : The libraries to link against
   if (OPENCL_FOUND)
      include_directories(${OPENCL_INCLUDE_DIRS})
   endif (OPENCL_FOUND)
endif (NOT DEFINED OPENCL_FOUND)

#find_package(ImageMagick COMPONENTS MagickCore MagickWand)
find_package(ImageMagick COMPONENTS MagickWand)
#   ImageMagick_FOUND                    - TRUE if all components are found.
#   ImageMagick_INCLUDE_DIRS             - Full paths to all include dirs.
#   ImageMagick_LIBRARIES                - Full paths to all libraries.
#   ImageMagick_<component>_FOUND        - TRUE if <component> is found.
#   ImageMagick_<component>_INCLUDE_DIRS - Full path to <component> include dirs.
#   ImageMagick_<component>_LIBRARIES    - Full path to <component> libraries.
if (ImageMagick_FOUND)
   if (ImageMagick_MagickWand_FOUND)
      add_definitions(-DHAVE_MAGICKWAND)
      add_definitions(-DMAGICKCORE_HDRI_ENABLE=0)
      add_definitions(-DMAGICKCORE_QUANTUM_DEPTH=16)
#     add_definitions(-DMAGICKCORE_HDRI_ENABLE=0)
#     add_definitions(-DMAGICKCORE_QUANTUM_DEPTH=16)
      include_directories(${ImageMagick_INCLUDE_DIRS})
   endif (ImageMagick_MagickWand_FOUND)

#   message("ImageMagick_FOUND                   ${ImageMagick_FOUND}")
#   message("ImageMagick_INCLUDE_DIRS            ${ImageMagick_INCLUDE_DIRS}")
#   message("ImageMagick_LIBRARIES               ${ImageMagick_LIBRARIES}")
##  message("ImageMagick_MagickCore_FOUND        ${ImageMagick_MagickCore_FOUND}")
##  message("ImageMagick_MagickCore_INCLUDE_DIRS ${ImageMagick_MagickCore_INCLUDE_DIRS}")
##  message("ImageMagick_MagickCore_LIBRARIES    ${ImageMagick_MagickCore_LIBRARIES}")
#   message("ImageMagick_MagickWand_FOUND        ${ImageMagick_MagickWand_FOUND}")
#   message("ImageMagick_MagickWand_INCLUDE_DIRS ${ImageMagick_MagickWand_INCLUDE_DIRS}")
#   message("ImageMagick_MagickWand_LIBRARIES    ${ImageMagick_MagickWand_LIBRARIES}")
endif (ImageMagick_FOUND)

if (NOT DEFINED MPI_FOUND)
#  find_package(MPI)
endif (NOT DEFINED MPI_FOUND)
if (MPI_FOUND)
   include_directories(${MPI_INCLUDE_PATH})
endif (MPI_FOUND)

if (NOT DEFINED OPENMP_FOUND)
   if (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} VERSION_GREATER 3.8.9)
     cmake_policy(PUSH)
     # Suppress warnings in FindOpenMP
     cmake_policy(SET CMP0054 OLD)
   endif (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} VERSION_GREATER 3.8.9)
   find_package( OpenMP)
   if (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} VERSION_GREATER 3.8.9)
     cmake_policy(POP)
   endif (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} VERSION_GREATER 3.8.9)
endif (NOT DEFINED OPENMP_FOUND)

if (NOT PRECISION_TYPE)
   set (PRECISION_TYPE full_precision)
endif (NOT PRECISION_TYPE)

if (PRECISION_TYPE STREQUAL "mixed_precision")
   set (MIXED_PRECISION on)
   set (PRECISION_TYPE mixed_precision)
   add_definitions(-DMIXED_PRECISION)
endif (PRECISION_TYPE STREQUAL "mixed_precision")

if (PRECISION_TYPE STREQUAL "full_precision")
   set (FULL_PRECISION on)
   set (PRECISION_TYPE full_precision)
   add_definitions(-DFULL_PRECISION)
endif (PRECISION_TYPE STREQUAL "full_precision")

if (PRECISION_TYPE STREQUAL "minimum_precision")
   set (MINIMUM_PRECISION on)
   set (PRECISION_TYPE minimum_precision)
   add_definitions(-DMINIMUM_PRECISION)
endif (PRECISION_TYPE STREQUAL "minimum_precision")

if (NOT GRAPHICS_TYPE)
   set (GRAPHICS_TYPE none)
endif (NOT GRAPHICS_TYPE)

if (GRAPHICS_TYPE STREQUAL "OpenGL")
   find_package(OpenGL)
endif (GRAPHICS_TYPE STREQUAL "OpenGL")

# If autodetect does not work try setting MPEHOME or
# MPE_INCLUDE_DIR and MPE_LIBRARIES
if (GRAPHICS_TYPE STREQUAL "MPE")
   find_package(MPE)
   find_package(X11)
endif (GRAPHICS_TYPE STREQUAL "MPE")

#message("MPI_LIBRARIES ${MPI_LIBRARIES}")
#message("MPI_INCLUDE_PATH ${MPI_INCLUDE_PATH}")
#message("GRAPHICS_TYPE is ${GRAPHICS_TYPE}")
#message("MPE_INCLUDE is ${MPE_INCLUDE_DIR}")

if (OPENGL_FOUND)
   set (HAVE_GRAPHICS on)
   set (HAVE_OPENGL on)
   include_directories(${OPENGL_INCLUDE_DIR})
   find_package(GLUT)
endif (OPENGL_FOUND)

if (NOT OPENGL_FOUND)
   if (MPE_FOUND)
      set (HAVE_GRAPHICS on)
      set (HAVE_MPE on)
      include_directories(${MPE_INCLUDE_DIR})
      if (X11_FOUND)
         include_directories(${X11_INCLUDE_DIR})
      endif (X11_FOUND)
      set (MPE_LIBS ${MPE_LIBRARIES})
      set (MPE_NOMPI_LIBS ${MPE_NOMPI_LIBRARIES})
      set (X11_LIBS ${X11_LIBRARIES})
   endif (MPE_FOUND)
endif (NOT OPENGL_FOUND)

########### global settings ###############
set(H_SRCS graphics.h display.h)
set(C_SRCS graphics.c display.c)
#set(F_SRCS fdisplay.F90 fgraphics.F90)

set(graphics_LIB_SRCS ${F_SRCS} ${C_SRCS} ${H_SRCS})

########### graphics target ###############

add_library(graphics STATIC ${graphics_LIB_SRCS})

set_target_properties(graphics PROPERTIES VERSION ${GRAPHICS_VERSION} SOVERSION ${GRAPHICS_VERSION_MAJOR})

set_target_properties(graphics PROPERTIES Fortran_MODULE_DIRECTORY ./CMakeFiles/graphics.dir)
target_link_libraries(graphics ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
target_link_libraries(graphics ${MPE_NOMPI_LIBS} ${X11_LIBS})
if (ImageMagick_MagickWand_FOUND)
   target_link_libraries(graphics ${ImageMagick_LIBRARIES})
endif (ImageMagick_MagickWand_FOUND)

install(TARGETS graphics DESTINATION lib)

########### pgraphics target ###############
if (MPI_FOUND)
   add_library(pgraphics STATIC ${graphics_LIB_SRCS})

   set_target_properties(pgraphics PROPERTIES VERSION ${GRAPHICS_VERSION} SOVERSION ${GRAPHICS_VERSION_MAJOR})
   set_target_properties(pgraphics PROPERTIES COMPILE_DEFINITIONS HAVE_MPI)

   set_target_properties(pgraphics PROPERTIES Fortran_MODULE_DIRECTORY ./CMakeFiles/pgraphics.dir)
   target_link_libraries(pgraphics ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
   target_link_libraries(pgraphics ${MPE_LIBS} ${X11_LIBS})
   target_link_libraries(pgraphics ${MPI_LIBRARIES})
   if (ImageMagick_MagickWand_FOUND)
      target_link_libraries(pgraphics ${ImageMagick_LIBRARIES})
   endif (ImageMagick_MagickWand_FOUND)

   install(TARGETS pgraphics DESTINATION lib)
endif (MPI_FOUND)

########### install files ###############

install(FILES  display.h DESTINATION include/graphics)
install(FILES  graphics.h DESTINATION include/graphics)

########### doc targets ################

if (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})
   add_custom_target(doc
                     COMMAND make pdfdoc
                     COMMAND make htmldoc
                     COMMAND make singlehtmldoc)
   set_target_properties(doc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/latex/Graphics.pdf
                      WORKING_DIRECTORY docs
                      COMMAND make latexpdf >& pdfdoc.out
                      DEPENDS ${PDFDOC_SRCS})
   add_custom_target(pdfdoc DEPENDS _build/latex/Graphics.pdf)
   set_target_properties(pdfdoc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/html/UsersGuide.html
                      WORKING_DIRECTORY docs
                      COMMAND make html >& htmldoc.out
                      DEPENDS ${DOC_SRCS})
   add_custom_target(htmldoc DEPENDS _build/html/UsersGuide.html)
   set_target_properties(htmldoc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/singlehtml/index.html
                      WORKING_DIRECTORY docs
                      COMMAND make singlehtml >& singlehtmldoc.out
                      DEPENDS ${DOC_SRCS})
   add_custom_target(singlehtmldoc DEPENDS _build/singlehtml/index.html)
   set_target_properties(singlehtmldoc PROPERTIES EXCLUDE_FROM_ALL TRUE)
else (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})
   add_custom_target(Graphics_doc
                     COMMAND make Graphics_pdfdoc
                     COMMAND make Graphics_htmldoc
                     COMMAND make Graphics_singlehtmldoc)
   set_target_properties(Graphics_doc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/latex/Graphics.pdf
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs
                      COMMAND make latexpdf >& pdfdoc.out
                      DEPENDS ${PDFDOC_SRCS})
   add_custom_target(Graphics_pdfdoc DEPENDS _build/latex/Graphics.pdf)
   set_target_properties(Graphics_pdfdoc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/html/UsersGuide.html
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs
                      COMMAND make html >& htmldoc.out
                      DEPENDS ${DOC_SRCS})
   add_custom_target(Graphics_htmldoc DEPENDS _build/html/UsersGuide.html)
   set_target_properties(Graphics_htmldoc PROPERTIES EXCLUDE_FROM_ALL TRUE)

   add_custom_command(OUTPUT _build/singlehtml/index.html
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs
                      COMMAND make singlehtml >& singlehtmldoc.out
                      DEPENDS ${DOC_SRCS})
   add_custom_target(Graphics_singlehtmldoc DEPENDS _build/singlehtml/index.html)
   set_target_properties(Graphics_singlehtmldoc PROPERTIES EXCLUDE_FROM_ALL TRUE)
endif (${CMAKE_PROJECT_NAME} MATCHES ${PROJECT_NAME})

########### clean files ################
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
   "install_manifest.txt;docs/_build;docs/htmldoc.out;docs/pdfdoc.out;docs/singlehtmldoc.out")

