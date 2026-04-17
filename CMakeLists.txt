cmake_minimum_required(VERSION 3.20)

project(pgbam VERSION 0.1.0 LANGUAGES CXX C)

include(CTest)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_POLICY_VERSION_MINIMUM 3.5)
option(PGBAM_BUILD_TESTS "Build pgbam tests" ON)

set(WITH_GOOGLETEST OFF CACHE BOOL "Disable vendored googletest in sdsl-lite" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "Disable external project tests" FORCE)

find_package(PkgConfig REQUIRED)
find_package(OpenSSL REQUIRED)
find_package(Threads REQUIRED)
find_package(OpenMP)

set(PGBAM_OPENMP_INCLUDE_DIRS)
set(PGBAM_OPENMP_COMPILE_OPTIONS)
if(NOT OpenMP_CXX_FOUND AND APPLE)
  set(PGBAM_LIBOMP_PREFIX "/opt/homebrew/opt/libomp")
  if(EXISTS "${PGBAM_LIBOMP_PREFIX}/include/omp.h" AND EXISTS "${PGBAM_LIBOMP_PREFIX}/lib/libomp.dylib")
    add_library(OpenMP::OpenMP_CXX INTERFACE IMPORTED)
    set_target_properties(OpenMP::OpenMP_CXX PROPERTIES
      INTERFACE_COMPILE_OPTIONS "-Xpreprocessor;-fopenmp"
      INTERFACE_INCLUDE_DIRECTORIES "${PGBAM_LIBOMP_PREFIX}/include"
      INTERFACE_LINK_LIBRARIES "${PGBAM_LIBOMP_PREFIX}/lib/libomp.dylib"
    )
    set(OpenMP_CXX_FOUND TRUE)
    list(APPEND PGBAM_OPENMP_INCLUDE_DIRS "${PGBAM_LIBOMP_PREFIX}/include")
    list(APPEND PGBAM_OPENMP_COMPILE_OPTIONS "-Xpreprocessor" "-fopenmp")
  endif()
endif()

pkg_check_modules(HTSLIB REQUIRED IMPORTED_TARGET htslib)
pkg_check_modules(ZSTD REQUIRED IMPORTED_TARGET libzstd)

# patchelf is only needed on Linux at install time to produce a self-contained
# package.  If it is not already on PATH we build it from source so the user
# never has to install it manually.
if(UNIX AND NOT APPLE)
  find_program(PATCHELF_EXECUTABLE patchelf)
  if(NOT PATCHELF_EXECUTABLE)
    message(STATUS "patchelf not found in PATH — will build from source")
    include(ExternalProject)
    set(_patchelf_prefix "${CMAKE_BINARY_DIR}/patchelf-host")
    ExternalProject_Add(patchelf_host
      URL      https://github.com/NixOS/patchelf/releases/download/0.18.0/patchelf-0.18.0.tar.bz2
      CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR>
      BUILD_COMMAND     make
      INSTALL_COMMAND   make install
      PREFIX            "${_patchelf_prefix}"
    )
    set(PATCHELF_EXECUTABLE "${_patchelf_prefix}/bin/patchelf")
  endif()
endif()

add_subdirectory(external/sdsl-lite EXCLUDE_FROM_ALL)
add_subdirectory(external/libhandlegraph EXCLUDE_FROM_ALL)

set(PGBAM_PARALLEL_LIBS Threads::Threads)
if(OpenMP_CXX_FOUND)
  list(APPEND PGBAM_PARALLEL_LIBS OpenMP::OpenMP_CXX)
endif()

add_library(gbwt STATIC
  external/gbwt/src/algorithms.cpp
  external/gbwt/src/bwtmerge.cpp
  external/gbwt/src/cached_gbwt.cpp
  external/gbwt/src/dynamic_gbwt.cpp
  external/gbwt/src/fast_locate.cpp
  external/gbwt/src/files.cpp
  external/gbwt/src/gbwt.cpp
  external/gbwt/src/internal.cpp
  external/gbwt/src/metadata.cpp
  external/gbwt/src/support.cpp
  external/gbwt/src/test.cpp
  external/gbwt/src/utils.cpp
  external/gbwt/src/variants.cpp
)
target_include_directories(gbwt
  PUBLIC
    external/gbwt/include
    external/sdsl-lite/include
    ${PGBAM_OPENMP_INCLUDE_DIRS}
)
target_compile_options(gbwt PRIVATE ${PGBAM_OPENMP_COMPILE_OPTIONS})
target_link_libraries(gbwt
  PUBLIC
    sdsl
    PkgConfig::ZSTD
    ${PGBAM_PARALLEL_LIBS}
)

add_library(gbwtgraph STATIC
  external/gbwtgraph/src/algorithms.cpp
  external/gbwtgraph/src/cached_gbwtgraph.cpp
  external/gbwtgraph/src/gbwtgraph.cpp
  external/gbwtgraph/src/gbz.cpp
  external/gbwtgraph/src/gfa.cpp
  external/gbwtgraph/src/index.cpp
  external/gbwtgraph/src/internal.cpp
  external/gbwtgraph/src/minimizer.cpp
  external/gbwtgraph/src/naive_graph.cpp
  external/gbwtgraph/src/path_cover.cpp
  external/gbwtgraph/src/subgraph.cpp
  external/gbwtgraph/src/utils.cpp
)
target_include_directories(gbwtgraph
  PUBLIC
    external/gbwtgraph/include
    external/gbwt/include
    external/sdsl-lite/include
    external/libhandlegraph/src/include
    ${PGBAM_OPENMP_INCLUDE_DIRS}
)
target_compile_options(gbwtgraph PRIVATE ${PGBAM_OPENMP_COMPILE_OPTIONS})
target_link_libraries(gbwtgraph
  PUBLIC
    gbwt
    handlegraph_static
    OpenSSL::Crypto
    PkgConfig::ZSTD
    ${PGBAM_PARALLEL_LIBS}
)

set(PGBAM_CORE_SOURCES
  src/annotate.cpp
  src/cli.cpp
  src/decode.cpp
  src/fingerprint.cpp
  src/gaf.cpp
  src/gbwt_backend.cpp
  src/sidecar.cpp
)

add_library(pgbam_core STATIC ${PGBAM_CORE_SOURCES})
target_include_directories(pgbam_core PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(pgbam_core
  PUBLIC
    gbwtgraph
    PkgConfig::HTSLIB
    OpenSSL::Crypto
    Threads::Threads
)

add_library(pgbam_shared SHARED ${PGBAM_CORE_SOURCES})
set_target_properties(pgbam_shared PROPERTIES
  OUTPUT_NAME pgbam
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
  EXPORT_NAME pgbam
)
target_include_directories(pgbam_shared PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(pgbam_shared
  PRIVATE
    gbwtgraph
    PkgConfig::HTSLIB
    OpenSSL::Crypto
    Threads::Threads
)
if(APPLE)
  add_custom_command(TARGET pgbam_shared POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E create_symlink
      "$<TARGET_FILE_NAME:pgbam_shared>"
      "$<TARGET_FILE_DIR:pgbam_shared>/libpgbam.so"
    VERBATIM
  )
endif()

add_executable(pgbam src/main.cpp)
target_link_libraries(pgbam PRIVATE pgbam_shared)
if(TARGET patchelf_host)
  add_dependencies(pgbam patchelf_host)
endif()
if(APPLE)
  set_target_properties(pgbam PROPERTIES INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
  set_target_properties(pgbam_shared PROPERTIES INSTALL_RPATH "@loader_path")
elseif(UNIX)
  set_target_properties(pgbam PROPERTIES INSTALL_RPATH "\$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
  set_target_properties(pgbam_shared PROPERTIES INSTALL_RPATH "\$ORIGIN")
endif()

if(PROJECT_IS_TOP_LEVEL AND PGBAM_BUILD_TESTS)
  enable_testing()
  add_executable(pgbam_tests
    tests/test_main.cpp
  )
  target_link_libraries(pgbam_tests PRIVATE pgbam_core)
  add_test(NAME pgbam_tests COMMAND pgbam_tests)
endif()

install(TARGETS pgbam_shared pgbam
  EXPORT pgbamTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Bundle non-system runtime dylibs and fix their load paths so the install
# tree is self-contained and does not point at package-manager prefixes.
#
# Both branches use $<TARGET_FILE:pgbam_shared> (the build-tree library) as
# the source for GET_RUNTIME_DEPENDENCIES.  The installed copy has already had
# its RPATH rewritten to @loader_path / $ORIGIN by the time install(CODE) runs,
# so it can no longer resolve its own deps; the build-tree copy still carries
# the original absolute paths and can.
if(APPLE)
  install(CODE "
    cmake_policy(SET CMP0011 NEW)
    cmake_policy(SET CMP0057 NEW)

    set(_lib_dir \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}\")
    set(_pgbam_lib \"\${_lib_dir}/libpgbam.dylib\")

    if(NOT EXISTS \"\${_pgbam_lib}\")
      message(FATAL_ERROR \"Bundle fixup: libpgbam.dylib not found at \${_pgbam_lib}\")
    endif()

    file(GET_RUNTIME_DEPENDENCIES
      LIBRARIES \"$<TARGET_FILE:pgbam_shared>\"
      RESOLVED_DEPENDENCIES_VAR   _resolved
      UNRESOLVED_DEPENDENCIES_VAR _unresolved
      PRE_EXCLUDE_REGEXES  \"^/usr/lib\" \"^/System/Library\"
      POST_EXCLUDE_REGEXES \"^/usr/lib\" \"^/System/Library\"
    )

    if(_unresolved)
      message(WARNING \"Bundle fixup: unresolved deps: \${_unresolved}\")
    endif()

    # Copy each dep, resolving symlinks so no dangling link ends up in the tree.
    foreach(_dep_path IN LISTS _resolved)
      get_filename_component(_dep_name \"\${_dep_path}\" NAME)
      set(_dest \"\${_lib_dir}/\${_dep_name}\")
      if(NOT EXISTS \"\${_dest}\")
        message(STATUS \"Bundling: \${_dep_name}\")
        file(REAL_PATH \"\${_dep_path}\" _real_path)
        execute_process(
          COMMAND \"${CMAKE_COMMAND}\" -E copy_if_different \"\${_real_path}\" \"\${_dest}\"
        )
        file(CHMOD \"\${_dest}\" PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ)
      endif()
      # Rewrite the bundled lib's own install name to an rpath-relative token.
      execute_process(
        COMMAND install_name_tool -id \"@rpath/\${_dep_name}\" \"\${_dest}\"
        ERROR_QUIET
      )
    endforeach()

    # Fix load-path entries in libpgbam.dylib and every bundled dep.
    set(_all_libs \"\${_pgbam_lib}\")
    foreach(_dep_path IN LISTS _resolved)
      get_filename_component(_dep_name \"\${_dep_path}\" NAME)
      list(APPEND _all_libs \"\${_lib_dir}/\${_dep_name}\")
    endforeach()

    foreach(_lib IN LISTS _all_libs)
      foreach(_dep_path IN LISTS _resolved)
        get_filename_component(_dep_name \"\${_dep_path}\" NAME)
        execute_process(
          COMMAND install_name_tool
            -change \"\${_dep_path}\" \"@loader_path/\${_dep_name}\"
            \"\${_lib}\"
          ERROR_QUIET
        )
      endforeach()
    endforeach()

    # Ensure libpgbam.dylib has @loader_path in its RPATH list.
    execute_process(
      COMMAND install_name_tool -add_rpath \"@loader_path\" \"\${_pgbam_lib}\"
      ERROR_QUIET
    )

    message(STATUS \"Runtime library bundling complete\")
  ")
elseif(UNIX)
  install(CODE "
    cmake_policy(SET CMP0011 NEW)
    cmake_policy(SET CMP0057 NEW)

    set(_patchelf \"${PATCHELF_EXECUTABLE}\")
    if(NOT EXISTS \"\${_patchelf}\")
      message(WARNING \"patchelf not available — skipping runtime library bundling.\")
      return()
    endif()

    set(_lib_dir \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}\")

    # Find the installed versioned .so (not a symlink).
    file(GLOB _pgbam_libs \"\${_lib_dir}/libpgbam.so*\")
    set(_pgbam_lib \"\")
    foreach(_f IN LISTS _pgbam_libs)
      if(NOT IS_SYMLINK \"\${_f}\")
        set(_pgbam_lib \"\${_f}\")
        break()
      endif()
    endforeach()
    if(NOT _pgbam_lib)
      message(FATAL_ERROR \"Bundle fixup: no libpgbam.so found in \${_lib_dir}\")
    endif()

    file(GET_RUNTIME_DEPENDENCIES
      LIBRARIES \"$<TARGET_FILE:pgbam_shared>\"
      RESOLVED_DEPENDENCIES_VAR   _resolved
      UNRESOLVED_DEPENDENCIES_VAR _unresolved
      PRE_EXCLUDE_REGEXES
        \"^/lib\" \"^/lib64\"
        \"^/usr/lib\" \"^/usr/lib64\"
        \"linux-vdso\" \"ld-linux\"
      POST_EXCLUDE_REGEXES
        \"^/lib\" \"^/lib64\"
        \"^/usr/lib\" \"^/usr/lib64\"
        \"linux-vdso\" \"ld-linux\"
    )

    if(_unresolved)
      message(WARNING \"Bundle fixup: unresolved deps: \${_unresolved}\")
    endif()

    # Copy each dep, resolving symlinks so no link into the source prefix ends
    # up in the install tree.
    foreach(_dep_path IN LISTS _resolved)
      get_filename_component(_dep_name \"\${_dep_path}\" NAME)
      set(_dest \"\${_lib_dir}/\${_dep_name}\")
      if(NOT EXISTS \"\${_dest}\")
        message(STATUS \"Bundling: \${_dep_name}\")
        file(REAL_PATH \"\${_dep_path}\" _real_path)
        execute_process(
          COMMAND \"${CMAKE_COMMAND}\" -E copy_if_different \"\${_real_path}\" \"\${_dest}\"
        )
        file(CHMOD \"\${_dest}\"
          PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                      GROUP_READ GROUP_EXECUTE
                      WORLD_READ WORLD_EXECUTE
        )
      endif()
      # Set \$ORIGIN RPATH on each bundled lib so its own transitive deps are
      # found in the same directory.
      execute_process(
        COMMAND \"\${_patchelf}\" --set-rpath \"\$ORIGIN\" \"\${_dest}\"
        ERROR_QUIET
      )
    endforeach()

    # libpgbam.so already has \$ORIGIN from INSTALL_RPATH, but re-stamp it to
    # be safe after any patchelf manipulation of the file.
    execute_process(
      COMMAND \"\${_patchelf}\" --set-rpath \"\$ORIGIN\" \"\${_pgbam_lib}\"
      ERROR_QUIET
    )

    message(STATUS \"Runtime library bundling complete\")
  ")
endif()

install(DIRECTORY include/pgbam
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  FILES_MATCHING PATTERN "*.hpp"
)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/pgbamConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/pgbamConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/pgbamConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pgbam
)

install(EXPORT pgbamTargets
  FILE pgbamTargets.cmake
  NAMESPACE pgbam::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pgbam
)

install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/pgbamConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/pgbamConfigVersion.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pgbam
)

if(APPLE)
  install(CODE "
    execute_process(
      COMMAND \"${CMAKE_COMMAND}\" -E create_symlink
              \"libpgbam.dylib\"
              \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/libpgbam.so\"
    )
  ")
endif()
