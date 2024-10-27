# generated from ament/cmake/core/templates/nameConfig.cmake.in

# prevent multiple inclusion
if(_nmpcsine_CONFIG_INCLUDED)
  # ensure to keep the found flag the same
  if(NOT DEFINED nmpcsine_FOUND)
    # explicitly set it to FALSE, otherwise CMake will set it to TRUE
    set(nmpcsine_FOUND FALSE)
  elseif(NOT nmpcsine_FOUND)
    # use separate condition to avoid uninitialized variable warning
    set(nmpcsine_FOUND FALSE)
  endif()
  return()
endif()
set(_nmpcsine_CONFIG_INCLUDED TRUE)

# output package information
if(NOT nmpcsine_FIND_QUIETLY)
  message(STATUS "Found nmpcsine: 1.0.0 (${nmpcsine_DIR})")
endif()

# warn when using a deprecated package
if(NOT "" STREQUAL "")
  set(_msg "Package 'nmpcsine' is deprecated")
  # append custom deprecation text if available
  if(NOT "" STREQUAL "TRUE")
    set(_msg "${_msg} ()")
  endif()
  # optionally quiet the deprecation message
  if(NOT ${nmpcsine_DEPRECATED_QUIET})
    message(DEPRECATION "${_msg}")
  endif()
endif()

# flag package as ament-based to distinguish it after being find_package()-ed
set(nmpcsine_FOUND_AMENT_PACKAGE TRUE)

# include all config extra files
set(_extras "ament_cmake_export_include_directories-extras.cmake")
foreach(_extra ${_extras})
  include("${nmpcsine_DIR}/${_extra}")
endforeach()
