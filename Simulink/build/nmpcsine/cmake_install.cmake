# Install script for directory: C:/Users/mpilo/Documents/Work/Research Project/Repo/EEE4022S-Project-2024-NXSMPI001/Simulink/src/nmpcsine

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "C:/Users/mpilo/Documents/Work/Research Project/Repo/EEE4022S-Project-2024-NXSMPI001/Simulink/install")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/nmpcsine" TYPE EXECUTABLE FILES "C:/Users/mpilo/Documents/Work/Research Project/Repo/EEE4022S-Project-2024-NXSMPI001/Simulink/build/nmpcsine/nmpcSine.exe")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/ament_index/resource_index/package_run_dependencies" TYPE FILE FILES "C:/Users/mpilo/Documents/Work/Research Project/Repo/EEE4022S-Project-2024-NXSMPI001/Simulink/build/nmpcsine/ament_cmake_index/share/ament_index/resource_index/package_run_dependencies/nmpcsine")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/ament_index/resource_index/parent_prefix_path" TYPE FILE FILES "C:/Users/mpilo/Documents/Work/Research Project/Repo/EEE4022S-Project-2024-NXSMPI001/Simulink/build/nmpcsine/ament_cmake_index/share/ament_index/resource_index/parent_prefix_path/nmpcsine")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/nmpcsine/environment" TYPE FILE FILES "C:/Program Files/MATLAB/R2024b/sys/ros2/win64/ros2/share/ament_cmake_core/cmake/environment_hooks/environment/ament_prefix_path.bat")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/nmpcsine/environment" TYPE FILE FILES "C:/Users/mpilo/Documents/Work/Research Project/Repo/EEE4022S-Project-2024-NXSMPI001/Simulink/build/nmpcsine/ament_cmake_environment_hooks/ament_prefix_path.dsv")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/nmpcsine/environment" TYPE FILE FILES "C:/Program Files/MATLAB/R2024b/sys/ros2/win64/ros2/share/ament_cmake_core/cmake/environment_hooks/environment/path.bat")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/nmpcsine/environment" TYPE FILE FILES "C:/Users/mpilo/Documents/Work/Research Project/Repo/EEE4022S-Project-2024-NXSMPI001/Simulink/build/nmpcsine/ament_cmake_environment_hooks/path.dsv")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/nmpcsine" TYPE FILE FILES "C:/Users/mpilo/Documents/Work/Research Project/Repo/EEE4022S-Project-2024-NXSMPI001/Simulink/build/nmpcsine/ament_cmake_environment_hooks/local_setup.bat")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/nmpcsine" TYPE FILE FILES "C:/Users/mpilo/Documents/Work/Research Project/Repo/EEE4022S-Project-2024-NXSMPI001/Simulink/build/nmpcsine/ament_cmake_environment_hooks/local_setup.dsv")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/nmpcsine" TYPE FILE FILES "C:/Users/mpilo/Documents/Work/Research Project/Repo/EEE4022S-Project-2024-NXSMPI001/Simulink/build/nmpcsine/ament_cmake_environment_hooks/package.dsv")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/ament_index/resource_index/packages" TYPE FILE FILES "C:/Users/mpilo/Documents/Work/Research Project/Repo/EEE4022S-Project-2024-NXSMPI001/Simulink/build/nmpcsine/ament_cmake_index/share/ament_index/resource_index/packages/nmpcsine")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/nmpcsine/cmake" TYPE FILE FILES "C:/Users/mpilo/Documents/Work/Research Project/Repo/EEE4022S-Project-2024-NXSMPI001/Simulink/build/nmpcsine/ament_cmake_export_include_directories/ament_cmake_export_include_directories-extras.cmake")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/nmpcsine/cmake" TYPE FILE FILES
    "C:/Users/mpilo/Documents/Work/Research Project/Repo/EEE4022S-Project-2024-NXSMPI001/Simulink/build/nmpcsine/ament_cmake_core/nmpcsineConfig.cmake"
    "C:/Users/mpilo/Documents/Work/Research Project/Repo/EEE4022S-Project-2024-NXSMPI001/Simulink/build/nmpcsine/ament_cmake_core/nmpcsineConfig-version.cmake"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/nmpcsine" TYPE FILE FILES "C:/Users/mpilo/Documents/Work/Research Project/Repo/EEE4022S-Project-2024-NXSMPI001/Simulink/src/nmpcsine/package.xml")
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "C:/Users/mpilo/Documents/Work/Research Project/Repo/EEE4022S-Project-2024-NXSMPI001/Simulink/build/nmpcsine/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
