### CMakeLists.txt; coding: utf-8 ---

FILE(GLOB modules "${CMAKE_CURRENT_SOURCE_DIR}/*.cmake")
INSTALL(FILES ${modules} DESTINATION share/cmake/modules COMPONENT Devel)
