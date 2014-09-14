#  -*- mode: cmake -*-

set(TPL_MISSING FALSE)

if(USE_GSL AND NOT GSL_FOUND)
set(TPL_MISSING TRUE)
endif()
