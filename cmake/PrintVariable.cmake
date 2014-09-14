# -*- mode: cmake -*-

function(PRINT_VAR VAR_NAME)
    message("=> " "${VAR_NAME}=${${VAR_NAME}}")
endfunction()    

