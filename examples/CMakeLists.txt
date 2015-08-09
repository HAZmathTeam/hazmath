## build examples

file( GLOB HAZMAT_EXMPLSRC RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" *.c *.f *.f90 )

##message(STATUS: "HAZMAT_EXMPLSRC=${HAZMAT_EXMPLSRC} ")

foreach( srcfile ${HAZMAT_EXMPLSRC} )
 	string( REGEX MATCH "\\.c$|\\.f90$|\\.cpp$" Z ${srcfile} )
    if("${Z}" STREQUAL ".c")
        string( REGEX REPLACE "\\.c$" "-c.ex" execfile ${srcfile} )
    elseif("${Z}" STREQUAL ".f90")
        string( REGEX REPLACE "\\.f90$" "-f.ex" execfile ${srcfile} )
    elseif("${Z}" STREQUAL ".cpp")
        string( REGEX REPLACE "\\.cpp$" "-cpp.ex" execfile ${srcfile} )
    endif()
#message(STATUS: "EXE=${executable} from ${srcfile}")
    if(execfile)
      add_executable( ${execfile} ${srcfile} )
      target_link_libraries( ${execfile} hazmat ${SUITESPARSE_LIBRARIES} ${ADD_STDLIBS} ) 
      install(TARGETS ${execfile} DESTINATION ${HAZMAT_SOURCE_DIR}/examples)
    endif(execfile)
endforeach( srcfile ${HAZMAT_EXMPLSRC} )

