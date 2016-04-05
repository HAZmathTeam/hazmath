#
#************************ The X-Windows Interface ********************** 
#
#       Select an X-Windows interface: xgui or xgui_stubs.
#
#       If X-Windows is not needed, one can select xgui_stubs,
#       allowing the creation of an executable without
#       the need to load the X-Windows libraries. However,
#       the X-Windows GUI and associated graphics options
#       are completely disabled.
#

X11_INTERFACE = xgui
#X11_INTERFACE = xgui_stubs

#
#       Set the path for X-Windows include files and libraries.
#       This is not needed if xgui_stubs was selected.
#

ifeq (${X11_INTERFACE}, xgui)
	X11_INCLUDE = 
	X11_LIBRARY =  -lXm -lXt -lX11      
else
	X11_INCLUDE = 
	X11_LIBRARY =      
endif


#
#********************* The Socket Graphics Interface ******************* 
#
#       Select a Socket Graphics interface: mgvio or mgvio_stubs.
#
#       If SG is not needed, one can select mgvio_stubs,
#       allowing the creation of an executable without
#       the need to load the VIO library. However,
#       the SG and BH graphics options are completely disabled.
#
#       The VIO/SG software package is available on the home page
#       of Michael Holst:
#       http://ccom.ucsd.edu/~mholst/
#       The package provides both the SG OpenGL display tool and
#       the maloc library.
#

VIO_INTERFACE = mgvio
#VIO_INTERFACE = mgvio_stubs

#
#       Set the path for VIO include files and libraries.
#       This is not needed if mgvio_stubs was selected.
#

ifeq (${VIO_INTERFACE}, mgvio)
	VIO_INCLUDE =
	VIO_LIBRARY = -lmaloc        
else
	VIO_INCLUDE =
	VIO_LIBRARY =
endif

#
#***************************** C and FORTRAN *************************** 
#                                              
#       Define C and Fortran compilers and compile options
#

CC =  gcc -O -Wall
F90 = gfortran -O

#
#*********************************************************************** 
#
#       The rest of the Makefile is independent of environment
#

INCS =  ${X11_INCLUDE} ${VIO_INCLUDE}
LIBS =  ${X11_LIBRARY} ${VIO_LIBRARY}

HB = atest.o solver.o driver.o\
	${X11_INTERFACE}.o ${VIO_INTERFACE}.o

.f.o:
	${F90} ${INCS} -c $< -o $@

.c.o:
	${CC} ${INCS} -c $< -o $@

Hb: ${HB}
	${F90} ${HB}  ${LIBS} -o Hb

${HB} : mg0.o

clean:
	rm -f *.o fort.* core* *.mod journl*.jnl* *.out jnltmp*.jnl* Hb*
