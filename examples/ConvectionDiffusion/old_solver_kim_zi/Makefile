#
# Makefile for convection diffusion program
#
FC = gfortran
#FFLAGS = -g #-O
FFLAGS = -O3

main0 = cd0.o mg_sym.o input_cd0.o
main1 = cd1.o mgm_power.o input_cd1.o dgmres.o power_modify.o laplace.o\
        mg_sym.o krylov_depend.o
main2 = cd2.o mgm.o input_cd2.o
main3 = elliptic.o
main4 = sld.o mgm.o input_data4.o
mainkz1 = kz1.o subdom.o mgm.o input_kz1.o dgmres.o laplace.o dcg.o krylov_depend.o

tgmain = tg1.o input_tg1.o formsub.o mg_sym.o laplace.o update.o

ddmain1 = dd1.o input_dd1.o formsub.o mg_sym.o laplace.o update.o
ddmain2 = dd2.o input_dd2.o formsub.o mg_sym.o laplace.o update.o

mainsol    = solv.o the_mg.o
mainsol_m  = solv_multi.o the_mg.o
mainsol_m_g  = solv_multi_gmres.o the_mg_gmres.o dgmres.o krylov_depend.o

maintgeig  = tgeig.o input_tgeig.o

mainlap   = lap_mg.o input_lap.o

maintgmg  = tgmg.o input_tgmg.o mg_sym.o laplace.o update.o wcycle_lvl.o

maintmactn  = tgmg_actn.o input_tgmg.o lap_mg_actn.o update.o wcycle_lvl.o


OBJSlap = newtools.o\
	quad.o\
	rw.o\
	func.o

OBJStgeig = mtrx.o\
	rw.o\
	gs.o\
        assmbl.o\
        permu.o\
        transp.o\
        sgauss.o\
	mg_sym.o\
	laplace.o\
	power.o\
        quad.o\
        prolng.o\
        eltstf.o\
        eltstf_eafe.o\
        eltstf_sld.o\
        eltstf_supg.o\
        func.o\
        newtools.o\
        mattools.o

OBJS0 = rw.o\
        gs.o\
        newtools.o\
        mattools.o\
        permu.o\
        transp.o\
        sgauss.o

OBJS = mtrx.o\
       assmbl.o\
       rw.o\
       permu.o\
       transp.o\
       eltstf.o\
       eltstf_eafe.o\
       eltstf_sld.o\
       eltstf_supg.o\
       sgauss.o\
       gs.o\
       newtools.o\
       mattools.o\
       func.o\
       quad.o\
       prolng.o

OBJS1 = mtrx_sym.o\
        input_ell.o\
	mgm_sym.o\
        assmbl.o\
        rw.o\
        transp.o\
        sgauss.o\
        gs.o\
        newtools.o\
        mattools.o\
        func_sym.o\
        quad.o\
        prolng.o

kz1   : $(mainkz1) $(OBJS)
	$(FC) $(FFLAGS) -o kz1 $(mainkz1) $(OBJS) -llapack

solv     : $(mainsol) $(OBJS0)
	$(FC) $(FFLAGS) -o solv $(mainsol) $(OBJS0)

solv_multi    : $(mainsol_m) $(OBJS0)
	$(FC) $(FFLAGS) -o solv_multi $(mainsol_m) $(OBJS0)

solv_multi_g    : $(mainsol_m_g) $(OBJS0)
	$(FC) $(FFLAGS) -o solv_multi_g $(mainsol_m_g) $(OBJS0)

cd0   : $(main0) $(OBJS)
	$(FC) $(FFLAGS) -o cd0 $(main0) $(OBJS)

cd1   : $(main1) $(OBJS)
	$(FC) $(FFLAGS) -o cd1 $(main1) $(OBJS)

cd2   : $(main2) $(OBJS)
	$(FC) $(FFLAGS) -o cd2 $(main2) $(OBJS)

elliptic : $(main3) $(OBJS1)
	   $(FC) $(FFLAGS) -o elliptic $(main3) $(OBJS1)

sld   : $(main4) $(OBJS)
	$(FC) $(FFLAGS) -o sld $(main4) $(OBJS)

dd1    : $(ddmain1) $(OBJS)
	$(FC) $(FFLAGS) -o dd1 $(ddmain1) $(OBJS)

dd2    : $(ddmain2) $(OBJS)
	$(FC) $(FFLAGS) -o dd1 $(ddmain2) $(OBJS)

tg1    : $(tgmain) $(OBJS)
	$(FC) $(FFLAGS) -o tg1 $(tgmain) $(OBJS)

tgeig  : $(maintgeig) $(OBJStgeig)
	$(FC) $(FFLAGS) -o tgeig $(maintgeig) $(OBJStgeig)

lap_mg  : $(mainlap) $(OBJSlap)
	$(FC) $(FFLAGS) -o lap_mg $(mainlap) $(OBJSlap)

tgmg    : $(maintgmg) $(OBJS)
	$(FC) $(FFLAGS) -o tgmg $(maintgmg) $(OBJS)

tmactn    : $(maintmactn) $(OBJS)
	$(FC) $(FFLAGS) -o tmactn $(maintmactn) $(OBJS)


#$(main0):
#$(OBJS):


clean:
#	rm dd1 cd0 cd1 elliptic *.o
	rm *.o

