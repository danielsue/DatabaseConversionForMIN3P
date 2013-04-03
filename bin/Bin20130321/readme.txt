So as to support the more than 5 products in reaction, the output file for complex.dbs and gases.dbs are minor changed. Please modify the following source codes and use the new compiled executable file. 


readsspc.F90
	101  format(6x,i1,3x,5(a12,1x,f7.3,1x))
to
	101  format(6x,i2,2x,20(a12,1x,f7.3,1x))

readgses.f90
	100  format(a12,2x,2f10.4,31x,f9.4/6x,i1,3x,5(a12,1x,f7.3,1x))
to
	100  format(a12,2x,2f10.4,31x,f9.4/6x,i2,2x,20(a12,1x,f7.3,1x))