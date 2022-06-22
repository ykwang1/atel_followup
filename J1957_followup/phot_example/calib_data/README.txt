How to determine corrections to ZTF photometry.

Determining the three separately photometric corrections 
(for the magnitude bias, field offset, and the focal plane
spatial structure), requires the following measurement 
parameters from the ZTF pipeline:
ZTF_mag, ZTF filter, CCD number, quadrant number, RC number, CCD X, CCD Y, Field

Note, here we use the CCD and quadrant numbers as well as RC number.
The CCD and quadrant in the files could simply be swapped with 
RC numbers if desired. Also, the "g" and "r" filter bands could 
be swapped with the usual integer values (i.e. g=1,r=2).

Individual Corrections.

The ZTF magnitude bias is corrected using values from simple 
fits that depend only on the measured ZTF magnitude.

The field offset corrections are simply a look-up table.
If no value is found (for example because the fields is 
newly added) the value should be assumed to be zero and 
a "no value" error flag set.

The focal plane structure correction requires determining 
a value at a given location based on the spline fit. The
fit coefficients are in 128 separate files (1 per quadrant
for g and r). These parameters could all simply be loaded 
into a single data structure.

To see the size of the corrections for and set of paramaters
you can enter values in Input_test and then run test.awk.

Implimentation pseudo code.

	#Do the necessary steps to get the required parameters
	#(magztf, Filter, CCD, Quad, RC, X, Y)

	#Load the magnitude bias fit coefficients from file Mag_Bias.txt
	magfile=open("Mag_Bias.txt");
	magfile.read(CCD,Quad,Fitler,a,b,c);
	...
	#Get the fit coefficients for the given parameters
	a,b,c = magbias_correction(CCD,Quad,Filter);
	if(magztf < 19.25) {
	magbias_cor= a + b*magztf + c*magztf*magztf;
	}
	else { #Use 19.25 since fainter corrections are unknown.
	magbias_cor=a + b*19.25 + c*19.25*19.25;
	}

	#-------------------------------#

	#Load the field offset corrections for every quadrant 
	#from file Field_Corrections.txt
	fieldfile=open("Field_Corrections.txt")
	fieldfile.read(Field,RC,Filter,Correction)
	...
	#Get the correction for this field, RC, Filter combination.
	field_cor = field_correction(Field,RC,Fitler)

        #------------------------------#

	#Load all the spline fit files and determine the spatial structure corrections
	#For example, for CCD=07, Quadrant=2, Fitler="g"
	#the file is labelled ztf_zg_c07_q2.out.txt.spline
	splinefile=open("ztf_z"+str(Filter)+"_c"+str(CCD)+"_q"+str(Quad)+".out.txt.spline")
	splinefile.read(coeff[0-440],Xknot[0-23], Yknot[0-23])
	...
	#Get the value of the correction for the specified parameters
	spline_cor=spline_correction(X,Y,Filt,CCD,Quad);
	#Example code on how to determine spline correction value for a given 
	#CCD X and Y is given in spine.C (for C++) and spline.py (for Python)

	#------------------------------#

	#Apply all the corrections to the input ZTF mag
	magztf_new = magztf + magbias_cor + field_cor + spline_cor;
Fin.
