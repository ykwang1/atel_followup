awk '{
    if(FNR ==1) {
	i=i+1;
	#if(i==1) {
	#magcor=0;
	#fieldcor=0;
	#}
    }
    if(i== 1 && FNR > 1) {
    #Check and load the input test file values 
    err=0;

    if(NF < 7 || NF > 7) {
	print "***** Error with number of inputs ****";
	err=1;
    }
    if($1 < 10 || $1 > 22) {
        print "***** Error with input magnitude "$1" ****";
        err=1;
    }
    if($2 != "g" && $2 != "r") {
        print "***** Error with Filter ID "$2" ****";
        err=1;
    }
    if(length($7) != 6 ) {
	print "***** Error with Field number "$7" ****";
	err=1;
    }    
    ccdi = int($3);
    if(ccdi < 0 || ccdi > 16 || length($3) < 2) {
	print "***** Error using CCD value "$3" ****";
	err=1;
    }
    if($4 < 1 || $4 > 4) {
	print "***** Error using Quadrant number "$4" ****";
	err=1;
    }
    if($5 < 0 || $5 > 3100) {
	print "***** Error with X coord "$5" ****";
	err=1;
    }
    if($6 < 0 || $6 > 3100) {
	print "***** Error with Y coord "$6" ****";
	err=1;
    }
    if(err > 0) {
	print "on input line "FNR"\n\n\n"
        exit;
    }

    m=m+1;
    #ZTF mag
    mag_inp[m]= $1;
    #ZTF filter
    filt_inp[m]=$2;
    #two digit CCD ID
    CCD_inp[m]=$3;
    #Quadrant number
    Quad_inp[m]=$4;
    #Find the RC number
    RC_inp[m] = (CCD_inp[m]-1)*4 + Quad_inp[m] -1;
    #CCD X coordinate
    X_inp[m]=$5;
    #CCD Y coordinate
    Y_inp[m]=$6;
    #Field six digit code
    Field_inp[m]=$7;
    #Arrays to hold the corrections
    magcor[m] =0;
    fieldcor[m]=0;
    }
    else if(i==2 && FNR > 1) {
	#Read the magnitude bias fit parameters
	j=j+1
	ccd[j] = $1;
	quad[j] = $2;
	filt[j] = $3;
	#coeffs
	a[j] = $4;
	b[j] = $5;
	c[j] = $6;
	for(l=1; l < m+1; l=l+1) {
	    if(CCD_inp[l] == ccd[j] && Quad_inp[l] == quad[j] && filt_inp[l] == filt[j]) {
	     magcor[l]= a[j] + b[j]*mag_inp[l] + c[j]*mag_inp[l]*mag_inp[l]; 
	     #break;
	    }
	}
    }
    else if(i==3 && FNR > 1) {
	#Read the field offset values
	k=k+1
	field[k] =$1;
	RC[k] = $2;
        Filt_field[k]= $3;
	Cor[k] = $4;
	for(l=1; l < m+1; l=l+1) {
	    if(Field_inp[l] == field[k] && RC_inp[l]==RC[k] && filt_inp[l] == Filt_field[k]) {
		fieldcor[l] = Cor[k];
		#break;
	    }
	}
    }
    else if(i ==4) {
	#Select the correct spline fit file.
	print "#Field CCD Quad Filter X Y Mag_in Mag_out Spline_cor Bias_cor Field_cor";
	filenum = 0;
	last_quad=0;
	for(l=1; l < m+1; l=l+1) {
	    if(Quad_inp[l] != last_quad) {#Different quadrant than previous
		filenum = filenum+1;
		splinefile ="ztf_z"filt_inp[l]"_c"CCD_inp[l]"_q"Quad_inp[l]".out.txt.spline"
		splineinp = "spline_inp"filenum;
		last_quad=Quad_inp[l];
		print "./spline_test.awk "splinefile" "splineinp >> "toexe";
	    }
	    print mag_inp[l], X_inp[l], Y_inp[l], magcor[l], fieldcor[l], Field_inp[l], CCD_inp[l], Quad_inp[l], filt_inp[l] >> splineinp;
	}
    }

}' Input_test Mag_Bias.txt Field_Corrections.txt tmp
#Find the values of the spline correction.
chmod u+x toexe
./toexe
#Remove the junk
rm spline_inp*
rm toexe

#The input test file parameters are
#ZTF mag, ZTF filter, CCD number, quadrant number, CCD X, CCD Y, Field
