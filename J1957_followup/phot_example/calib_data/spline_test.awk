awk '{
    if(FNR ==1) {
	ii=ii+1
    }
    if(ii ==1) {
	File=FILENAME;
	if($1=="coeff") {	 
	coeff[k]= $4;
	#print "coeff "coeff[k];
	k=k+1;
	}
	if($1=="X") {	 
	tx[l]= $6;
	#print "tx "tx[l];
	l=l+1;
	}
	if($1=="Y") {	 
	ty[m]= $6;
	#print "ty "ty[m];
	m=m+1;
	}
    }
    else if(ii ==2) {
	mag_in=$1;
	x_in=$2;
	y_in=$3;
	magcor_in=$4;
	fieldcor_in=$5;
	field = $6;
	ccd=$7
	quad=$8;
	filt=$9;

    for(i=0; i < 25; i=i+1) {
	h[i] = 0.0;
	hh[i] = 0.0;
	w_x[i] = 0.0;
	w_y[i] = 0.0;
	#print i;
    }
    #bicubic spline, 25 knots and 440 coefficients.
    nx = 24;
    ny = 24;
    #kx = 3;
    #ky = 3;
    kx = 2;
    ky = 2;

   i=0 
   j=0 
   li=0 
   lj=0 
   lx=0 
   ky1=0 
   nky1=0 
   ly=0 
   i1=0 
   j1=0 
   l2=0;
   f=0.0;
   temp=0.0;

    kx1 = kx+1;
    nkx1 = nx-kx1;
    l = kx1;
    l1 = l+1;
 while ((x_in >= tx[l1-1]) && (l != nkx1))
    {
        l = l1;
        l1 = l+1;
    }
    h[0] = 1.0;
    for (j = 1; j < kx+1; j++)
    {
        for (i = 0; i < j; i++)
        {
            hh[i] = h[i];
        }
        h[0] = 0.0;
        for (i = 0; i < j; i++)
        {
            li = l+i;
            lj = li-j;
            if (tx[li] != tx[lj])
            {
                f = hh[i] / (tx[li] - tx[lj]);
                h[i] = h[i] + f * (tx[li] - x_in);
                h[i+1] = f * (x_in - tx[lj]);
            }
            else
            {
                h[i+1-1] = 0.0;
            }
        }
    }
    lx = l-kx1;
    for (j = 0; j < kx1; j++)
    {
        w_x[j] = h[j];
    }

    ky1 = ky+1;
    nky1 = ny-ky1;
    l = ky1;
    l1 = l+1;

    while ((y_in >= ty[l1-1]) && (l != nky1))
    {
        l = l1;
        l1 = l+1;
    }
    
    h[0] = 1.0;
    for (j = 1; j < ky+1; j++)
    {
        for (i = 0; i < j; i++)
        {
            hh[i] = h[i];
        }
        h[0] = 0.0;
        for (i = 0; i < j; i++)
        {
            li = l+i;
            lj = li-j;
            if (ty[li] != ty[lj])
            {
                f = hh[i] / (ty[li] - ty[lj]);
                h[i] = h[i] + f * (ty[li] - y_in);
                h[i+1] = f * (y_in - ty[lj]);
            }
            else
            {
                h[i+1-1] = 0.0;
            }
        }
    }

    ly = l-ky1;
    for (j = 0; j < ky1; j++)
    {
        w_y[j] = h[j];
    }

    l = lx*nky1;
    for (i1 = 0; i1 < kx1; i1++)
    {
        h[i1] = w_x[i1];
    }
        
    l1 = l+ly;
    temp = 0.0;
    for (i1 = 0; i1 < kx1; i1++)
    {
        l2 = l1;
        for (j1 = 0; j1 < ky1; j1++)
        {
            l2 = l2+1;
            temp = temp + coeff[l2-1] * h[i1] * w_y[j1];
        }
        l1 = l1+nky1;
    }
    mag_out = mag_in + magcor_in + fieldcor_in + temp;
    #if(FNR == 1){
    #print "#File Mag_out Mag_in Spline_cor Bias_cor Field_cor";
    #}
    print   field, ccd, quad, filt, x_in, y_in, mag_in, mag_out, temp, magcor_in, fieldcor_in;
    }
}' $1 $2
