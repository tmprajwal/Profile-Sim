int  spin(int *ievent)       /* spin propagation; in this subroutine particle moves in X, not Z direction */
int ievent;
{

double steph,stepx,deltat,dx,dx1,dy,dy1,dz,dz1;
double bx,by,bz,b,bx1,by1,bz1,b1;
double gam=20000.;   /* Gyro of he3, in rad/sec/G */
double bav,thetaB,thetaS,adiab;
int i,k,l,m,j;
double u[3],v[3],w[3],dum[3],mod,M[3][3],Mt[3][3],S[3],S1[3],S2[3],pol;

	
	steph=2;               	/* step in magnetic field map=.2cm */
	stepx=.2;   		/* 1 mm*/
	deltat=stepx/vz; 	/*in msec*/
	x=0; y=ystart; z=zstart;
	for(i=0; i<500; i++)
	{
		
		k=x/steph;              /*  x-index in field table */
		l=(y+20)/steph;		m=(z+20)/steph;
		dx=x-k*steph; dx1=steph-dx;
		dy=y+20-l*steph; dy1=steph-dy;
		dz=z+20-m*steph; dz1=steph-dz;

				/* calculating fields - inside a space cell - 8 corners of a cube*/
		bx=hx[k][l][m]*dx1*dy1*dz1+hx[k+1][l][m]*dx*dy1*dz1+hx[k][l+1][m]*dx1*dy*dz1;
		bx=bx+hx[k+1][l+1][m]*dx*dy*dz1+hx[k][l][m+1]*dx1*dy1*dz+hx[k+1][l][m+1]*dx*dy1*dz;
		bx=(bx+hx[k][l+1][m+1]*dx1*dy*dz+hx[k+1][l+1][m+1]*dx*dy*dz)/steph/steph/steph;

		by=hy[k][l][m]*dx1*dy1*dz1+hy[k+1][l][m]*dx*dy1*dz1+hy[k][l+1][m]*dx1*dy*dz1;
		by=by+hy[k+1][l+1][m]*dx*dy*dz1+hy[k][l][m+1]*dx1*dy1*dz+hy[k+1][l][m+1]*dx*dy1*dz;
		by=(by+hy[k][l+1][m+1]*dx1*dy*dz+hy[k+1][l+1][m+1]*dx*dy*dz)/steph/steph/steph;

		bz=hz[k][l][m]*dx1*dy1*dz1+hz[k+1][l][m]*dx*dy1*dz1+hz[k][l+1][m]*dx1*dy*dz1;
		bz=bz+hz[k+1][l+1][m]*dx*dy*dz1+hz[k][l][m+1]*dx1*dy1*dz+hz[k+1][l][m+1]*dx*dy1*dz;
		bz=(bz+hz[k][l+1][m+1]*dx1*dy*dz+hz[k+1][l+1][m+1]*dx*dy*dz)/steph/steph/steph;
		b=sqrt(bx*bx+by*by+bz*bz);
		if (i<1)
		{
			S[0]=bx/b; S[1]=by/b; S[2]=bz/b;  /* at the start spin is along B */
		}

    	x=x+vz*deltat;
		y=y+vy*deltat;
		z=z+vx*deltat;
		xquad=x-40;      /* xquad=0 at the exit of the quadrupole */



			;
		   /* now the same at the end of the step */
		k=x/steph;              /*  x-index in field table */
		l=(y+20)/steph;		m=(z+20)/steph;
		dx=x-k*steph; dx1=steph-dx;
		dy=y+20-l*steph; dy1=steph-dy;
		dz=z+20-m*steph; dz1=steph-dz;

				/* calculating fields - inside a space cell - 8 corners of a cube*/
		bx1=hx[k][l][m]*dx1*dy1*dz1+hx[k+1][l][m]*dx*dy1*dz1+hx[k][l+1][m]*dx1*dy*dz1;
		bx1=bx1+hx[k+1][l+1][m]*dx*dy*dz1+hx[k][l][m+1]*dx1*dy1*dz+hx[k+1][l][m+1]*dx*dy1*dz;
		bx1=(bx1+hx[k][l+1][m+1]*dx1*dy*dz+hx[k+1][l+1][m+1]*dx*dy*dz)/steph/steph/steph;

		by1=hy[k][l][m]*dx1*dy1*dz1+hy[k+1][l][m]*dx*dy1*dz1+hy[k][l+1][m]*dx1*dy*dz1;
		by1=by1+hy[k+1][l+1][m]*dx*dy*dz1+hy[k][l][m+1]*dx1*dy1*dz+hy[k+1][l][m+1]*dx*dy1*dz;
		by1=(by1+hy[k][l+1][m+1]*dx1*dy*dz+hy[k+1][l+1][m+1]*dx*dy*dz)/steph/steph/steph;

		bz1=hz[k][l][m]*dx1*dy1*dz1+hz[k+1][l][m]*dx*dy1*dz1+hz[k][l+1][m]*dx1*dy*dz1;
		bz1=bz1+hz[k+1][l+1][m]*dx*dy*dz1+hz[k][l][m+1]*dx1*dy1*dz+hz[k+1][l][m+1]*dx*dy1*dz;
		bz1=(bz1+hz[k][l+1][m+1]*dx1*dy*dz+hz[k+1][l+1][m+1]*dx*dy*dz)/steph/steph/steph;
		b1=sqrt(bx1*bx1+by1*by1+bz1*bz1);

		thetaB=acos((bx*bx1+by*by1+bz*bz1)/b/b1);   /* angle of B field rotation, rad */
		bav=(b+b1)/2;
		thetaS=deltat/1000*gam*bav;        /* angle of spin rotation, rad */
		adiab=thetaB/thetaS;
			/* spin evaluation */
		u[0]=bx/b; u[1]=by/b; u[2]=bz/b;     /* making field vectors unitary */
		dum[0]=bx1/b1; dum[1]=by1/b1; dum[2]=bz1/b1;
		/* new coordinate system u,v,w: U - along B, B1 iz in UV plane, W is ortogonal to B and B1 */
		vectorproduct(u,dum,w);		/*finding w*/
	    mod=sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]); 
	    w[0]=w[0]/mod; w[1]=w[1]/mod; w[2]=w[2]/mod;   /*scaling to unitary vector w*/
	    vectorproduct(w,u,v);                          /*finding v*/
		for(j=0; j<3; j++)
	  {
		  M[j][0]=u[j];  M[j][1]=v[j];  M[j][2]=w[j];       /* rotation matrix */
		  Mt[0][j]=u[j];   Mt[1][j]=v[j];   Mt[2][j]=w[j];   /* inverse rotation matrix */
	  }
		
		rotate(S,M,S1);  /* translate Spin into uvw system*/
		S2[0]=S1[0];
		S2[1]=S1[1]*cos(thetaS)+S1[2]*sin(thetaS);  /* spin rotation in B */
		S2[2]=S1[2]*cos(thetaS)-S1[1]*sin(thetaS);	
		rotate(S2,Mt,S);  /* return back to global coordinate system */
		pol=u[0]*S[0]+u[1]*S[1]+u[2]*S[2];     /* define polarization as angle betwwen S and B */
		polmem[ievent]=pol*100.;
		
/*fprintf(out, "%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f \n",x,y,z,vx,vy,vz,bx,by,bz,adiab,pol); */


	/*	fprintf(out, "%f, %f, %f,%f,%f, %f, %f,%f, %f, %f, %f, %f, %f, %f \n",x,y,z,adiab,pol,b,bx,by,bz,S[0],S[1],S[2],thetaB,thetaS);*/

	}

	
/*		fprintf(out, "%f, %f, %f,%f,%f, %f, %f,%f, %f, %f, %f, %f, %f, %f \n",x,y,z,adiab,pol,b,bx,by,bz,S[0],S[1],S[2],thetaB,thetaS);*/


/*		if(pol>0.6)
		{
			fclose(out);
			if ((out = fopen("test.csv", "wb")) == NULL)
			{  printf( "Cannot open output file \n"); }
			return 0;
		}
		fclose(out);
		return 1;
		


		fprintf(out, "%f, %f, %f ,%f, %f, %f \n",x,y,z,ymem,zmem,pol);*/

return 0;
}
