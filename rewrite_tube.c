
#include   <stdio.h>
#include   <stdlib.h>
#include   <math.h>
//double pi=3.14159;

double randgen(double min,double max)   /* generates random number between min and max */
//double min,max;
{
double res;
int iran;

  iran=rand();
  res=min+double(max-min)*iran/RAND_MAX;
  return res;
}

double func(double x,double y,double pi, double r, double vp, double vz)   /* gaussian distribution */
//double x,y;
{
double A, sigma, r2,res;

//  sigma=3.5;
	sigma=3;
  A=1.0/2./pi/sigma/sigma;
  r2=x*x+y*y;
//  res=A*exp(-r2/2/sigma/sigma);
//  res=exp(-sqrt(r2));
//	res=abs(x/r);
	res=vp/vz;
  return res;
}

double maxfunc(double x,double y,double r,double pi, double vp, double vz)   /* finds approx. max of a function within a circle at x,y with radius r  */
//double x,y,r;
{
double step,res,max,x1,y1;
int i,j;

  step=r/20.;
  max=0;
  for (i=0; i<40; i++)   
 {
     for (j=0; j<40; j++)   
     {
		 x1=x-r+step*i;
		 y1=y-r+step*j;
         res=func(x1,y1,pi,r,vp,vz);
		 if(res>max) max=res;
	 }
  
 }
  res=max*1.2;  /* to insure that max is really max) */
  return res;
}

double integr(double x0,double y0,double r,double pi)   /* integration over the circle */
//double x0,y0,r;
{
double max,total,x1,y1,z1,c1,c2,res,vx,vy,vz,vp,vr;
int i,ntot,ngood,n;

	double *xi= new double[10000000];
	double *yi= new double[10000000];
	double *zi= new double[10000000];
	double *xqstore= new double[10000000];
	double *yqstore= new double[10000000];
	double *zqstore= new double[10000000];
	double *xustore= new double[10000000];
	double *yustore= new double[10000000];
	double *zustore= new double[10000000];
	double *xdstore= new double[10000000];
	double *ydstore= new double[10000000];
	double *zdstore= new double[10000000];
	double *vqx= new double[10000000];
	double *vqy= new double[10000000];
	double *vqz= new double[10000000];
	double *vux= new double[10000000];
	double *vuy= new double[10000000];
	double *vuz= new double[10000000];
	double *vdx= new double[10000000];
	double *vdy= new double[10000000];
	double *vdz= new double[10000000];
	double *pol= new double[10000000];
	
	int polcount=0;

//  max=maxfunc(x0,y0,r,pi,vp,vz);    /* max of the function in the circle*/
	max=1;
  total=pi*r*r*max;
	ntot=0; ngood=0;
//*
	ifstream input;
	input.open("profiles_PB_small_step_VLN_theta0p1.txt");

	while(!input.eof())
	{
		input>>xqstore[polcount]>>yqstore[polcount]>>zqstore[polcount]>>vqx[polcount]>>vqy[polcount]>>vqz[polcount]>>pol[polcount]>>xustore[polcount]>>yustore[polcount]>>zustore[polcount]>>vux[polcount]>>vuy[polcount]>>vuz[polcount]>>pol[polcount]>>xdstore[polcount]>>ydstore[polcount]>>zdstore[polcount]>>vdx[polcount]>>vdy[polcount]>>vdz[polcount]>>pol[polcount];
/*
		xstore[polcount]=xqstore[polcount];
		ystore[polcount]=yqstore[polcount];
		zstore[polcount]=zqstore[polcount];
//*/
/*
		xstore[polcount]=xustore[polcount];
		ystore[polcount]=yustore[polcount];
		zstore[polcount]=zustore[polcount];
//*/
///*
		x1=xdstore[polcount];
		y1=ydstore[polcount];
		z1=zdstore[polcount];
		vx=vdx[polcount];
		vy=vdy[polcount];
		vz=vdz[polcount];
//*/
//*
		vp=sqrt(vx*vx+vy*vy);
		vr=sqrt(vx*vx+vy*vy+vz*vz);
		 if((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)<r*r)
		 {
			 ntot++;
/*
//			 c1=randgen(0.0,max);
			c1=randgen(0.0,0.04);
			c1=randgen(0.0,12e-15);
//			 c2=func(x1,y1,pi,r);
			c2=func(x1,y1,pi,r,vp,vz);
/*
			if(x0!=0)
			{c2=double(abs(x0)/982.0);}
			if(x0==0)
			{c2=0.4;}
/*
			if(vz!=0)
			{c2=abs(vx[polcount])/vz[polcount];}
			if(vz==0)
			{c2=1;}
//*/
//*
//			 if(c2>c1) ngood++;
		 }
		polcount++;
}

input.close();
//*/
//  max=maxfunc(x0,y0,r,pi);    /* max of the function in the circle*/
//  total=pi*r*r*max;      /* total square */
/*  n=1000000; ntot=0; ngood=0;
  for (i=0; i<n; i++)
	 {
//		 x1=randgen(x0-r,x0+r);
//		 y1=randgen(y0-r,y0+r);

		 x1=randgen(0-1e-8,1e-8);
		 y1=randgen(0-1e-8,1e-8);
		 if((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)<r*r)
		 {
			 ntot++;
//			 c1=randgen(0.0,max);
//			 c2=func(x1,y1,pi,r);
//			 if(c2>c1) ngood++;
		 }
	 }
//*/  
  
//  res=total*ngood/ntot;
	res=ntot;
//cout<<ntot<<" "<<ngood<<" "<<double(abs(x0)/980.0)<<endl;
  return res;
}
		


void rewrite_tube()
//int argc;
//char **argv;
{
double  x,y,r,aaa=2,bbb=2;
int i;

double pi=3.14159;

FILE *out;


     
    if ((out = fopen("dum1.csv", "wb")) == NULL)
    {  printf( "Cannot open output file \n"); } 
  
	  srand( (unsigned)time( NULL ) ); /* initializing random generator */
      rand();
	  x=0; y=0; r=3.5;
	  for (i=0; i<80; i++)
	  {
		     x=-40+i;
			 aaa=integr(x,y,r,pi)/5.44;
//			 bbb=func(x,y,pi,r)*1.;
			 printf(" x=%f aaa=%f  bbb=%f \n   ",x,aaa,bbb);
			 fprintf(out," %f  %f  %f  %f \n   ",x,y,aaa,bbb);
	  }
	  
	  fclose(out);
/*   scanf("%lf",aaa);*/
   return;


    
}


  
