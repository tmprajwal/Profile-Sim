
#include   <stdio.h>
#include   <stdlib.h>
#include   <math.h>
double pi=3.14159;

double randgen(min,max)   /* generates random number between min and max */
double min,max;
{
double res;
int iran;

  iran=rand();
  res=min+(max-min)*iran/32767.;
  return res;
}

double func(x,y)   /* gaussian distribution */
double x,y;
{
double A, sigma, r2,res;

  sigma=3.5;
  A=1.0/2./pi/sigma/sigma;
  r2=x*x+y*y;
  res=A*exp(-r2/2/sigma/sigma);
  res=exp(-sqrt(r2));
  return res;
}

double maxfunc(x,y,r)   /* finds approx. max of a function within a circle at x,y with radius r  */
double x,y,r;
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
         res=func(x1,y1);
		 if(res>max) max=res;
	 }
  
 }
  res=max*1.2;  /* to insure that max is really max) */
  return res;
}

double integr(x0,y0,r)   /* integration over the circle */
double x0,y0,r;
{
double max,total,x1,y1,c1,c2,res;
int i,ntot,ngood,n;

  max=maxfunc(x0,y0,r);    /* max of the function in the circle*/
  total=pi*r*r*max;      /* total square */
  n=1000000; ntot=0; ngood=0;
  for (i=0; i<n; i++)
	 {
		 x1=randgen(x0-r,x0+r);
		 y1=randgen(y0-r,y0+r);
		 if((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)<r*r)
		 {
			 ntot++;
			 c1=randgen(0.0,max);
			 c2=func(x1,y1);
			 if(c2>c1) ngood++;
		 }
	 }
  
  
  res=total*ngood/ntot;
  return res;
}
		


void main(argc,argv)
int argc;
char **argv;
{
double  x,y,r,aaa=2,bbb;
int i;

FILE *out;


     
    if ((out = fopen("dum1.csv", "wb")) == NULL)
    {  printf( "Cannot open output file \n"); } 
  
	  srand( (unsigned)time( NULL ) ); /* initializing random generator */
      rand();
	  x=0; y=0; r=3.5;
	  for (i=0; i<80; i++)
	  {
		     x=-40+i;
			 aaa=integr(x,y,r)/5.44;
			 bbb=func(x,y)*1.;
			 printf(" x=%f aaa=%f  bbb=%f \n   ",x,aaa,bbb);
			 fprintf(out," %f, %f,  %f \n   ",x,aaa,bbb);
	  }
	  
	  fclose(out);
/*   scanf("%lf",aaa);*/
   return;


    
}


  