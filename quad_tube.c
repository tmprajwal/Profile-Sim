//This code is designed to simulate helium-3 atoms passing through the quadrupole tube as a background gas.
//The background helium-3 gas will have thermalized to room temperature after bouncing off of cryostat components.
//


#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<algorithm>

#define pi 3.14159		//Pi, you know from math and circles and such
#define c 3.0e8			//Speed of light (rounded) in m/sec
#define kbol 8.62e-11		//Boltzman consant, MeV/K
#define mhe3 2809.4		//He3 mass, MeV/c^2
#define mkg 5.008237e-27	//He3 mass, kg
#define mu 1.07e-26		//He3 magnetic moment, J/T
#define tnoz 1.0		//Nozzle temperature, K
#define ntot 1000000		//Number of particles
//#define ntot 100000000
//#define ntot 50
//#define ntot 100
//#define ntot 100000
#define stept 3e-6		//Time steps, in seconds
#define gradq 112.5		//Gradient in quads, T/m (9 kG on the tip)



//Since some of these functions call each other (in there definitions). It is a good idea to declare each function before defining them by writing out all of the headers. A function can be called as long as it has been declared. If it hasn't been defined yet, the code will search for the definition. This is the easiest way to ensure that no errors occur due to a function being called before it is declared. MB
double getnew(double rnoz, double xnoz, double ynoz, double znoz, volatile double& x, volatile double& y, volatile double& z, volatile double& vx, volatile double& vy, volatile double& vz, double& xmem, double& ymem, double& zmem);

double ranr(double r0);

double randgen(double min, double max);

double ransincos(double max);

double ranmax();

int ranpol();

int movetube(int bounce, double xtube, double ytube, double ztube, double rtube, double ltube, double znoz, volatile double vx, volatile double vy, volatile double vz, volatile double& x, volatile double& y, volatile double& z, int& dterror, int& dtrerror, double sbounce, int magnumber);

double newtheta(volatile double& vx, volatile double& vy, volatile double& vz, volatile double& x, volatile double& y);

int movefree(volatile double& x, volatile double& y, volatile double& z, volatile double vx, volatile double vy, volatile double vz, double zRGA, double ztube, double ltube);




//double getnew(int& itest, double rnoz, double xnoz, double ynoz, double znoz, double& x, double& y, double& z, double& vx, double& vy, double& vz, int& npol, double& xmem, double& ymem, double& zmem, double raper, double zaper)	/* generates a new particle at(0,0,0) with random direction */
//Actually on the z coordinate is guaranteed to be zero. x and y are randomized to a certain point on the cross-sectional area of the nozzle. MB
/*{
	double theta,phi,thetamax;
	double v,r;
//	Values needed by this function


//	cout<<"one"<<endl;
	r=ranr(rnoz);		//Calls the ranr() function (defined in code) and sets "r" equal to the return value
//	cout<<"two"<<" "<<r<<endl;
	phi=randgen(0.,2*pi);	//Calls the randgen() function (defined in code) and sets phi to the return value. MB
//	cout<<"three"<<" "<<phi<<endl;
	x=r*sin(phi)+xnoz;	//Sets starting value for x at some random position on the nozzle. Just a conversion from spherical to cartesian coordinates and then adding the x coordinate location of the nozzle. MB
	y=r*cos(phi)+ynoz;	//Same thing as x, but for y. MB
	z=znoz;			//"z" is the direction which the nozzle is facing. All particles start at the same location in z. MB

	npol=ranpol();		//Uses ranpol() to randomly assign 1 or -1 to npol. MB
//	npol=-1;
//	cout<<"four"<<" "<<npol<<endl;
	return atan2(sqrt(vx*vx+vy*vy),vz);
}
*/


double ranr(double r0)		/* generates random number between 0 and r0 with y=r/r0 distribution */
{
double r, comp;		//variables defined just for use in generating this random number. Still need to double check that the distribution matches what Genya expected. MB
    do
	{
	r=randgen(0.,r0);
	comp=randgen(0.,r0);
//	The two lines above each generate a random number according to the randgen function (defined elsewhere in code) with a minimum value of 0 and a maximum value of r0. MB
	} while(comp>r);
//	The "do" loop is performed until r is equal to or greater than comp. MB
	return r;	//function returns the value of r. MB
}



double randgen(double min, double max)	/* generates random number between min and max */
{
double res;	//Variable just for use in this function

	res=min+(double(max)-min)*rand()/RAND_MAX;
//	Uses predefined rand() function and picks a random number between min and max. "res" is set to this value. "RAND_MAX" is the maximum value of the rand() function. MB
	return res;	//function returns the value "res" to the main program. MB
}



double ransincos(double max)	/* generates random number between 0 and max with sin*cos distribution */
{
	double theta, comp;	//Variables just for this function
	do
	{
		theta=randgen(0.,max);
//		Obtain random value for theta using randgen(). MB
		comp=randgen(0.,0.5);
//		Obtain random value for comp using randgen(). MB
	}while(comp>sin(theta)*cos(theta));
//	Perform do loop until while condition is met. Still need to double check that this random distribution is as expected. MB
	return theta;	//returns the value theta
}



double ranmax()		/* generates random number between 0 and 5*vmax with Maxwell distribution */
{
	double dum, vmax, v, pmax, pcur, comp, v1, v2, fvmax;
//	Variables just for this function
	dum=mhe3/(2*kbol*tnoz*c*c);
//	cout<<dum<<endl;
	v1=sqrt(1/dum);
//	cout<<v1<<endl;
//	v2=200;			//Scattering. Genya
	v2=0.001;		//No scattering. Genya
//	Scattering on residual Gas.
//	cout<<v2<<endl;
	vmax=v1*sqrt(0.75+sqrt(9/16+v2*v2/(v1*v1)));
//	cout<<vmax<<" "<<5*vmax<<endl;
//	maximum-probabilty velocity. Genya
	pmax=vmax*vmax*vmax*exp(-vmax*vmax/(v1*v1)-v2*v2/(vmax*vmax));
//	fvmax=5*vmax;
	do
	{
		v=randgen(0.,5*vmax);
		pcur=v*v*v*exp(-v*v/(v1*v1)-v2*v2/(v*v))/pmax;
		comp=randgen(0.,1.);
	}while(comp>pcur);

	return v;	
}



int movetube(int bounce, double xtube, double ytube, double ztube, double rtube, double ltube, double znoz, volatile double vx, volatile double vy, volatile double vz, volatile double& x, volatile double& y, volatile double& z, int& dterror, int& dtrerror, double sbounce, int magnumber)		/* moving the particle through skimmer */
{
	double dtr,dtr1,dtr2,dtz,xdum,ydum,zdum,rdum,vr;	//Some function parameters. MB
	xdum=x;					//X position before tube. MB
	ydum=y;					//Y position before tube. MB
	zdum=z;					//Z position before tube. MB
	rdum=sqrt(x*x+y*y);
	
	vr=sqrt(vx*vx+vy*vy);
	dtz=0;
	dtr=0;
	
	if(bounce>=0)
	{
		if(vr!=0)
		{

			dtr1=(sqrt(vr*vr*((rtube*rtube)-(rdum*rdum))+(x*vx+y*vy)*(x*vx+y*vy))-x*vx-y*vy)/(vr*vr);
			dtr2=(0-sqrt(vr*vr*((rtube*rtube)-(rdum*rdum))+pow(x*vx+y*vy,2))-x*vx-y*vy)/(vr*vr);
			
			if(dtr1<1e-5 && dtr1>0 && dtr2<1e-5 && dtr2>0)
			{
				cout<<"Weirdness:	x: "<<x<<"	y: "<<y<<"	z: "<<z<<"	vx: "<<vx<<"	vy: "<<vy<<"	vz: "<<vz<<endl<<"	R: "<<sqrt(x*x+y*y)<<"	bounce: "<<bounce<<endl;
			}
/*			
			if(dtr1<1e-7 && dtr1>=0)
			{dtr1=0;}
		
			if(dtr2<1e-7 && dtr2>=0)
			{dtr2=0;}
*/
//			cout<<dtr1<<" "<<dtr2<<endl;
		
			if(dtr1<1e-5 && dtr2<1e-5)
			{
				dtrerror++;
				cout<<"first "<<magnumber<<" "<<sbounce<<" "<<bounce<<" "<<rdum<<"	"<<x*vx<<"	"<<y*vy<<"	"<<dtr1<<"	"<<dtr2<<endl;
				return 3;
			}
		
			if(dtr1<=0 && dtr2>0)
			{
				dtr=dtr2;
			}
		
			if(dtr2<=0 && dtr1>0)
			{
				dtr=dtr1;
			}
		
			if(dtr1>0 && dtr2>0)
			{
				if(dtr1<dtr2)
				{
					dtr=dtr1;
				}
			
				if(dtr2<dtr1)
				{
					dtr=dtr2;
				}
			
			}
		
		}
		
	}


	if(vz>0)
	{
		dtz=(ztube+ltube-zdum)/vz;	//Time it takes particle to traverse the tube/channel (assuming it doesn't hit an obstacle). MB
	}
	
	if(vz<0)
	{
		dtz=(ztube-zdum)/vz;
	}
	
//	cout<<dtr<<" "<<dtz<<endl;
	
	
	if((dtr<dtz && dtr!=0) || (dtz<=0 && dtr>0))
	{
//		cout<<x<<" "<<y<<" "<<z<<endl;
		x=x+vx*dtr;				//X position after movetube.
		y=y+vy*dtr;				//Y position after movetube.
		z=z+vz*dtr;				//Z position after movetube.
//		cout<<x<<" "<<y<<" "<<z<<endl;
		return 0;
	}
	
	else if((dtz<=dtr && dtz!=0) || (dtr<=0 && dtz>0))
	{
		x=x+vx*dtz;				//X position after movetube.
		y=y+vy*dtz;				//Y position after movetube.
		z=z+vz*dtz;				//Z position after movetube.
		return -1;
	}
	
	if(dtz==0 && dtr==0)
	{
		dterror++;
		return 2;
	}
	
	return 10;
}



int movefree(volatile double& x, volatile double& y, volatile double& z, volatile double vx, volatile double vy, volatile double vz, double zRGA, double ztube, double ltube)
//Move particle freely for a small distance after exiting the quad. Checks natural spread of particles. MB
{
	double dt,xdum,ydum;	//Some function parameters. MB
	z=zRGA;		//Setting starting z position for particle. MB
	dt=(zRGA-ztube-ltube)/vz;	//Time it takes particle to traverse the skimmer (assuming it doesn't hit an obstacle). MB
	xdum=x;			//X position before skimmer. MB
	ydum=y;			//Y position before skimmer. MB
	x=x+vx*dt;		//X position after skimmer (assuming it hasn't hit an obstacle). MB
	y=y+vy*dt;		//Y position after skimmer (assuming it hasn't hit an obstacle). MB
	return 1;
}



double newtheta(volatile double& vx, volatile double& vy, volatile double& vz, volatile double& x, volatile double& y)
{
	double theta,phi,phiprime,thetamax;
	double v=sqrt(vx*vx+vy*vy+vz*vz);
	
	phi=atan2(y,x);
	
	thetamax=pi/2;	//Thetamax is pi/2.
	theta=ransincos(thetamax);
	phiprime=randgen(0.,2*pi);

	vx=v*(cos(phi)*cos(pi-theta)-sin(phi)*cos(phiprime)*sin(pi-theta));
	vy=v*(sin(phi)*cos(pi-theta)+cos(phi)*cos(phiprime)*sin(pi-theta));
	vz=v*sin(phiprime)*sin(pi-theta);
	
//	cout<<theta*180/pi<<"	"<<phiprime*180/pi<<"	"<<phi*180/pi<<"	"<<atan(vy/vx)*180/pi<<endl;
	
	return theta;
}



int moveholeby(int bounce, double xtube, double ytube, double ztube, double rtube, double ltube, double znoz, volatile double vx, volatile double vy, volatile double vz, volatile double& x, volatile double& y, volatile double& z, int& dterror, int& dtrerror, double rhole, double lhole, int magnet, double zcent)
{
	double dtr,dtr1,dtr2,dty,xdum,ydum,zdum,rdum,vr;	//Some function parameters. MB
	xdum=x;					//X position before tube. MB
	ydum=y;					//Y position before tube. MB
	zdum=z;					//Z position before tube. MB
	rdum=sqrt(x*x+(z-zcent)*(z-zcent));
	
	double phasex=0;
	double phasey=0;
	double phasez=0;
	vr=sqrt(vx*vx+vz*vz);
	dty=0;
	dtr=0;
	
	if(vr!=0)
	{

		dtr1=(sqrt(vr*vr*((rhole*rhole)-(rdum*rdum))+(x*vx+(z-zcent)*vz)*(x*vx+(z-zcent)*vz))-x*vx-(z-zcent)*vz)/(vr*vr);
		dtr2=(0-sqrt(vr*vr*((rhole*rhole)-(rdum*rdum))+pow(x*vx+(z-zcent)*vz,2))-x*vx-(z-zcent)*vz)/(vr*vr);
		
		if(dtr1<1e-5 && dtr1>0 && dtr2<1e-5 && dtr2>0)
		{
			cout<<"Weirdness:	x: "<<x<<"	y: "<<y<<"	z: "<<z-zcent<<"	vx: "<<vx<<"	vy: "<<vy<<"	vz: "<<vz<<endl<<"	R: "<<sqrt(x*x+y*y)<<"	r: "<<sqrt(y*y+(z-zcent)*(z-zcent))<<"	bounce: "<<bounce<<endl;
		}
/*		
		if(dtr1<1e-5 && dtr1>=0)
		{
			dtr1=0;
//			cout<<dtr1<<" "<<dtr2<<endl;
		}
		
		if(dtr2<1e-5 && dtr2>=0)
		{
			dtr2=0;
//			cout<<dtr1<<" "<<dtr2<<endl;
		}
*/
//		cout<<dtr1<<" "<<dtr2<<endl;
	
		if(dtr1<0 && dtr2<0)
		{
			if(magnet==1)
			{return -3;}
			
			dtrerror++;
			cout<<"second "<<bounce<<" "<<rdum<<"	"<<x*vx<<"	"<<(z-zcent)*vz<<"	"<<dtr1<<"	"<<dtr2<<endl;
			return 3;
		}
		
		if(dtr1<=1e-5 && dtr2>1e-5)
		{
			dtr=dtr2;
		}
		
		if(dtr2<=1e-5 && dtr1>1e-5)
		{
			dtr=dtr1;
		}
		
		if((dtr1>1e-5 && dtr2>1e-5) || (dtr1<=1e-5 && dtr1>0 && dtr2<=1e-5 && dtr2>0))
		{
			if(dtr1<dtr2)
			{
				dtr=dtr1;
			}
			
			if(dtr2<dtr1)
			{
				dtr=dtr2;
			}
			
		}
		
	}
	
	
	
	if(y>0)
	{
		if(vy>0)
		{
			phasex=x+vx*dtr;
//			dty=(rtube+lhole-ydum)/vy;	//Time it takes particle to traverse the tube/channel (assuming it doesn't hit an obstacle). MB
			dty=((rtube+lhole)*cos(asin((phasex-0)/(rtube+lhole)))-ydum)/vy;
			if(dty<=dtr)
			{
				dty=(rtube+lhole-ydum)/vy;
				
				if(dty<0)
				{
					cout<<"ERROR"<<"Past Magnet:	x: "<<x<<"	y: "<<y<<"	z: "<<z-zcent<<"	vx: "<<vx<<"	vy: "<<vy<<"	vz: "<<vz<<endl<<"	R: "<<sqrt(x*x+y*y)<<"	r: "<<sqrt(y*y+(z-zcent)*(z-zcent))<<"	bounce: "<<bounce<<endl;
					return 4;
				}
				
				x=x+vx*dty;				//X position after movetube.
				y=y+vy*dty;				//Y position after movetube.
				z=z+vz*dty;				//Z position after movetube.
				return 1;
			}
			
			phasex=0;
			
		}
	
		if(vy<0)
		{
			phasex=x+vx*dtr;
//			dty=(rtube*cos(asin((phasex-0)/rtube))+lhole-ydum)/vy;
			dty=(rtube*cos(asin((phasex-0)/rtube))-ydum)/vy;
			if(magnet==1)
			{
				if(dtr<dty)
				{
					phasey=y+vy*dtr;
					if(phasey>(rtube+lhole)*cos(asin((phasex-0)/(rtube+lhole))))
					{
						return -3;
					}
					
					phasey=0;
				}
				
			}
			
			phasex=0;
		}
		
	}
	
	
	if(y<0)
	{
		if(vy<0)
		{
			phasex=x+vx*dtr;
			dty=(0-rtube-lhole-ydum)/vy;	//Time it takes particle to traverse the tube/channel (assuming it doesn't hit an obstacle). MB
			dty=(0-(rtube+lhole)*cos(asin((phasex-0)/(rtube+lhole)))-ydum)/vy;
			if(dty<=dtr)
			{
				dty=(0-rtube-lhole-ydum)/vy;
				x=x+vx*dty;				//X position after movetube.
				y=y+vy*dty;				//Y position after movetube.
				z=z+vz*dty;				//Z position after movetube.
				return 1;
			}
			
			phasex=0;
			
		}
	
		if(vy>0)
		{
			phasex=x+vx*dtr;
//			dty=(0-rtube*cos(asin((phasex-0)/rtube))-lhole-ydum)/vy;
			dty=(0-rtube*cos(asin((phasex-0)/rtube))-ydum)/vy;
			if(magnet==1)
			{
				if(dtr<dty)
				{
					phasey=y+vy*dtr;
					if(phasey<0-(rtube+lhole)*cos(asin((phasex-0)/(rtube+lhole))))
					{
						return -3;
					}
					
					phasey=0;
				}
				
			}
			
			phasex=0;
		}
		
	}
	
	
	if((dtr<dty && dtr!=0) || (dty==0 && dtr>0))
	{
//		cout<<dtr<<" "<<x<<" "<<y<<" "<<z-zcent<<endl;
		x=x+vx*dtr;				//X position after movetube.
		y=y+vy*dtr;				//Y position after movetube.
		z=z+vz*dtr;				//Z position after movetube.
//		cout<<x<<" "<<y<<" "<<z-zcent<<endl;
		return 0;
	}
	
	else if((dty<=dtr && dty!=0) || (dtr==0 && dty!=0))
	{
/*		
		if(y>0)
		{
			dty=(sqrt(rtube*rtube-rhole*rhole)-ydum)/vy;
		}
		
		if(y<0)
		{
			dty=(0-sqrt(rtube*rtube-rhole*rhole)-ydum)/vy;
		}
		
		x=x+vx*dty;				//X position after movetube.
		y=y+vy*dty;				//Y position after movetube.
		z=z+vz*dty;				//Z position after movetube.
*/		
		return -2;
	}
	
	if(dty==0 && dtr==0)
	{
		dterror++;
		return 2;
	}
	
	return 10;
		
}



int moveholebx(int bounce, double xtube, double ytube, double ztube, double rtube, double ltube, double znoz, volatile double vx, volatile double vy, volatile double vz, volatile double& x, volatile double& y, volatile double& z, int& dterror, int& dtrerror, double rhole, double lhole, int magnet, double zcent)
{
	double dtr,dtr1,dtr2,dtx,xdum,ydum,zdum,rdum,vr;	//Some function parameters. MB
	xdum=x;					//X position before tube. MB
	ydum=y;					//Y position before tube. MB
	zdum=z;					//Z position before tube. MB
	rdum=sqrt(y*y+(z-zcent)*(z-zcent));
	
	double phasex=0;
	double phasey=0;
	double phasez=0;
	vr=sqrt(vy*vy+vz*vz);
	dtx=0;
	dtr=0;
	
	if(vr!=0)
	{

		dtr1=(sqrt(vr*vr*((rhole*rhole)-(rdum*rdum))+(y*vy+(z-zcent)*vz)*(y*vy+(z-zcent)*vz))-y*vy-(z-zcent)*vz)/(vr*vr);
		dtr2=(0-sqrt(vr*vr*((rhole*rhole)-(rdum*rdum))+pow(y*vy+(z-zcent)*vz,2))-y*vy-(z-zcent)*vz)/(vr*vr);
		
		if(dtr1<1e-5 && dtr1>0 && dtr2<1e-5 && dtr2>0)
		{
			cout<<"Weirdness:	x: "<<x<<"	y: "<<y<<"	z: "<<z-zcent<<"	vx: "<<vx<<"	vy: "<<vy<<"	vz: "<<vz<<endl<<"	R: "<<sqrt(x*x+y*y)<<"	r: "<<sqrt(y*y+(z-zcent)*(z-zcent))<<"	bounce: "<<bounce<<endl;
		}
/*		
		if(dtr1<1e-5 && dtr1>=0)
		{dtr1=0;}
		
		if(dtr2<1e-5 && dtr2>=0)
		{dtr2=0;}
*/
//		cout<<dtr1<<" "<<dtr2<<endl;
	
		if(dtr1<0 && dtr2<0)
		{
			if(magnet==1)
			{return -3;}
			
			dtrerror++;
			cout<<"third "<<bounce<<" "<<rdum<<"	"<<x*vx<<"	"<<z*vz<<"	"<<dtr1<<"	"<<dtr2<<endl;
			return 3;
		}
		
		if(dtr1<=1e-5 && dtr2>1e-5)
		{
			dtr=dtr2;
		}
		
		if(dtr2<=1e-5 && dtr1>1e-5)
		{
			dtr=dtr1;
		}
		
		if((dtr1>1e-5 && dtr2>1e-5) || (dtr1<=1e-5 && dtr1>0 && dtr2<=1e-5 && dtr2>0))
		{
			if(dtr1<dtr2)
			{
				dtr=dtr1;
			}
			
			if(dtr2<dtr1)
			{
				dtr=dtr2;
			}
			
		}
		
	}
	
	
	
	if(x>0)
	{
		if(vx>0)
		{
			phasey=y+vy*dtr;
//			dtx=(rtube+lhole-xdum)/vx;	//Time it takes particle to traverse the tube/channel (assuming it doesn't hit an obstacle). MB
			dtx=((rtube+lhole)*cos(asin((phasey-0)/(rtube+lhole)))-xdum)/vx;
			if(dtx<=dtr)
			{
//				dtx=(rtube+lhole-xdum)/vx;
				x=x+vx*dtx;				//X position after movetube.
				y=y+vy*dtx;				//Y position after movetube.
				z=z+vz*dtx;				//Z position after movetube.
				return 1;
			}
			
			phasey=0;
			
		}
	
		if(vx<0)
		{
			phasey=y+vy*dtr;
//			dtx=(rtube*cos(asin((phasey-0)/rtube))+lhole-xdum)/vx;
			dtx=(rtube*cos(asin((phasey-0)/rtube))-xdum)/vx;
			if(magnet==1)
			{
				if(dtr<dtx)
				{
					phasex=x+vx*dtr;
					if(phasex>(rtube+lhole)*cos(asin((phasey-0)/(rtube+lhole))))
					{
						return -3;
					}
					
					phasex=0;
				}
				
			}
			
			phasey=0;
		}
		
	}
	
	
	if(x<0)
	{
		if(vx<0)
		{
			phasey=y+vy*dtr;
//			dtx=(0-rtube-lhole-xdum)/vx;	//Time it takes particle to traverse the tube/channel (assuming it doesn't hit an obstacle). MB
			dtx=(0-(rtube+lhole)*cos(asin((phasey-0)/(rtube+lhole)))-xdum)/vx;
			if(dtx<=dtr)
			{
//				dtx=(0-rtube-lhole-xdum)/vx;
				x=x+vx*dtx;				//X position after movetube.
				y=y+vy*dtx;				//Y position after movetube.
				z=z+vz*dtx;				//Z position after movetube.
				return 1;
			}
			
			phasey=0;
			
		}
	
		if(vx>0)
		{
			phasey=y+vy*dtr;
//			dtx=(0-rtube*cos(asin((phasey-0)/rtube))-lhole-xdum)/vx;
			dtx=(0-rtube*cos(asin((phasey-0)/rtube))-xdum)/vx;
			if(magnet==1)
			{
				if(dtr<dtx)
				{
					phasex=x+vx*dtr;
					if(phasex<0-(rtube+lhole)*cos(asin((phasey-0)/(rtube+lhole))))
					{
						return -3;
					}
					
					phasex=0;
				}
				
			}
			
			phasey=0;
		}
		
	}
	
	
	if((dtr<dtx && dtr!=0) || (dtx==0 && dtr>0))
	{
//		cout<<dtr<<" "<<x<<" "<<y<<" "<<z-zcent<<endl;
		x=x+vx*dtr;				//X position after movetube.
		y=y+vy*dtr;				//Y position after movetube.
		z=z+vz*dtr;				//Z position after movetube.
//		cout<<x<<" "<<y<<" "<<z-zcent<<endl;
		return 0;
	}
	
	else if((dtx<=dtr && dtx!=0) || (dtr==0 && dtx!=0))
	{
/*
		if(x>0)
		{
			dtx=(sqrt(rtube*rtube-rhole*rhole)-xdum)/vx;
		}
		
		if(x<0)
		{
			dtx=(0-sqrt(rtube*rtube-rhole*rhole)-xdum)/vx;
		}
		
		x=x+vx*dtx;				//X position after movetube.
		y=y+vy*dtx;				//Y position after movetube.
		z=z+vz*dtx;				//Z position after movetube.
*/
		return -2;
	}
	
	if(dtx==0 && dtr==0)
	{
		dterror++;
		return 2;
	}
	
	return 10;
		
}



double newthetabmy(volatile double& vx, volatile double& vy, volatile double& vz, volatile double& x, volatile double& y)
{
//	double theta,phi,phiprime,thetamax;
	double theta,phiprime,thetamax;
	double v=sqrt(vx*vx+vy*vy+vz*vz);
	
//	phi=atan2(y,x);
	
	thetamax=pi/2;	//Thetamax is pi/2.
	theta=ransincos(thetamax);
	phiprime=randgen(0.,2*pi);
	
	if(y>0)
	{
		vx=v*sin(theta)*sin(phiprime);
		vy=0-v*cos(theta);
		vz=v*sin(theta)*cos(phiprime);
	}
	
	if(y<0)
	{
		vx=v*sin(theta)*sin(phiprime);
		vy=v*cos(theta);
		vz=v*sin(theta)*cos(phiprime);
	}
	
//	cout<<theta*180/pi<<"	"<<phiprime*180/pi<<"	"<<phi*180/pi<<"	"<<atan(vy/vx)*180/pi<<endl;
	
	return theta;
}



double newthetabmx(volatile double& vx, volatile double& vy, volatile double& vz, volatile double& x, volatile double& y)
{
//	double theta,phi,phiprime,thetamax;
	double theta,phiprime,thetamax;
	double v=sqrt(vx*vx+vy*vy+vz*vz);
	
//	phi=atan2(y,x);
	
	thetamax=pi/2;	//Thetamax is pi/2.
	theta=ransincos(thetamax);
	phiprime=randgen(0.,2*pi);
	
	if(x>0)
	{
		vx=0-v*cos(theta);
		vy=v*sin(theta)*cos(phiprime);
		vz=v*sin(theta)*sin(phiprime);
	}
	
	if(x<0)
	{
		vx=v*cos(theta);
		vy=v*sin(theta)*cos(phiprime);
		vz=v*sin(theta)*sin(phiprime);
	}
	
//	cout<<theta*180/pi<<"	"<<phiprime*180/pi<<"	"<<phi*180/pi<<"	"<<atan(vy/vx)*180/pi<<endl;
	
	return theta;
}



double newthetabhy(volatile double& vx, volatile double& vy, volatile double& vz, volatile double& x, volatile double& y, volatile double& z, double zcent)
{
	double theta,phi,phiprime,thetamax;
	double v=sqrt(vx*vx+vy*vy+vz*vz);
	
	phi=atan2(x,z-zcent);
//	phi=atan2(z,x);
	
	thetamax=pi/2;	//Thetamax is pi/2.
	theta=ransincos(thetamax);
	phiprime=randgen(0.,2*pi);

	vx=v*(sin(phi)*cos(pi-theta)+cos(phi)*cos(phiprime)*sin(pi-theta));
	vy=v*sin(phiprime)*sin(pi-theta);
	vz=v*(cos(phi)*cos(pi-theta)-sin(phi)*cos(phiprime)*sin(pi-theta));
/*	
	vx=v*(sin(phi)*cos(pi-theta)+cos(phi)*cos(phiprime)*sin(pi-theta));
	vy=v*sin(pi-theta)*sin(phi);
	vz=v*(cos(phi)*cos(pi-theta)-sin(phi)*cos(phiprime)*sin(pi-theta));
*/	
/*
	phi=atan2(y,x);
	
	thetamax=pi/2;	//Thetamax is pi/2.
	theta=ransincos(thetamax);
	phiprime=randgen(0.,2*pi);

	vx=v*(cos(phi)*cos(pi-theta)-sin(phi)*cos(phiprime)*sin(pi-theta));
	vy=v*(sin(phi)*cos(pi-theta)+cos(phi)*cos(phiprime)*sin(pi-theta));
	vz=v*sin(phiprime)*sin(pi-theta);
/*/
//	cout<<theta*180/pi<<"	"<<phiprime*180/pi<<"	"<<phi*180/pi<<"	"<<atan(vy/vx)*180/pi<<endl;
	
	return theta;
}



double newthetabhx(volatile double& vx, volatile double& vy, volatile double& vz, volatile double& x, volatile double& y, volatile double& z, double zcent)
{
	double theta,phi,phiprime,thetamax;
	double v=sqrt(vx*vx+vy*vy+vz*vz);
	
	phi=atan2(z-zcent,y);
	
	thetamax=pi/2;	//Thetamax is pi/2.
	theta=ransincos(thetamax);
	phiprime=randgen(0.,2*pi);


	vx=v*sin(phiprime)*sin(pi-theta);
	vy=v*(cos(phi)*cos(pi-theta)-sin(phi)*cos(phiprime)*sin(pi-theta));
	vz=v*(sin(phi)*cos(pi-theta)+cos(phi)*cos(phiprime)*sin(pi-theta));
/*	
	vx=v*sin(phiprime)*sin(pi-theta);
	vy=v*(cos(phi)*cos(pi-theta)-sin(phi)*cos(phiprime)*sin(pi-theta));
	vz=v*(sin(phi)*cos(pi-theta)+cos(phi)*cos(phiprime)*sin(pi-theta));
*/
/*	//Original NewTheta Code start.
	phi=atan2(y,x);
	
	thetamax=pi/2;	//Thetamax is pi/2.
	theta=ransincos(thetamax);
	phiprime=randgen(0.,2*pi);

	vx=v*(cos(phi)*cos(pi-theta)-sin(phi)*cos(phiprime)*sin(pi-theta));
	vy=v*(sin(phi)*cos(pi-theta)+cos(phi)*cos(phiprime)*sin(pi-theta));
	vz=v*sin(phiprime)*sin(pi-theta);
	//Original NewTheta Code stop.
*/
//	cout<<theta*180/pi<<"	"<<phiprime*180/pi<<"	"<<phi*180/pi<<"	"<<atan(vy/vx)*180/pi<<endl;
	
	return theta;
}



void quad_tube()
{
	volatile double x, y, z, vx, vy, vz;
//	particle position and velocity in cartesian coordinates
	double xnoz, ynoz, znoz, rnoz, xtube, ytube, ztube, rtube, ltube;
//	nozzle and tube positions (in cartesian coordinates) and radii.
	double zRGA;
//	Inclusion by me for the movefree function. MB
	int polcount=0;
//	Iterator for polarization array. MB
	double xmem, ymem, zmem;
//	x and y velocity of particle. zmem is the polar speed of the particle. Three are saved here to be compared with end velocities.
	double rhole, lhole, alphamax, hgap, lsec, zcent;
	double hoffset=10.9;
	double zoffset;
	double thetay, dtr1, dtr2, vr, rdum;
	bool mag=false;
	int magnumber=0;

	int ktot, kpass, kfail, kpump, khole;
	int sc=0;
	int numy=0;
	int dterror=0;
	int dtrerror=0;
	int flag, check, radcheck;
	int sbounce=0;
	int magnet=0;
	int sfail=0;
	int spass=0;
	
	int scl=2200792;
	scl=250832;
	
	int *bounce= new int[ntot];
	int *hbounce= new int[ntot];

	Double_t *xstore= new Double_t[ntot];
	Double_t *ystore= new Double_t[ntot];
	Double_t *zstore= new Double_t[ntot];
	Double_t *vxstore= new Double_t[ntot];
	Double_t *vystore= new Double_t[ntot];
	Double_t *vzstore= new Double_t[ntot];
	Double_t *exittheta= new Double_t[ntot];
	Double_t *exitphi= new Double_t[ntot];
	
	Double_t *vxsc= new Double_t[5000000];
	Double_t *vysc= new Double_t[5000000];
	Double_t *vzsc= new Double_t[5000000];
	Double_t *thetasc= new Double_t[5000000];
	Double_t *phisc= new Double_t[5000000];
	Double_t *xsc= new Double_t[5000000];
	Double_t *ysc= new Double_t[5000000];
	Double_t *zsc= new Double_t[5000000];
	
	double exitvrho;
//	I want to store the position of all the particles that make it through the quad. This way I can reconstruct the predicted profile. MB
	
	
	xnoz=0.0;   ynoz=0.0;   znoz=0.0;	rnoz=6.337/2;						//Initial position for channel.
	xtube=0.0;  ytube=0.0;  ztube=0.0;  rtube=6.3;		ltube=1288;		//Changed ltube to be an even multiple of the section length.
	rhole=2.5;	lhole=1;	alphamax=0.7545;
	hgap=12.4;	lsec=161;

	zRGA=ztube+ltube+40;

/*	
	rnoz=5;
	
	rnoz=0.0508/2; rtube=rnoz; ltube=5.08;		//Ghost Ranch Dimensions
*/	

	ltube=1295/2;

//	rtube=rnoz;	

	ofstream output;
	output.open("qt_profile.txt");
	//Create variable to hold output file and create/open/rewrite an output file.
	
	ofstream errorout;
	errorout.open("qt_error_file.txt");
	
	ofstream track;
	track.open("qt_tracker.txt");
	
	ifstream input;
	input.open("qt_initial_profile_1.txt");
	
	ofstream bout;
	bout.open("qt_bounce.txt");
	
	while(!input.eof())
	{
		input>>xsc[sc]>>ysc[sc]>>zsc[sc]>>vxsc[sc]>>vysc[sc]>>vzsc[sc];
		
		thetasc[sc]=atan2(sqrt(vxsc[sc]*vxsc[sc]+vysc[sc]*vysc[sc]),vzsc[sc]);
		phisc[sc]=atan2(vysc[sc],vxsc[sc]);
		sc++;
	}

	input.close();
	
	
	srand((unsigned)time(NULL));	//Uses internal clock to fetch seed for random number generator.
	rand();
	
	double avthetas=0;
	double avthetaf=0;
	double avthetad=0;
	
	
	
	ktot=0; kpass=0; kfail=0; kpump=0; khole=0;
		
	for(int i=0; i<ntot; i++)
	{
//		avthetas+=getnew(rnoz,xnoz,ynoz,znoz,x,y,z,vx,vy,vz,xmem,ymem,zmem);
		bounce[i]=0;
		flag=0;
		
		double thetamax=0.0506;
//		thetamax=pi/2;
			
		do
		{
			x=xsc[numy];
			y=ysc[numy];
			z=0;
			vx=vxsc[numy]/1000;
			vy=vysc[numy]/1000;
			vz=vzsc[numy]/1000;
			thetay=thetasc[numy];
			radcheck=1;
			if(sqrt(x*x+y*y)>=rtube-0.000001*rtube)
			{
				radcheck=0;
			}
			
			numy++;
			if(numy>scl)
			{numy=0;}
		}while(thetay>thetamax && radcheck==0);
		
//		cout<<"NEW: "<<sqrt(x*x+y*y)<<endl;
		
//		track<<endl<<i<<" "<<x<<" "<<y<<" "<<z;

		magnumber=0;
		
		do
		{			
			flag=movetube(bounce[i],xtube,ytube,ztube,rtube,ltube,znoz,vx,vy,vz,x,y,z,dterror,dtrerror,sbounce,magnumber);
			if(mag==true)
			{
				sbounce++;
			}
//			cout<<flag<<endl;
			
//			track<<endl<<i<<" "<<x<<" "<<y<<" "<<z;

			if(flag==10)
			{
				cout<<"MAJOR ERROR Initial"<<endl;
				break;
			}
			
			if(flag==0-1 && z>=ztube+ltube-0.000001*ltube)
			{
				exitvrho=sqrt(vx*vx+vy*vy);
				exittheta[kpass]=atan2(exitvrho,vz);
				exitphi[kpass]=atan2(vy,vx);
//*				
				output<<endl<<i<<" "<<x<<" "<<y<<" "<<z<<" "<<vx<<" "<<vy<<" "<<vz<<" "<<exittheta[kpass]<<" "<<exitphi[kpass]<<" "<<bounce[i];
				
				xstore[kpass]=x;
				ystore[kpass]=y;
				zstore[kpass]=z;
				vxstore[kpass]=vx;
				vystore[kpass]=vy;
				vzstore[kpass]=vz;
//				cout<<sqrt(x*x+y*y)<<"	"<<x<<"	"<<y<<"	"<<z<<" "<<exittheta[kpass]*180/pi<<" "<<exitphi[kpass]*180/pi<<"	"<<bounce[i]<<" "<<kpass<<endl;
//*/				
				kpass++;
			}
			
			if(z<=ztube+0.00000000001 && bounce[i]>0)
			{
				flag=0-2;
				kfail++;
//				cout<<kfail<<"	"<<kpass<<"	"<<bounce[i]<<endl;
			}
			
			if(flag>1)
			{
				errorout<<endl<<"Error: "<<dterror+dtrerror<<":"<<endl;
				if(flag==2)
				{
					errorout<<"Infinite Velocity: "<<dterror<<":"<<endl<<i<<" "<<x<<" "<<y<<" "<<z<<" "<<vx<<" "<<vy<<" "<<vz<<" "<<bounce[i];	
				}
				
				if(flag==3)
				{
					errorout<<"Negative Time: "<<dtrerror<<":"<<endl<<i<<" "<<x<<" "<<y<<" "<<z<<" "<<vx<<" "<<vy<<" "<<vz<<" "<<bounce[i];
					cout<<"Hit Side of Hole Wall:	x: "<<x<<"	y: "<<y<<"	z: "<<z-zcent<<"	vx: "<<vx<<"	vy: "<<vy<<"	vz: "<<vz<<endl<<"	R: "<<sqrt(x*x+y*y)<<"	r1: "<<sqrt(x*x+(z-zcent)*(z-zcent))<<"	r2: "<<sqrt(y*y+(z-zcent)*(z-zcent))<<"	bounce: "<<bounce[i]<<endl;
				}
			
				cout<<"Error"<<endl;
				break;
			}
			
			if(sqrt(x*x+y*y)>rtube+0.000001*rtube)
			{
				dterror++;
				errorout<<endl<<"Error: "<<dterror+dtrerror<<":"<<endl;
				errorout<<"Exceeded tube radius: "<<dterror<<":"<<endl<<i<<" "<<x<<" "<<y<<" "<<z<<" "<<vx<<" "<<vy<<" "<<vz<<" "<<bounce[i]<<" "<<hbounce[i]<<" "<<check;
				cout<<"fourth "<<sqrt(x*x+y*y)<<"	"<<rtube<<endl;
				break;
			}
			
			if((z>=0 && z<lsec) || (z>=lsec && z<2*lsec) || (z>=2*lsec && z<3*lsec) || (z>=3*lsec && z<4*lsec) || (z>=4*lsec && z<5*lsec) || (z>=6*lsec && z<7*lsec) || (z>=7*lsec && z<8*lsec))
			{
				zoffset=0;
				for(int j=0;j<8;j++)
				{
					zcent=double(int(z/lsec)*lsec)+hoffset+zoffset+j*hgap+rhole;
					if(z>=zcent-rhole && z<=zcent+rhole)
					{
						if((y>=rtube*cos(asin((x-0)/rtube)) && y<=rtube) || (y>=0-rtube && y<=0-rtube*cos(asin((x-0)/rtube))))
						{
							if(x>=0-rhole && x<=rhole)
							{
								if(sqrt(x*x+pow(z-zcent,2))<=rhole)
								{
									khole++;
									hbounce[i]=0;
									magnumber++;
									do{
										check=moveholeby(hbounce[i],xtube,ytube,ztube,rtube,ltube,znoz,vx,vy,vz,x,y,z,dterror,dtrerror,rhole,lhole,magnet,zcent);
										hbounce[i]=hbounce[i]+1;
										magnet=0;
										mag=true;
										sbounce=0;
										if(check==10)
										{
											cout<<"MAJOR ERROR"<<endl<<endl<<endl<<endl<<endl<<endl;
											flag=0-7;
											break;
										}
										
										if(check==-3)
										{
											flag=0-5;
											kpump++;
											kfail++;
//											cout<<"Pumped"<<endl;
											break;
										}
										
										if(check==-2)
										{
											flag=4;
											flag=movetube(bounce[i],xtube,ytube,ztube,rtube,ltube,znoz,vx,vy,vz,x,y,z,dterror,dtrerror,sbounce,magnumber);
//											cout<<"Exiting Pump Hole:	vx: "<<vx<<"	vy: "<<vy<<"	vz: "<<vz<<"	r: "<<sqrt(x*x+y*y)<<endl;
//											cout<<"Exiting Pump Hole"<<endl;
//											cout<<"Exiting pump hole:	x: "<<x<<"	y: "<<y<<"	z: "<<z-zcent<<"	vx: "<<vx<<"	vy: "<<vy<<"	vz: "<<vz<<endl<<"	R: "<<sqrt(x*x+y*y)<<"	r: "<<sqrt(x*x+(z-zcent)*(z-zcent))<<"	bounce: "<<bounce[i]<<endl;
										}
										
										if(check==1)
										{
											avthetad=newthetabmy(vx,vy,vz,x,y);
											bounce[i]=bounce[i]+1;
											magnet=1;
//											mag=true;
//											cout<<"Hit Magnet"<<endl;
										}
										
										if(check>1)
										{
											errorout<<endl<<"Hole Error: "<<dterror+dtrerror<<":"<<endl;
											if(check==2)
											{
												errorout<<"Infinite Velocity: "<<dterror<<":"<<endl<<i<<" "<<x<<" "<<y<<" "<<z-zcent<<" "<<vx<<" "<<vy<<" "<<vz<<" "<<bounce[i]<<" "<<hbounce[i];	
											}
				
											if(check==3)
											{
												errorout<<"Negative Time: "<<dtrerror<<":"<<endl<<i<<" "<<x<<" "<<y<<" "<<z-zcent<<" "<<vx<<" "<<vy<<" "<<vz<<" "<<bounce[i]<<" "<<hbounce[i];
											}
											
											if(check==4)
											{
												errorout<<"Past Magnet: "<<dtrerror<<":"<<endl<<i<<" "<<x<<" "<<y<<" "<<z-zcent<<" "<<vx<<" "<<vy<<" "<<vz<<" "<<bounce[i]<<" "<<hbounce[i];
											}
											
											flag=0-6;
											break;
										}
										
										if(check==0)
										{
											avthetad=newthetabhy(vx,vy,vz,x,y,z,zcent);
											bounce[i]=bounce[i]+1;
/*
											rdum=sqrt(x*x+(z-zcent)*(z-zcent));
											vr=sqrt(vx*vx+vz*vz);
											dtr1=(sqrt(vr*vr*((rhole*rhole)-(rdum*rdum))+(x*vx+(z-zcent)*vz)*(x*vx+(z-zcent)*vz))-x*vx-(z-zcent)*vz)/(vr*vr);
											dtr2=(0-sqrt(vr*vr*((rhole*rhole)-(rdum*rdum))+pow(x*vx+(z-zcent)*vz,2))-x*vx-(z-zcent)*vz)/(vr*vr);
											cout<<dtr1<<" "<<dtr2<<endl;
//*/
//											cout<<"Hit Side of Hole Wall:	x: "<<x<<"	y: "<<y<<"	z: "<<z-zcent<<"	vx: "<<vx<<"	vy: "<<vy<<"	vz: "<<vz<<endl<<"	R: "<<sqrt(x*x+y*y)<<"	r: "<<sqrt(x*x+(z-zcent)*(z-zcent))<<"	bounce: "<<bounce[i]<<endl;
//											cout<<"Hit Side of Hole Wall"<<endl;
										}
										
									}while(check>=0);
//									cout<<"broken"<<endl<<endl;
									break;
									
								}
								
							}
							
						}
						
						if((x>=rtube*cos(asin((y-0)/rtube)) && x<=rtube) || (x>=0-rtube && x<=0-rtube*cos(asin((y-0)/rtube))))
						{
							if(y>=0-rhole && y<=rhole)
							{
								if(sqrt(y*y+(z-zcent)*(z-zcent))<=rhole)
								{
									khole++;
									hbounce[i]=0;
									magnumber++;
									do{
										check=moveholebx(hbounce[i],xtube,ytube,ztube,rtube,ltube,znoz,vx,vy,vz,x,y,z,dterror,dtrerror,rhole,lhole,magnet,zcent);
										magnet=0;
										mag=true;
										hbounce[i]=hbounce[i]+1;
										if(check==10)
										{
											cout<<"MAJOR ERROR 2"<<endl<<endl<<endl<<endl<<endl<<endl;
											flag=0-7;
											break;
										}
										
										if(check==-3)
										{
											flag=0-5;
											kpump++;
											kfail++;
//											cout<<"pumped"<<endl;
											break;
										}
										
										if(check==-2)
										{
											flag=4;
											flag=movetube(bounce[i],xtube,ytube,ztube,rtube,ltube,znoz,vx,vy,vz,x,y,z,dterror,dtrerror,sbounce,magnumber);
//											cout<<"Exiting pump hole:	vx: "<<vx<<"	vy: "<<vy<<"	vz: "<<vz<<endl;
//											cout<<"Exiting pump hole"<<endl;
//											cout<<"Exiting pump hole:	x: "<<x<<"	y: "<<y<<"	z: "<<z-zcent<<"	vx: "<<vx<<"	vy: "<<vy<<"	vz: "<<vz<<endl<<"	R: "<<sqrt(x*x+y*y)<<"	r: "<<sqrt(y*y+(z-zcent)*(z-zcent))<<"	bounce: "<<bounce[i]<<endl;
										}
										
										if(check==1)
										{
											avthetad=newthetabmx(vx,vy,vz,x,y);
											bounce[i]=bounce[i]+1;
											magnet=1;
//											mag=true;
//											cout<<"Hit bottom of magnet"<<endl;
										}
										
										if(check>1)
										{
											errorout<<endl<<"Hole Error: "<<dterror+dtrerror<<":"<<endl;
											if(check==2)
											{
												errorout<<"Infinite Velocity: "<<dterror<<":"<<endl<<i<<" "<<x<<" "<<y<<" "<<z-zcent<<" "<<vx<<" "<<vy<<" "<<vz<<" "<<bounce[i]<<" "<<hbounce[i];	
											}
				
											if(check==3)
											{
												errorout<<"Negative Time: "<<dtrerror<<":"<<endl<<i<<" "<<x<<" "<<y<<" "<<z-zcent<<" "<<vx<<" "<<vy<<" "<<vz<<" "<<bounce[i]<<" "<<hbounce[i];
											}
											
											flag=0-6;
											break;
										}
										
										if(check==0)
										{
											avthetad=newthetabhx(vx,vy,vz,x,y,z,zcent);
											bounce[i]=bounce[i]+1;
//											cout<<"Hit side of hole wall:	x: "<<x<<"	y: "<<y<<"	z: "<<z-zcent<<"	vx: "<<vx<<"	vy: "<<vy<<"	vz: "<<vz<<endl<<"	R: "<<sqrt(x*x+y*y)<<"	r: "<<sqrt(y*y+(z-zcent)*(z-zcent))<<"	bounce: "<<bounce[i]<<endl;
//											cout<<"Hit side of hole wall"<<endl;
										}
										
									}while(check>=0);
//									cout<<"broken"<<endl<<endl;
									break;
									
								}
								
							}
							
						}
						
					}
					
				}
				
			}
			
			if(flag==0)
			{
				avthetad=newtheta(vx,vy,vz,x,y);
				bounce[i]=bounce[i]+1;
			}
			
		}while(flag>=0);
		
		if(kpass>spass)
		{bout<<endl<<"Pass "<<1<<" "<<bounce[i];}
		
		if(kfail>sfail)
		{bout<<endl<<"Fail "<<0<<" "<<bounce[i];}
		
		spass=kpass;
		sfail=kfail;
		
	}
	
	
	
	cout<<"Length/radius: "<<ltube/rtube<<endl;
	cout<<"Percent through channel: "<<double(kpass)/double(kpass+kfail)*100.0<<endl;
	cout<<"Percent Rejected: "<<double(kfail-kpump)/double(kpass+kfail)*100.0<<endl;
	cout<<"Percent Pumped: "<<double(kpump)/double(kpass+kfail)*100<<endl;
//	cout<<"Percent entering hole: "<<double(khole)/double(kpass+kfail)*100<<endl;
//	cout<<"Percent entering hole that are pumped: "<<double(kpump)/double(khole)*100<<endl;
	cout<<kpass+kfail<<" "<<ntot<<endl;
	
	
	
}

