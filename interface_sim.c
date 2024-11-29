//This code is designed to simulate particles passing through a single tube, in order to extract their exit distribution.
//The original use for this code was to help understand helium-3 atoms at 1 Kelvin passing through a microchannel plate.
//Some of the parameters will be tailored to this cause and may need to be changed if this is used for anything else.

//Libraries
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
//#define ntot 1000000		//Number of particles
#define ntot 10000000
//#define ntot 50
//#define ntot 1000
//#define ntot 100000
#define stept 3e-6		//Time steps, in seconds
#define gradq 112.5		//Gradient in quads, T/m (9 kG on the tip)
//#define gradq 56.25
//#define gradq 99.2
//#define gradq 50
#define a1 20
#define a2 0.65
//Definitions are available to all functions. However, they generally cannot/should not be altered within the code. This makes them a good substitute for global variables, but only for variables that are expected to remain constant throughout the code.


//Since some of these functions call each other (in their definitions). It is a good idea to declare each function before defining them by writing out all of the headers. A function can be called as long as it has been declared. If it hasn't been defined yet, the code will search for the definition. This is the easiest way to ensure that no errors occur due to a function being called before it is declared. MB
double getnew(double rnoz, double xnoz, double ynoz, double znoz, double& x, double& y, double& z, double& vx, double& vy, double& vz, double& xmem, double& ymem, double& zmem);

double ranr(double r0);

double randgen(double min, double max);

double ransincos(double max);

double ranmax();

int movetube(int bounce, double xtube, double ytube, double ztube, double rtube, double ltube, double znoz, double vx, double vy, double vz, double& x, double& y, double& z, int& dterror, int& dtrerror);

double newtheta(double& vx, double& vy, double& vz, double& x, double & y);
//int moveaper(double& x, double& y, double& z, double vx, double vy, double vz, double xaper, double yaper, double zaper, double raper, double zskim);

//int move4pol(double& x, double& y, double& z, double& vx, double& vy, double &vz, double xtube, double ytube, double ztube, double rtube, double ltube, double zaper, double npol, std::ofstream& out2);

int movefree(double& x, double& y, double& z, double vx, double vy, double vz, double zRGA, double ztube, double ltube);
//Inside of the parenthesis for each declaration are the parameters that each function expects as an input. This is how I chose to pass needed variables from function to function. I would have prefered global variables for some of these, but I am using root and its most recent versions do not treat global variables like c++ does. This method ensures that it will work properly in both c++ and root. If a parameter is initialized with an ampersand after the type, this passes the actual variable to the function and allows it to modify the variable before passing it back. If there is no ampersand, the variable is simply copied to the function; any modifications made to those variables will be lost once the function ends.
//All of these headers will reappear right before the definition of the function itself. MB


double getnew(double rnoz, double xnoz, double ynoz, double znoz, double& x, double& y, double& z, double& vx, double& vy, double& vz, double& xmem, double& ymem, double& zmem)	/* generates a new particle at(0,0,0) with random direction */
//Actually on the z coordinate is guaranteed to be zero. x and y are randomized to a certain point on the cross-sectional area of the nozzle. MB
{
	double theta,phi,thetamax;
	double v,r;
//	Values needed by this function


//	cout<<"one"<<endl;
//	r=randgen(0.,rnoz);
	r=ranr(rnoz);		//Calls the ranr() function (defined in code) and sets "r" equal to the return value
//	cout<<"two"<<" "<<r<<endl;
//	cout<<r<<endl;
	phi=randgen(0.,2*pi);	//Calls the randgen() function (defined in code) and sets phi to the return value. MB
//	cout<<"three"<<" "<<phi<<endl;
	x=r*sin(phi)+xnoz;	//Sets starting value for x at some random position on the nozzle. Just a conversion from spherical to cartesian coordinates and then adding the x coordinate location of the nozzle. MB
	y=r*cos(phi)+ynoz;	//Same thing as x, but for y. MB
	z=znoz;			//"z" is the direction which the nozzle is facing. All particles start at the same location in z. MB
//	thetamax=atan((raper+rnoz)/(zaper-znoz));
//	Sets some value for thetamax base on the geometery. MB
	thetamax=pi/2;	//Thetamax is pi/2.
//	thetamax=0.3;		//Geometery doesn't know what the hell it's talking about. This value is, apparently, better. MB
//	thetamax=0.01;		//Test value that I'm using, 100% into quad. MB
//	thetamax=0.001;
	theta=ransincos(thetamax);
//	theta=randgen(0.,thetamax);
//	theta=ransin(thetamax);
//	Time to find a random value for theta that is less than thetamax. Uses the ransincos() function defined in the code. MB
//	cout<<"four"<<" "<<theta<<endl;
	v=ranmax()*1000;		//sets v to the return value of ranmax(). MB
//	cout<<v<<endl;
//	cout<<"five"<<" "<<v<<endl;
	phi=randgen(0.,2*pi);	//Get new value for phi, for velocity. I guess that phi isn't needed later. MB
//	cout<<"six"<<" "<<phi<<endl;
	vx=v*sin(theta)*cos(phi);
	vy=v*sin(theta)*sin(phi);
	vz=v*cos(theta);
//	Break up v into its constituent cartesian coordinate parts.
	xmem=vx;
	ymem=vy;
	zmem=sqrt(vx*vx+vy*vy);
//	Store initial velocities in x, y, and z-mem. Not sure why zmem isn't equal to vz. Will have to look over code more. MB
	return theta;
}


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
//	Variable just for this function
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
//	cout<<RAND_MAX<<endl;
//	fvmax=5*vmax;
	do
	{
		v=randgen(0.,5*vmax);
		pcur=v*v*v*exp(-v*v/(v1*v1)-v2*v2/(v*v))/pmax;
		comp=randgen(0.,1.);
//		cout<<v<<" "<<pcur<<" "<<comp<<endl;
	}while(comp>pcur);
//	Again, I need to double check this distribution. MB

	return v;	
}


int movetube(int bounce, double xtube, double ytube, double ztube, double rtube, double ltube, double znoz, double vx, double vy, double vz, double& x, double& y, double& z, int& dterror, int& dtrerror)		/* moving the particle through skimmer */
{
	double dtr,dtr1,dtr2,dtz,xdum,ydum,zdum,rdum,vr;	//Some function parameters. MB
	xdum=x;					//X position before tube. MB
	ydum=y;					//Y position before tube. MB
	zdum=z;					//Z position before tube. MB
	rdum=sqrt(x*x+y*y);
	
	double xl1, xl2, L1, L2, phi, phir, dot, theta;
	
	vr=sqrt(vx*vx+vy*vy);
	phi=atan(vy/vx);
	phir=atan(y/x);
	theta=atan(vy/vx);
	dot=(x*vx+y*vy)/(vr*rdum);
	dtz=0;
	dtr=0;
	

	if(bounce>=0)
	{
		if(vr!=0)
		{
/*
			xl1=(x*sin(phi)-y*cos(phi))*sin(phi)+cos(phi)*sqrt((rtube*rtube)-pow(y*cos(phi)-x*sin(phi),2));
			xl2=(x*sin(phi)-y*cos(phi))*sin(phi)-cos(phi)*sqrt((rtube*rtube)-pow(y*cos(phi)-x*sin(phi),2));
	
			L1=xl1-x;
			L2=xl2-x;
	
			dtr1=L1/vx;
			dtr2=L2/vx;
		
			if(dtr1<0 && dtr2<0)
			{
				dtrerror++;
				cout<<rdum<<" "<<dot<<" "<<dtr1<<" "<<dtr2<<endl;
				return 3;
			}
*/
//*
			dtr1=(sqrt(vr*vr*((rtube*rtube)-(rdum*rdum))+(x*vx+y*vy)*(x*vx+y*vy))-x*vx-y*vy)/(vr*vr);
			dtr2=(0-sqrt(vr*vr*((rtube*rtube)-(rdum*rdum))+pow(x*vx+y*vy,2))-x*vx-y*vy)/(vr*vr);
//			dtr2=(0-sqrt(vr*vr*rtube*rtube-pow(x*vy-y*vx,2))-x*vx-y*vy)/(vr*vr);

//			cout<<dtr1<<" "<<dtr2<<endl;
		
			if(dtr1<0 && dtr2<0)
			{
				dtrerror++;
				cout<<"bounce"<<" "<<rdum<<"	"<<x*vx<<"	"<<y*vy<<"	"<<dtr1<<"	"<<dtr2<<endl;
				return 3;
			}
//*/		
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
//*/		
		}
		
	}

/*
	if(bounce>0)
	{
		
		if(vr!=0)
		{
			dtr=(2*rtube/vr)*cos(theta);
		}
		
	}
*/
	if(vz>0)
	{
		dtz=(ztube+ltube-zdum)/vz;	//Time it takes particle to traverse the tube/channel (assuming it doesn't hit an obstacle). MB
	}
	
	if(vz<0)
	{
		dtz=(ztube-zdum)/vz;
	}
	
	
	if(dtr<dtz)
	{
		x=x+vx*dtr;				//X position after movetube.
		y=y+vy*dtr;				//Y position after movetube.
		z=z+vz*dtr;				//Z position after movetube.
		return 0;
	}
	
	else if(dtz<=dtr)
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


int movefree(double& x, double& y, double& z, double vx, double vy, double vz, double zRGA, double ztube, double ltube)
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

double newtheta(double& vx, double& vy, double& vz, double& x, double& y)
{
	double theta,phi,phiprime,thetamax;
	double v=sqrt(vx*vx+vy*vy+vz*vz);
	
	phi=atan2(y,x);
/*	
	double phiprimemin=phi+pi/2;
	double phiprimemax=phi+3*pi/2;
*/	
	thetamax=pi/2;	//Thetamax is pi/2.
	theta=ransincos(thetamax);
	phiprime=randgen(0.,2*pi);
//	phiprime=randgen(phiprimemin,phiprimemax);

	vx=v*(cos(phi)*cos(pi-theta)-sin(phi)*cos(phiprime)*sin(pi-theta));
	vy=v*(sin(phi)*cos(pi-theta)+cos(phi)*cos(phiprime)*sin(pi-theta));
	vz=v*sin(phiprime)*sin(pi-theta);
	
//	cout<<theta*180/pi<<"	"<<phiprime*180/pi<<"	"<<phi*180/pi<<"	"<<atan(vy/vx)*180/pi<<endl;
	
	return theta;
}


void interface_sim()
{
	double x, y, z, vx, vy, vz;
//	particle position and velocity in cartesian coordinates
	double xnoz, ynoz, znoz, rnoz, xtube, ytube, ztube, rtube, ltube;
//	nozzle and tube positions (in cartesian coordinates) and radii.
	double zRGA;
//	Inclusion by me for the movefree function. MB
	int polcount=0;
//	Iterator for polarization array. MB
	double xmem, ymem, zmem;
//	x and y velocity of particle. zmem is the polar speed of the particle. Three are saved here to be compared with end velocities.
	int ktot, kpass, kfail;
	int dterror=0;
	int dtrerror=0;
	int flag;
	
	int *bounce= new int[ntot];

	Double_t *xstore= new Double_t[ntot];
	Double_t *ystore= new Double_t[ntot];
	Double_t *zstore= new Double_t[ntot];
	Double_t *vxstore= new Double_t[ntot];
	Double_t *vystore= new Double_t[ntot];
	Double_t *vzstore= new Double_t[ntot];
	Double_t *exittheta= new Double_t[ntot];
	Double_t *exitphi= new Double_t[ntot];
	
	double exitvrho;
//	I want to store the position of all the particles that make it through the quad. This way I can reconstruct the predicted profile. MB
	
	
	xnoz=0.0;   ynoz=0.0;   znoz=0.0;	rnoz=0.005;						//Initial position for channel.
	xtube=0.0;  ytube=0.0;  ztube=0.0;  rtube=rnoz;		ltube=0.5;		//Single channel stuff.

	zRGA=ztube+ltube+40;

	ltube=50*rnoz;
/*	
	rnoz=5;
	
	rnoz=0.0508/2; rtube=rnoz; ltube=5.08;		//Ghost Ranch Dimensions
//*/	
//	rtube=rnoz;

	rnoz=6.337/2;	rtube=6.3;	ltube=1295/2;	ltube=1295*0.1;
	
	ofstream output;
	output.open("sc_profile.txt");
	//Create variable to hold output file and create/open/rewrite an output file.

	ofstream putty;
	putty.open("sc_initial_profile.txt");
	
	ofstream errorout;
	errorout.open("SC_error_file.txt");
	
	ofstream track;
	track.open("sc_tracker.txt");
	
	
	srand((unsigned)time(NULL));	//Uses internal clock to fetch seed for random number generator.
	rand();
	
	double avthetas=0;
	double avthetaf=0;
	double avthetad=0;
	int storenum;
	
	
	
	ktot=0; kpass=0; kfail=0;
		
	for(int i=0; i<ntot; i++)
	{
		avthetas+=getnew(rnoz,xnoz,ynoz,znoz,x,y,z,vx,vy,vz,xmem,ymem,zmem);
		bounce[i]=0;
		flag=0;
		
		storenum=i;
		
//		cout<<i<<" "<<kpass+kfail<<endl;
		
//		track<<endl<<i<<" "<<x<<" "<<y<<" "<<z;
		
		do
		{
//			cout<<sqrt(x*x+y*y+z*z)<<" "<<vx<<" "<<vy<<" "<<vz<<endl;
			
			flag=movetube(bounce[i],xtube,ytube,ztube,rtube,ltube,znoz,vx,vy,vz,x,y,z,dterror,dtrerror);
			
//			track<<endl<<i<<" "<<x<<" "<<y<<" "<<z;
			if(flag==10)
			{
				cout<<"MAJOR ERROR"<<endl;
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
				storenum++;
			}
			
			if(z<=ztube+0.0000000000001 && bounce[i]>0)
			{
				flag=0-2;
				kfail++;
				storenum++;
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
				}
			
				break;
			}
			
			if(sqrt(x*x+y*y)>rtube+0.000001*rtube)
			{
				dterror++;
				errorout<<endl<<"Error: "<<dterror+dtrerror<<":"<<endl;
				errorout<<"Exceeded tube radius: "<<dterror<<":"<<endl<<i<<" "<<x<<" "<<y<<" "<<z<<" "<<vx<<" "<<vy<<" "<<vz<<" "<<bounce[i]<<" "<<flag;
				cout<<sqrt(x*x+y*y)<<"	"<<rtube<<endl;
				break;
			}
			
			if(flag==0)
			{
				avthetad=newtheta(vx,vy,vz,x,y);
				bounce[i]=bounce[i]+1;
			}
			
		}while(flag>=0);
		
		if(storenum<i+1)
		{
			cout<<kpass+kfail<<" "<<i+1<<":"<<endl<<"flag: "<<flag<<"	x: "<<x<<"	y: "<<y<<"	z: "<<z<<"	vx: "<<vx<<"	vy: "<<vy<<"	vz: "<<vz<<endl<<"	R: "<<sqrt(x*x+y*y+z*z)<<endl;
		}
		
/*		
		if(flag==(0-1))
		{
			movefree(x,y,z,vx,vy,vz,zRGA,ztube,ltube);
			output<<endl<<i<<" "<<x<<" "<<y<<" "<<z<<" "<<vx<<" "<<vy<<" "<<vz<<" "<<exittheta[kpass]<<" "<<exitphi[kpass]<<" "<<bounce[i];
		}
*/			
	}
	
	cout<<"Length/radius: "<<ltube/rtube<<endl;
	cout<<"Percent through channel: "<<double(kpass)/double(kpass+kfail)*100.0<<endl;
	cout<<"Percent Rejected: "<<double(kfail)/double(kpass+kfail)*100.0<<endl;
	cout<<kpass<<" "<<kfail<<" "<<kpass+kfail<<endl;
	
	output.close();
	putty.close();
	errorout.close();
	track.close();
	
	
	
	
	
	
	
}




















