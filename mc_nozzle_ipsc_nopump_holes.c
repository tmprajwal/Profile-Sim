//This code is designed to generate helium-3 atoms with a randomly assigned polarizations and veocities. Velocity assignments can be weighted in different ways in order to more accurately simulate a gas being expelled from the nozzle.

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
//#define ntot 1000000000		//Number of particles
//#define ntot 1000000000
//#define ntot 50
//#define ntot 1
#define ntot 100000
//#define stept 3e-6		//Time steps, in seconds
#define stept 3e-6
#define gradq 112.5		//Gradient in quads, T/m (9 kG on the tip)
//#define gradq 56.25
//#define gradq 99.2
//#define gradq 50
#define a1 20
#define a2 0.65
//Definitions are available to all functions. However, they generally cannot/should not be altered within the code. This makes them a good substitute for global variables, but only for variables that are expected to remain constant throughout the code.


//Since some of these functions call each other (in there definitions). It is a good idea to declare each function before defining them by writing out all of the headers. A function can be called as long as it has been declared. If it hasn't been defined yet, the code will search for the definition. This is the easiest way to ensure that no errors occur due to a function being called before it is declared. MB
double getnew(int& itest, double rnoz, double xnoz, double ynoz, double znoz, double& x, double& y, double& z, double& vx, double& vy, double& vz, int& npol, double& xmem, double& ymem, double& zmem, double raper, double zaper);

double ranr(double r0);

double randgen(double min, double max);

double ransincos(double max);

double ransin(double max);

double ranmax();

int ranpol();

int moveskim(double xskim, double yskim, double zskim, double rskim, double znoz, double vx, double vy, double vz, double& x, double& y, double& z);

int moveaper(double& x, double& y, double& z, double vx, double vy, double vz, double xaper, double yaper, double zaper, double raper, double zskim);

int move4pol(double& x, double& y, double& z, double& vx, double& vy, double &vz, double xtube, double ytube, double ztube, double rtube, double ltube, double zaper, double npol, std::ofstream& out2);

int movetube(double& x, double& y, double& z, double& vx, double& vy, double& vz, double xtube, double ytube, double ztube, double rtube, double ltube, double zaper, double npol, std::ofstream& out2);

int movefree(double& x, double& y, double& z, double vx, double vy, double vz, double zRGA, double ztube, double ltube);

double fitg(double *xf, double *par);

double fitgc(double *xf, double *par);

double ranmc(double max);
//Inside of the parenthesis for each declaration are the parameters that each function expects as an input. This is how I chose to pass needed variables from function to function. I would have prefered global variables for some of these, but I am using root and its most recent versions do not treat global variables like c++ does. This method ensures that it will work properly in both c++ and root. If a parameter is initialized with an ampersand after the type, this passes the actual variable to the function and allows it to modify the variable before passing it back. If there is no ampersand, the variable is simply copied to the function; any modifications made to those variables will be lost once the function ends.
//All of these headers will reappear right before the definition of the function itself. MB

double ransin(double max)     /* generates random number between 0 and max with sin distribution */
{
double theta,comp;
    do
	{
	theta=randgen(0.,max);
	comp=randgen(0.,1.);
	} while(comp>sin(theta));
	return theta;
}

double getnew(int& itest, double rnoz, double xnoz, double ynoz, double znoz, double& x, double& y, double& z, double& vx, double& vy, double& vz, int& npol, double& xmem, double& ymem, double& zmem, double raper, double zaper)	/* generates a new particle at(0,0,0) with random direction */
//Actually on the z coordinate is guaranteed to be zero. x and y are randomized to a certain point on the cross-sectional area of the nozzle. MB
{
	double theta,phi,thetamax;
	double v,r;
//	Values needed by this function


//	cout<<"one"<<endl;
//	r=randgen(0.,rnoz);
	r=ranr(rnoz);		//Calls the ranr() function (defined in code) and sets "r" equal to the return value
//	cout<<"two"<<" "<<r<<endl;
	phi=randgen(0.,2*pi);	//Calls the randgen() function (defined in code) and sets phi to the return value. MB
//	cout<<"three"<<" "<<phi<<endl;
	x=r*sin(phi)+xnoz;	//Sets starting value for x at some random position on the nozzle. Just a conversion from spherical to cartesian coordinates and then adding the x coordinate location of the nozzle. MB
	y=r*cos(phi)+ynoz;	//Same thing as x, but for y. MB
	z=znoz;			//"z" is the direction which the nozzle is facing. All particles start at the same location in z. MB
/*
//	thetamax=atan((raper+rnoz)/(zaper-znoz));
//	Sets some value for thetamax base on the geometery. MB
	thetamax=pi/2;	//Thetamax is pi/2.
//	thetamax=0.3;		//Geometery doesn't know what the hell it's talking about. This value is, apparently, better. MB
//	thetamax=0.01;		//Test value that I'm using, 100% into quad. MB
//	thetamax=0.001;
//	theta=ransincos(thetamax);
	theta=ranmc(thetamax);
//	theta=randgen(0.,thetamax);
//	theta=ransin(thetamax);
//	Time to find a random value for theta that is less than thetamax. Uses the ransincos() function defined in the code. MB
//	cout<<"four"<<" "<<theta<<endl;
	v=ranmax();		//sets v to the return value of ranmax(). MB
//	cout<<"five"<<" "<<v<<endl;
	phi=randgen(0.,2*pi);	//Get new value for phi, for velocity. I guess that phi isn't needed later. MB
//	cout<<"six"<<" "<<phi<<endl;
	vx=v*sin(theta)*cos(phi);
	vy=v*sin(theta)*sin(phi);
	vz=v*cos(theta);
//	Break up v into its constituent cartesian coordinate parts.
//*/
	npol=ranpol();		//Uses ranpol() to randomly assign 1 or -1 to npol. MB
//	npol=-1;
//	cout<<"seven"<<" "<<npol<<endl;
/*
	xmem=vx;
	ymem=vy;
	zmem=sqrt(vx*vx+vy*vy);
//	Store initial velocities in x, y, and z-mem. Not sure why zmem isn't equal to vz. Will have to look over code more. MB
	return theta;
*/
	return atan2(sqrt(vx*vx+vy*vy),vz);
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

double ranmc(double max)
{
	double theta, comp;
	do
	{
		theta=randgen(0.0,max);
		comp=randgen(0.0,1.0);
	}while(comp>pow(1+pow(2*theta/a1,2)*(pow(2,1/a2)-1),0-a2));
	return theta;
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


int ranpol()	/* generates randomly 1 or -1  */
{
	double r, pol;
	r=randgen(0.,1.);	//randomly generate value between 0 and 1 using randgen() and set r equal to this value. MB
	if(r>0.5)
	{pol=1;}
	else if(r<=0.5)
	{pol=-1;}
	return pol;
}


int moveskim(double xskim, double yskim, double zskim, double rskim, double znoz, double vx, double vy, double vz, double& x, double& y, double& z)		/* moving the particle through skimmer */
{
	double dt,xdum,ydum;	//Some function parameters. MB
	z=zskim;		//Setting starting z position for particle. MB
	dt=(zskim-znoz)/vz;	//Time it takes particle to traverse the skimmer (assuming it doesn't hit an obstacle). MB
	xdum=x;			//X position before skimmer. MB
	ydum=y;			//Y position before skimmer. MB
	x=x+vx*dt;		//X position after skimmer (assuming it hasn't hit an obstacle). MB
	y=y+vy*dt;		//Y position after skimmer (assuming it hasn't hit an obstacle). MB
	if(sqrt((x-xskim)*(x-xskim)+(y-yskim)*(y-yskim))>=rskim)
	{return 0;}
//	If particle would have hit skimmer, return 0. MB
	else if(sqrt((x-xskim)*(x-xskim)+(y-yskim)*(y-yskim))<rskim)
	{return 1;}
//	If particle did not hit skimmer, return 1. MB
	return 0;
}


int moveaper(double& x, double& y, double& z, double vx, double vy, double vz, double xaper, double yaper, double zaper, double raper, double zskim)		/* moving the particle through aperture */
{
	double dt;
	z=zaper;		//Set particle z position. MB
	dt=(zaper-zskim)/vz;	//Time it takes particle to traverse second skimmer (aperature). MB
	x=x+vx*dt;
	y=y+vy*dt;
	if(sqrt((x-xaper)*(x-xaper)+(y-yaper)*(y-yaper))>=raper)
	{return -1;}
	else if(sqrt((x-xaper)*(x-xaper)+(y-yaper)*(y-yaper))<raper)
	{return 1;}
	return -1;
//	This is all very similar to moveskim(). See moveskim() for more comments. MB
}


int move4pol(double& x, double& y, double& z, double& vx, double& vy, double& vz, double xtube, double ytube, double ztube, double rtube, double ltube, double zaper, double npol, std::ofstream& out2)		/* moving the particle through 6poles */
{
	double dt;
//	Time step. MB
	double r, rlocal;
//	Particle position and offset from tube center. MB
	double accx, accy, accz, delvx, delvy, delvz;
//	Acceleration in x, y, and z. And change in velocity in x, y, and z. MB
	double bx, by, bz, btot;
//	Magnetic field in x, y, z, and the total. MB
	double dbxdx, dbxdy, dbxdz, dbydx, dbydy, dbydz, dbzdx, dbzdy, dbzdz;
//	Magnetic field position based partial differentials. MB
	double dz, zmin;
//	Step in z and ? MB
	int flag1;
//	A new flag for tracking if particle hits anything. MB
	double ystart, zstart;
//*
	double accpx, accpy, Fp, SA, rHe, IP;
	IP=2e-4;
//	IP=0;
	rHe=31e-12;
	SA=2*pi*rHe*rHe;
	Fp=IP*SA;
//*/
	z=ztube;
	dt=(ztube-zaper)/vz;
	x=x+vx*dt;
	y=y+vy*dt;
//	Moving particle into the quad entrace. MB
	zmin=2000;
//	out2<<x<<" "<<y<<" "<<z<<" "<<vx<<" "<<vy<<" "<<vz<<" "<<npol<<endl;
	if(sqrt((x-xtube)*(x-xtube)+(y-ytube)*(y-ytube))>rtube)
	{return -2;}	//Checking to see if the particle hit the entrance of the quad tube/pipe. MB
	do
	{
		x=x+vx*stept*1000;
		y=y+vy*stept*1000;
		z=z+vz*stept*1000;
//		Move each particle a little bit forward. Extra factor of 1000 comes from the fact that vi is in m/s while ri is in mm. When calculating dt in each function it comes out in units of mm*s/m which, when multiplied by vi yields a value in mm. MB
//		out2<<x<<" "<<y<<" "<<z<<" "<<vx<<" "<<vy<<" "<<vz<<" "<<npol<<endl;
		dz=sqrt((z-(ztube+ltube-40))*(z-(ztube+ltube-40)));
		if(dz<zmin)
		{
			ystart=y;
			zstart=x;
			zmin=dz;
/*memorizing x,y at 4 cm before exit. Genya*/
//			I don't understand what's supposed to be happening here yet.
		}

		r=sqrt(x*x+y*y);
		rlocal=sqrt((x-xtube)*(x-xtube)+(y-ytube)*(y-ytube)); /* includes tube offset. Genya */	
		//Particle distance from the center of the quad tube's cross-sectional area. MB
		if(rlocal>=rtube)
		{return -3;}	//Checking to see if the particle hit the tube. MB
		bx=gradq*y/1000;
		by=gradq*x/1000;
		bz=0;
/*
		bx+=randgen(0,bx/10)*ranpol();
		by+=randgen(0,by/10)*ranpol();
		bz=randgen(0,(bx+by)/(2*100))*ranpol();
*/
		btot=sqrt(bx*bx+by*by+bz*bz);
//		Reevaluate Magnetic field values at new positions. MB
		dbxdx=0;//randgen(0,gradq/10)*ranpol();
		dbxdy=gradq;//+randgen(0,gradq/10)*ranpol();
		dbxdz=0;//randgen(0,gradq/10)*ranpol();
		dbydx=gradq;//+randgen(0,gradq/10)*ranpol();
		dbydy=0;//randgen(0,gradq/10)*ranpol();
		dbydz=0;//randgen(0,gradq/10)*ranpol();
		dbzdx=0;//randgen(0,gradq/10)*ranpol();
		dbzdy=0;//randgen(0,gradq/10)*ranpol();
		dbzdz=0;//randgen(0,gradq/10)*ranpol();
//		Set values for the magnetic field partial differentials. MB
		accx=(mu/mkg)*((bx*dbxdx+by*dbydx+bz*dbzdx)/btot)*npol;
		accy=(mu/mkg)*((bx*dbxdy+by*dbydy+bz*dbzdy)/btot)*npol;
		accz=(mu/mkg)*((bx*dbxdz+by*dbydz+bz*dbzdz)/btot)*npol;
//		Reevaluate acceleration due to magnetic field and gradients. MB
/*
		if(x!=0 || y!=0)
		{
			if(r>=0.3)
			{
				accpx=(Fp/(mkg*r*r*r))*x;
				accpy=(Fp/(mkg*r*r*r))*y;
			}

			if(r<0.3)
			{
				accpx=(Fp/(mkg*r))*x;
				accpy=(Fp/(mkg*r))*y;
			}

		}

		if(x==0 && y==0)
		{
			accpx=0;
			accpy=0;
		}

		accx+=accpx;
		accy+=accpy;
//*/
//		cout<<npol<<endl;
		delvx=accx*stept;
		delvy=accy*stept;
		delvz=accz*stept;
//		Reevaluate the change in particle velocity due to magnetic field. MB
		vx=vx+delvx;
		vy=vy+delvy;
		vz=vz+delvz;
//		Change particle velocity accordingly. MB
	}while(z<ztube+ltube);
//	Continue this loop until particle traverse the full length of the quad. MB

	return 1;	//If particle doesn't hit or exit tube through pumping holes, return 1. MB	
}
int movetube(double& x, double& y, double& z, double& vx, double& vy, double& vz, double xtube, double ytube, double ztube, double rtube, double ltube, double zaper, double npol, std::ofstream& out2)			/* moving the particle through skimmer */


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
//		dtr2=(0-sqrt(vr*vr*rtube*rtube-pow(x*vy-y*vx,2))-x*vx-y*vy)/(vr*vr);

//		cout<<dtr1<<" "<<dtr2<<endl;
		
		if(dtr1<0 && dtr2<0)
		{
//			dtrerror++;
			cout<<" "<<rdum<<"	"<<x*vx<<"	"<<y*vy<<"	"<<dtr1<<"	"<<dtr2<<endl;
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
//		dterror++;
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
//	if(sqrt((x-xskim)*(x-xskim)+(y-yskim)*(y-yskim))>=rskim)
//	{return 0;}
//	If particle would have hit skimmer, return 0. MB
//	else if(sqrt((x-xskim)*(x-xskim)+(y-yskim)*(y-yskim))<rskim)
//	{return 1;}
//	If particle did not hit skimmer, return 1. MB
}



//Start of program main
void mc_nozzle_ipsc_nopump_holes()
{
	cout<<"hello"<<endl;
	
//	Parameters that need to be available for multiple different functions
//	double pi, c, kbol, mu, mhe3, mkg;
//	Pi, speed of light, boltzman constant, mass of helium-3 atom, that mass in kg, magnetic moment of helium-3
//	double tnoz;
//	nozzle temperature
	double x, y, z, vx, vy, vz;
//	particle position and velocity in cartesian coordinates
//	double x4, y4, z4, vx4, vy4, vz4;
//	for saving particle position and velocity at end of quad
	double xnoz, ynoz, znoz, rnoz, xskim, yskim, zskim, rskim;
//	nozzle and skimmer positions (in cartesian coordinates) and radii. 
	double xaper, yaper, zaper, raper, xtube, ytube, ztube, rtube, ltube, xdiaph, ydiaph, zdiaph, rdiaph;
//	positions and radii for the quadrupole aperture, tubing, and baffle diaphram
	double xpump, ypump, zpump, lpump, rpump, xpipe4, ypipe4, zpipe4, lpipe4, rpipe4, xcell, ycell, zcell, rcell, lcell;
//	Dimensions for main experiment components that come after the ABS. These will likely change.
	double zRGA;
//	Inclusion by me for the movefree function. MB
//	double stept;
	double stepz;
//	Step sizes for moving the particle in time and the z direction
//	double gradq;
//	Magnetic field gradient in quad
//	int ntot;
//	total number of particles.
	double hist[10000], xhist[305];
//	hist appears to save the speed for each individual particle. I cannot determine what xhist is supposed to be. It depends on multiple variables that are all initialized with no values. This would make all elements of xhist identically zero. Nevermind found where on of its independent variables is defined, in a commented out section of code. I must have ignored it before because it was commented out.
	int yhist[305];
//	yhist appears to be bins that start at zero and increase by one for every value in hist that is inside their bin. khist is used to keep track of bin width.
	int npol;
//	polarization of particle (either 1 or -1). Changes for each particle, not saved for each.
	Int_t *polarization= new Int_t[ntot];
//	int polarization[ntot];
//	Array for saving npol of each particle. MB
	int polcount=0;
//	Iterator for polarization array. MB
	double polmem[10000], xquad;
//	polmem is supposed to save the polarization state of each particle, but in current code, it appears to be defined as identically 1 for all values. xquad just shows how far along the quad the particle is. It shifts the zero point for the particle position to the end of the quad. Not yet used in this code. MB
	double xmem, ymem, zmem;
//	x and y velocity of particle. zmem is the polar speed of the particle. Three are saved here to be compared with end velocities.
	double ystart, zstart;
//	Starting positions for y and (confusingly) x.


	int flag;	//Used to track if an atom hits a component wall. Info used to abort tracking on that particle if it hits a wall. Last location of particle saved later.
	int ktot=0;	//Number of particles tested
	int kskim=0;	//Number of particles that hit skimmer
	int kaper=0;	//Number of particles that hit aperature (second skimmer)
	int kentry=0;	//Number of particles lost on entry to quad
	int kquad=0;	//Number of particles that make it past quad (both bad and good)
	int khole=0;	//Number of particles get pumped through the holes in the quad
	int kback=0;	//Number of particles pumped through entrance into quad (from Genya). However, its flag option can never arise in Genya's code. Not certain what this actually means
	int kleak=0;	//Number of particles that leak through the exit of quad (Genya). Again this flag option can never arise in code as written.
	int kgood=0;	//Number of particles that make it through the quad properly. With no collisions and are properly polarized.
	int kreverse=0;	//Number of particles that reflect off of quad tube and then exit the quad tube through the end near the nozzle.


	int itest=0;	//Iterator for Getnew() function. However, itest currently seems pointless in the code.
	int ihist=0;	//Iterator for hist. Hist collects abs(velocity) of each particle.
	int ipol=0;	//Iterator for spin function. Spin function appears to track spin of each particle. In current state, the function always returns 0 though. Some sections are commented out. this could be for a special case or debugging purpose. However, it will take some time to properly analyze.
	int i=0;	//Likely an iterator. Have to wait until I come across it naturally though. Can't exactly search it in a useful way.
	int flag1;	//Equals the returned value of the spin function. This is then tested to see if its larger than 0. Seems like a pointless extra step. Code could be more compact.
	double sum1;	//I don't even know what he was going for here. sum1 will never change from zero, despite code that clearly seems to expect that it will.
	double dumdum;	//Originally "dum" (Genya). However, a double array in the "spin()" function was also named dum. Confusing. dumdum appears to be the standard error (not deviation) of the average polarization.
	double polsum, polavg, polerr, del;
//	polsum is del squared, polavg is the average polarization, polerr appears to be the standard deviation of this average, and del is the individual deviation from the mean for a single particle. del is overwritten for every new particle, but not until it is included in dumdum.

	int scl=2200792;
	scl=250832;
/*
	This is all stuff that I'm writing. MB
*/
	Double_t *xstore= new Double_t[ntot];
	Double_t *ystore= new Double_t[ntot];
	Double_t *zstore= new Double_t[ntot];
	
	Double_t *vxsc= new Double_t[5000000];
	Double_t *vysc= new Double_t[5000000];
	Double_t *vzsc= new Double_t[5000000];
	Double_t *thetasc= new Double_t[5000000];
	Double_t *phisc= new Double_t[5000000];
	Double_t *xsc= new Double_t[5000000];
	Double_t *ysc= new Double_t[5000000];
	Double_t *zsc= new Double_t[5000000];
	Double_t *number= new Double_t[5000000];
	Double_t *bounce= new Double_t[5000000];

//	double xsc, ysc, zsc, number, bounce;
	
//	double xstore[ntot],ystore[ntot],zstore[ntot];
//	I want to store the position of all the particles that make it through the quad. This way I can reconstruct the predicted profile. MB
/*
	End of my section. MB
*/




/*	pi=3.14159;		//Pi, you know from math and circles and such
	c=3.0e8;		//Speed of light (rounded) in m/sec
	kbol=8.62e-11;        	//Boltzman consant, MeV/K
	mhe3=2809.4;            //He3 mass, MeV
	mkg=mhe3*1.783e-30;     //He3 mass, kg
	mu=1.07e-26;          	//He3 magnetic moment, J/T
	tnoz=1.3;               //Nozzle temperature, K
	ntot=1000000;           //Number of particles
	stept=3e-5;           	//Time steps, in seconds
	gradq=112.5;           	//Gradient in quads, T/m (9 kG on the tip)
*/
//	Geometry, in mm. Filling in some values from start of the main.
	xnoz=0.0;   ynoz=0.0;   znoz=0.0;       rnoz=.4;	rnoz=4.4;	//This is the real rnoz. Don't know what Genya had.
											rnoz=6.337/2;			//MC plate radius.
	xskim=0.0;  yskim=0.0;  zskim=22.6;     rskim=4.75;
	xaper=0.0;  yaper=0.0;  zaper=110.6;    raper=4.75;
	xtube=0.0;  ytube=0.0;  ztube=186.8;    rtube=6.3;	ltube=1295.;
  
	xdiaph=0.0; ydiaph=0.0; zdiaph=ztube+ltube+187.4;   rdiaph=11.1;

	xpump=0.0;  ypump=0.0;  zpump=zdiaph+529.2;         rpump=15.9;   lpump=164.9;
	xpipe4=0.0; ypipe4=0.0; zpipe4=zpump+lpump+133.5;   rpipe4=18.55; lpipe4=199.2;
	xcell=0.0;  ycell=0.0;  zcell=zpipe4+lpipe4+150.1;  rcell=22.25;  lcell=216.7;	//compression tube

	zRGA=ztube+ltube+50.8;
//	My addition for the movefree function. Number indicates how far RGA is from end of quad. Currently just an estimate. MB
//	rtube=10.08;		//Test maximum radius per field gradient.

//	First attempt at changing skimmer and aperture dimensions/locations.
/*
	rnoz=1;
	zskim=57.45;	rskim=1.9;
	zaper=141.8;	raper=3.6;
//*/
/*
	zskim=57.45;	rskim=3.4;
	zaper=141.8;	raper=4.3;
//*/
/*
//	For dimensions found in solidworks:
	zskim=57;	rskim=3.4;
	zaper=129;	raper=3.8;		//raper=3.75 yields QPass=0.246248% and Lost Q Face=0.007142%
//*/
/*
	zskim=92.45;	rskim=???;		//Best dimensions for poor mcplate fit; rskim>3.5; trying 4.0.
	zaper=176.8;	raper=5.793;
//*/	
//*
	ztube=214.6;
	zskim=57;	rskim=3.4;
	zaper=170.23; raper=4.5;
//*/
//	Adjusting raper:
//	raper=5.793;	// Does not appear to have any negative effect on fraction exiting quad.
//	Probably only blocks particles that would have hit the quad tube on the inside anyway with only a negligible amount of exceptions.

//	Adjusting rskim:
//	rskim=4.0;

//	Ghost Ranch nozzle dimensions
//	rnoz=0.45; zskim=40;

/*
	double zTAmin=ztube-zaper;
//	zskim=35;
//	zaper=130;
//	ztube=zaper+zTAmin;
//*/
/*
	FILE *in;
	FILE *out;
//	File IO! Genya's version
*/

	ofstream output;
	output.open("profile.txt");
	//Create variable to hold output file and create/open/rewrite an output file.

	ofstream putty;
	putty.open("initial_profile.txt");

	ofstream out2;
	out2.open("position.txt");
	
	ofstream steve;
	steve.open("exit_profile.txt");

	srand((unsigned)time(NULL));	//Uses internal clock to fetch seed for random number generator.
	rand();				//Genya had this here it calls the rand function but does not store the value. Perhaps the first value in the sequence wasn't random enough?

	ifstream input;
	int num=0;
	int inpi=0;
	int inpj=0;
	int inpk=0;
	//Create variable to hold input file and useful input counters. Me
	Double_t ***xcart= new Double_t**[1056];
	Double_t ***ycart= new Double_t**[1056];
	Double_t ***zcart= new Double_t**[1056];
//	double xcart[1056][21][21], ycart[1056][21][21], zcart[1056][21][21];
	Double_t *Hx= new Double_t[500000];
	Double_t *Hy= new Double_t[500000];
	Double_t *Hz= new Double_t[500000];
	Double_t *xf= new Double_t[500000];
	Double_t *yf= new Double_t[500000];
	Double_t *zf= new Double_t[500000];
//	double Hx[500000], Hy[500000], Hz[500000];
//	double xf[500000], yf[500000], zf[500000];
	//Create some useful variables for input data. Me
	Double_t ***hx= new Double_t**[1056];
	Double_t ***hy= new Double_t**[1056];
	Double_t ***hz= new Double_t**[1056];
//	double hx[1056][21][21],hy[1056][21][21],hz[1056][21][21];
	for(int i=0;i<1056;++i)
	{
		xcart[i]= new Double_t*[21];
		ycart[i]= new Double_t*[21];
		zcart[i]= new Double_t*[21];
		hx[i]= new Double_t*[21];
		hy[i]= new Double_t*[21];
		hz[i]= new Double_t*[21];
		
		for(int j=0;j<21;++j)
		{
			xcart[i][j]= new Double_t[21];
			ycart[i][j]= new Double_t[21];
			zcart[i][j]= new Double_t[21];
			hx[i][j]= new Double_t[21];
			hy[i][j]= new Double_t[21];
			hz[i][j]= new Double_t[21];
		}
		
	}
//	hfield map. To be read in from file.
	double avtheta=0;
//	Something I use to help check for errors
	input.open("Bcoil500.dat");
//	Open input file. Actually, this isn't used until after the quad. So, while the code is currently reading it out, its not yet being used. MB
	while(!input.eof()) //Read in data until the end of file is reached
	{
		input>>xf[num]>>yf[num]>>zf[num]>>Hx[num]>>Hy[num]>>Hz[num];
		//Put that row of data into these variables in this order. If the number of variables here equals the number of variables per line in the input file, this will read in exactly one row of data per loop. And the number of loops will equal the number of rows of data in the file. MB
//		cout<<xf[num]<<endl;
//		To be honest since this input data isn't being used yet, I'm not 100% certain that I will use this organizational scheme. So it seems pointless to explain it right now.
		if(num==0)
		{
			xcart[inpi][inpj][inpk]=xf[num];
			ycart[inpi][inpj][inpk]=yf[num];
			zcart[inpi][inpj][inpk]=zf[num];
			hx[inpi][inpj][inpk]=Hx[num];
			hy[inpi][inpj][inpk]=Hy[num];
			hz[inpi][inpj][inpk]=Hz[num];
//			cout<<xcart[inpi][inpj][inpk]<<endl;
			inpk++;
		}

		else if(num>0)
		{
			if(xf[num]!=xf[num-1])
			{
				inpj=0;
				inpk=0;
				inpi++;
				xcart[inpi][inpj][inpk]=xf[num];
				ycart[inpi][inpj][inpk]=yf[num];
				zcart[inpi][inpj][inpk]=zf[num];
				hx[inpi][inpj][inpk]=Hx[num];
				hy[inpi][inpj][inpk]=Hy[num];
				hz[inpi][inpj][inpk]=Hz[num];
//				cout<<xcart[inpi][inpj][inpk]<<endl;
			}

			else if(yf[num]!=yf[num-1])
			{
				inpk=0;
				inpj++;
				xcart[inpi][inpj][inpk]=xf[num];
				ycart[inpi][inpj][inpk]=yf[num];
				zcart[inpi][inpj][inpk]=zf[num];
				hx[inpi][inpj][inpk]=Hx[num];
				hy[inpi][inpj][inpk]=Hy[num];
				hz[inpi][inpj][inpk]=Hz[num];
//				cout<<xcart[inpi][inpj][inpk]<<endl;
			}

			else
			{
				xcart[inpi][inpj][inpk]=xf[num];
				ycart[inpi][inpj][inpk]=yf[num];
				zcart[inpi][inpj][inpk]=zf[num];
				hx[inpi][inpj][inpk]=Hx[num];
				hy[inpi][inpj][inpk]=Hy[num];
				hz[inpi][inpj][inpk]=Hz[num];
//				cout<<xcart[inpi][inpj][inpk]<<endl;
				inpk++;
			}

		}

		num++;
	}

	input.close();
//	Close the input file

	int sc=0;
	ifstream scprof;
	scprof.open("sc_profile_MCplate_VLN.txt");
//	scprof.open("sc_profile_GhostRanch_nozzle_VLN.txt");
	
	while(!scprof.eof())
	{
		scprof>>number[sc]>>xsc[sc]>>ysc[sc]>>zsc[sc]>>vxsc[sc]>>vysc[sc]>>vzsc[sc]>>thetasc[sc]>>phisc[sc]>>bounce[sc];
		sc++;
	}

	scprof.close();
	
	int numy=0;
	double thetay=0;

//	Here is where the code really begins. Currently this for loop doesn't actually do anything. I think Genya included it so that he could have the code run multiple tests overnight while iterating over some parameter. But there currently isn't any indication of which parameters that would have been. Perhaps he never got around to using it. Could be useful though. MB
	for(itest=0;itest<1;itest++)
	{
		ktot=0; kskim=0; kaper=0; kquad=0; kgood=0; kentry=0;  khole=0; kback=0; kleak=0;
//		kdiaph=0; kpump=0; kpipe4=0; kcell=0;
		//resets particle numbers to zero at start of for loop, but before start of do while loop.
		sum1=0;		//Again, probably useless.

		do
		{
			avtheta+=getnew(itest,rnoz,xnoz,ynoz,znoz,x,y,z,vx,vy,vz,npol,xmem,ymem,zmem,raper,zaper);
//			cout<<"third"<<endl;
			
			double thetamax=0.01;
			thetamax=pi/2;
			
			do
			{
				vx=vxsc[numy]/1000;
				vy=vysc[numy]/1000;
				vz=vzsc[numy]/1000;
				thetay=thetasc[numy];
				numy++;
				if(numy>scl)
				{numy=0;}
			}while(thetay>thetamax);
				
				
//			out2<<x<<" "<<y<<" "<<z<<" "<<vx<<" "<<vy<<" "<<vz<<" "<<npol<<endl;
//			Generate new particle
			
//			cout<<x<<" "<<y<<" "<<z<<endl;
//			putty<<endl<<x<<" "<<y<<" "<<z;

			flag=moveskim(xskim,yskim,zskim,rskim,znoz,vx,vy,vz,x,y,z);
//			Move the particle through the skimmer and track if/where it hits a wall. MB
//			out2<<x<<" "<<y<<" "<<z<<" "<<vx<<" "<<vy<<" "<<vz<<" "<<npol<<endl;
//			cout<<x<<" "<<y<<" "<<z<<endl;
//			putty<<endl<<x<<" "<<y<<" "<<z;

			if(flag>0)	//If particle did not hit skimmer. MB
			{
				flag=moveaper(x,y,z,vx,vy,vz,xaper,yaper,zaper,raper,zskim);
//				Move particle through aperture. MB
//				out2<<x<<" "<<y<<" "<<z<<" "<<vx<<" "<<vy<<" "<<vz<<" "<<npol<<endl;
			}
//			cout<<x<<" "<<y<<" "<<z<<endl;

			if(flag>0)	//If particle didn't hit aperture. MB
			{
				out2<<x<<" "<<y<<" "<<z<<" "<<vx<<" "<<vy<<" "<<vz<<endl;
				flag=move4pol(x,y,z,vx,vy,vz,xtube,ytube,ztube,rtube,ltube,zaper,npol,out2);
//				Mover particle through quad. MB
			}
//			cout<<x<<" "<<y<<" "<<z<<endl;

			if(flag==-3)
			{
				do
				{
//					cout<<sqrt(x*x+y*y+z*z)<<" "<<vx<<" "<<vy<<" "<<vz<<endl;
			
					flag=movetube(x,y,z,vx,vy,vz,xtube,ytube,ztube,rtube,ltube,zaper,npol,out2);
								
//					track<<endl<<i<<" "<<x<<" "<<y<<" "<<z;
					if(flag==10)
					{
						cout<<"MAJOR ERROR"<<endl;
						break;
					}
/*			
					if(flag==0-1 && z>=ztube+ltube-0.000001*ltube)
					{
						exitvrho=sqrt(vx*vx+vy*vy);
						exittheta[kpass]=atan2(exitvrho,vz);
						exitphi[kpass]=atan2(vy,vx);
//*				
//						output<<endl<<i<<" "<<x<<" "<<y<<" "<<z<<" "<<vx<<" "<<vy<<" "<<vz<<" "<<exittheta[kpass]<<" "<<exitphi[kpass]<<" "<<bounce[i];
/*				
						xstore[kpass]=x;
						ystore[kpass]=y;
						zstore[kpass]=z;
						vxstore[kpass]=vx;
						vystore[kpass]=vy;
						vzstore[kpass]=vz;
//*/
//						cout<<sqrt(x*x+y*y)<<"	"<<x<<"	"<<y<<"	"<<z<<" "<<exittheta[kpass]*180/pi<<" "<<exitphi[kpass]*180/pi<<"	"<<bounce[i]<<" "<<kpass<<endl;
//*/				
//						kpass++;
//						storenum++;
//					}
//*/			
					if(z<=ztube+0.0000000000001 && bounce[i]>0)
					{
						flag=0-2;
//						kfail++;
//						storenum++;
//						cout<<kfail<<"	"<<kpass<<"	"<<bounce[i]<<endl;
					}
			
					if(flag>1)
					{
//						errorout<<endl<<"Error: "<<dterror+dtrerror<<":"<<endl;
						if(flag==2)
						{
//							errorout<<"Infinite Velocity: "<<dterror<<":"<<endl<<i<<" "<<x<<" "<<y<<" "<<z<<" "<<vx<<" "<<vy<<" "<<vz<<" "<<bounce[i];	
							cout<<"Error Infinite Velocity:	"<<x<<" "<<y<<" "<<z<<" "<<vx<<" "<<vy<<" "<<vz<<endl;
						}
				
						if(flag==3)
						{
//							errorout<<"Negative Time: "<<dtrerror<<":"<<endl<<i<<" "<<x<<" "<<y<<" "<<z<<" "<<vx<<" "<<vy<<" "<<vz<<" "<<bounce[i];
							cout<<"Error Negative Time:	"<<x<<" "<<y<<" "<<z<<" "<<vx<<" "<<vy<<" "<<vz<<endl;
						}
			
						break;
					}
			
					if(sqrt(x*x+y*y)>rtube+0.000001*rtube)
					{
//						dterror++;
//						errorout<<endl<<"Error: "<<dterror+dtrerror<<":"<<endl;
//						errorout<<"Exceeded tube radius: "<<dterror<<":"<<endl<<i<<" "<<x<<" "<<y<<" "<<z<<" "<<vx<<" "<<vy<<" "<<vz<<" "<<bounce[i]<<" "<<flag;
						cout<<"Error radius:	"<<sqrt(x*x+y*y)<<"	"<<rtube<<endl;
						break;
					}
			
				}while(flag>=0);
				
				if(flag==-1){flag=1;}
				if(flag==-2){flag=0-6;}
				
			}

			if (flag >0)	//If particle didn't hit anything.
			{
				kgood++ ;	   /* went through q-pole ballistic */
				polarization[polcount]=npol;	//Particle made it through. Save polarization. MB
/*
				My Stuff. MB
*/
//				cout<<"4pole: "<<x<<" "<<y<<" "<<z<<endl;
				output<<endl<<x<<" "<<y<<" "<<z<<" "<<vx<<" "<<vy<<" "<<vz<<" "<<npol;
//				output<<endl<<x<<" "<<y<<" "<<z;
				steve<<endl<<x<<" "<<y<<" "<<z<<" "<<vx<<" "<<vy<<" "<<vz;
//				Save beam profile at end of quad. MB

				zRGA=ztube+ltube+50.8;		//Position of US RGA. MB
				movefree(x,y,z,vx,vy,vz,zRGA,ztube,ltube);
//				Move to US RGA. MB
//				out2<<x<<" "<<y<<" "<<z<<" "<<vx<<" "<<vy<<" "<<vz<<" "<<npol<<endl;
				xstore[polcount]=x;
				ystore[polcount]=y;
				zstore[polcount]=z;
//				Store position of particle. MB
				output<<" "<<x<<" "<<y<<" "<<z<<" "<<vx<<" "<<vy<<" "<<vz<<" "<<npol;
//				output<<" "<<x<<" "<<y<<" "<<z;
//				Save beam profile at US RGA. MB.
//				cout<<"free: "<<x<<" "<<y<<" "<<z<<endl;

				zRGA=ztube+ltube+982-zRGA;		//Position of DS RGA
				movefree(x,y,z,vx,vy,vz,zRGA,ztube,ltube);
//				Move to DS RGA. MB
//				out2<<x<<" "<<y<<" "<<z<<" "<<vx<<" "<<vy<<" "<<vz<<" "<<npol<<endl;
				output<<" "<<x<<" "<<y<<" "<<z<<" "<<vx<<" "<<vy<<" "<<vz<<" "<<npol;
//				output<<" "<<x<<" "<<y<<" "<<z;
//				Save beam profile at DS RGA. MB
/*
				End my stuff. MB
*/
//				cout<<polcount<<" "<<polarization[polcount]<<" "<<npol<<endl;
				polcount++;	//Increase iterator for next particle. MB
			}
			if (flag == 0) kskim++ ;	   /* hit skimmer */
			if (flag == -1) kaper++ ;   /* hit aperture */
			if (flag == -2) kentry++ ;  /* hit entrance into q-pole */
			if (flag == -3) khole++ ;   /* pumped through the pumping holes */
			if (flag == -4) kback++ ;   /* pumped through the entrance into q-pole */	//This hasn't been included yet. MB
			if (flag == -5) kleak++ ;   /* leaked through the exit of q-pole */	//This hasn't been included yet. MB
			if (flag == -6) kreverse++;
			kquad=khole+kback+kleak+kgood+kreverse;
//			Total number of particles that enter the quad. MB

//			if(flag>0)
//			{output<<endl<<x<<" "<<y<<" "<<z<<" "<<vx<<" "<<vy<<" "<<vz<<" "<<xmem<<" "<<ymem;}	//Probably just for error checking. MB
//			cout<<ktot<<endl;
			ktot++;		//We have a new particle; increase the count. MB
		}while(ktot<ntot);

		sum1=sum1/ktot;	//Again this will never be anything other than zero. MB
	}

	cout<<"Total: "<<ktot<<". Tnoz="<<tnoz<<endl;
//	Output total number of particle and nozzle temp. MB
	cout<<"Percent lost on first skimmer: "<<kskim*100.0/ktot<<"%"<<endl<<"Percent through first skimmer: "<<(ktot-kskim)*100.0/ktot<<"%"<<endl;
	cout<<"Percent lost on second skimmer: "<<kaper*100.0/ktot<<"%"<<endl<<"Percent through second skimmer: "<<(ktot-kskim-kaper)*100.0/ktot<<"%"<<endl;
	cout<<"Percent lost to Quad entryway: "<<kentry*100.0/ktot<<"%"<<endl<<"Percent \"rejected\" by Quad field: "<<khole*100.0/ktot<<"%"<<endl;
	cout<<"Percent lost to Quad: "<<(kentry+khole)*100.0/ktot<<"%"<<endl<<"Percent through Quad: "<<(ktot-kskim-kaper-kentry-khole)*100.0/ktot<<"%"<<endl;
//	Output percentage of particles lost to each component and the perctange that were transmitted through each component. MB
	cout<<"average theta: "<<avtheta/ktot<<endl;
//	Output average intial angle for particles. Helps for error checking and showing if the initial angle is artificially constrained for some other analysis. Artificially constraining the angle produces more particles that make it through the quad, ignoring the ones that do not. This allows for more statistics to be gathered on the output polarization by not wasting time on particles that are guaranteed to be lost anyway. The initial angle can be artificially constrained by changing thetamax in the getnew() function. MB

/*
	polsum=0;
	
	for(i=0;i<ipol;i++)
	{
		polsum=polsum+polmem[i];
	}

	polavg=polsum/ipol;
	polsum=0;

	for(i=0;i<ipol;i++)
	{
		del=polavg-polmem[i];
		polsum=polsum+del*del;
	}

	dumdum=sqrt(polsum/ipol);
	polerr=dumdum/sqrt(ipol);
*/	//This is for use after guiding particles into the cell. Not included in this code yet.  MB


	polsum=0;
	
	for(i=0;i<polcount;i++)
	{
		polsum=polsum+polarization[i];
//		cout<<i<<" "<<polarization[i]<<endl;
	}

	polavg=polsum/polcount;
	polsum=0;

	for(i=0;i<polcount;i++)
	{
		del=polavg-polarization[i];
		polsum=polsum+del*del;
	}

	dumdum=sqrt(polsum/polcount);
	polerr=dumdum/sqrt(polcount);
//	Checks polarization just outside the quad. MB



	cout<<"Pol="<<polavg<<" : "<<polerr<<endl;
	cout<<"At x="<<xquad<<endl;
	cout<<"Spread="<<dumdum<<endl;
//	Outputs results of whichever polarization checking code is being used. Make sure that the one not in use is always commented out. MB
	output.close();



/*
Everything from this point forward was written by me. The section that displays figures uses ROOT specific code. If you do not use ROOT, you will need to comment that out. MB
*/
	double xbins=20;
	double ybins=20;
	double xmax, ymax, xmin, ymin, xblip, yblip, x1, y1, x2, y2;
//	cout<<"Here"<<endl;
	double xgram[10000]={};
	double ygram[10000]={};
	double xbox[10000], ybox[10000];
	double check=0;
	double dcheck=0;
	xmax=xstore[0];
	xmin=xstore[0];

	for(int p=0;p<polcount;p++)
	{
		if(xstore[p]>xmax)
		{xmax=xstore[p];}
		if(xstore[p]<xmin)
		{xmin=xstore[p];}
	}

	xblip=(xmax-xmin)/xbins;
	x1=xmin;
	x2=xmin+xblip;
//	cout<<xblip<<" "<<x1<<" "<<x2<<endl;
//	cout<<xmin<<" "<<xmax<<endl;

	for(int p=0;p<xbins;p++)
	{
		xbox[p]=(x1+x2)/2;
		x1+=xblip;
		x2+=xblip;
//		cout<<xblip<<" "<<xmin<<" "<<xmax<<" "<<x1<<" "<<x2<<endl;
	}

	for(int p=0;p<polcount;p++)
	{
		x1=xmin;
		x2=xmin+xblip;

		for(int px=0;px<xbins;px++)
		{
//			cout<<x1<<" "<<x2<<" "<<xstore[p]<<endl;

			if(xstore[p]>=x1 && xstore[p]<=x2)
			{
				xgram[px]++;
				check++;
				break;
			}

//			if(px==0){dcheck++;}
			x1+=xblip;
			x2+=xblip;
		}

	}

	cout<<polcount<<" "<<check<<endl;



	TCanvas *ABS= new TCanvas("ABS","Predicted Beam Profile",10,10,1280,620);
	ABS->SetGridx(0);
	ABS->SetGridy(0);
	ABS->SetFillColor(10);
	ABS->GetFrame()->SetFillColor(-10);
	ABS->SetBorderMode(0);
	ABS->SetBorderSize(2);
//	ABS->Divide(2,2);
//	ABS->Divide(1,2);

	gStyle->SetOptFit();
	gStyle->SetTitleSize(.08,"t");


//		ABS->cd(1);
		auto gr0=new TGraphErrors(xbins,xbox,xgram,0,0);		

		TF1 *fitG=new TF1("fitg",fitg,-10,10,4);
//		fitG->SetParameter(0,200);
		fitG->SetParameter(0,2000);
		fitG->SetParameter(1,0.5);
		fitG->SetParameter(2,1);
		fitG->SetParameter(3,0);
		fitG->SetParName(0,"Peak");
		fitG->SetParName(1,"Center");
		fitG->SetParName(2,"Sigma");
		fitG->SetParName(3,"Background");

		TF1 *fitGC=new TF1("fitgc",fitgc,-10,10,6);
		fitGC->SetParameter(0,200);
		fitGC->SetParameter(1,0.5);
		fitGC->SetParameter(2,1);
		fitGC->SetParameter(3,0);
		fitGC->SetParameter(4,1/28);
		fitGC->SetParameter(5,0.5);
		fitGC->SetParName(0,"Peak");
		fitGC->SetParName(1,"Center");
		fitGC->SetParName(2,"Sigma");
		fitGC->SetParName(3,"Amplitude");
		fitGC->SetParName(4,"Wavelength");
		fitGC->SetParName(5,"Phase");

		gr0->SetMarkerStyle(20);
		gr0->SetMarkerColor(1);
		gr0->SetMarkerSize(1.0);
		gr0->SetLineWidth(1);
		gr0->SetLineColor(1);

		gr0->SetTitle("Predicted Beam Profile zRGA: 50.8");
//		gr0->GetTitle()->SetTitleSize(.08);
		gr0->GetXaxis()->SetTitle("RGA position (mm)");
		gr0->GetYaxis()->SetTitle("Number of Particles");
		gr0->GetXaxis()->SetTitleSize(.05);
		gr0->GetYaxis()->SetTitleSize(.05);
		gr0->GetXaxis()->CenterTitle(1);
		gr0->GetYaxis()->CenterTitle(1);
		gr0->GetXaxis()->SetTitleOffset(.91);
		gr0->GetYaxis()->SetTitleOffset(.9);
		gr0->GetXaxis()->SetLabelSize(.05);
		gr0->GetYaxis()->SetLabelSize(.05);
		gr0->Draw("ap");
		gr0->Fit("fitg","r");

		double peak, errpeak, center, errcenter, sigma, errsigma, back, errback;
		peak=fitG->GetParameter(0);
		center=fitG->GetParameter(1);
		sigma=fitG->GetParameter(2);
		back=fitG->GetParameter(3);
		errpeak=fitG->GetParError(0);
		errcenter=fitG->GetParError(1);
		errsigma=fitG->GetParError(2);
		errback=fitG->GetParError(3);

		double integral=fitG->Integral(35,90);
		double errintegral=fitG->IntegralError(35,90);




	cout<<endl<<"Area: "<<integral<<" : "<<errintegral<<" (torr*mm^2)"<<endl<<endl;

	cout<<"peak: "<<peak<<" : "<<errpeak<<" (torr)"<<endl;
	cout<<"center: "<<center<<" : "<<errcenter<<" (mm)"<<endl;
	cout<<"sigma: "<<sigma<<" : "<<errsigma<<" (mm)"<<endl;
	cout<<"background: "<<back<<" : "<<errback<<" (torr)"<<endl;


}


double fitg(double *xf, double *par){

	double zf=xf[0];

	double gauss=par[0]*exp(-0.5*pow(((zf-par[1])/par[2]),2))+par[3];
	return gauss;
}

double fitgc(double *xf, double *par){

	double zf=xf[0];

	double gauss=par[0]*exp(-0.5*pow(((zf-par[1])/par[2]),2))+par[3]*cos(par[4]*zf-par[5]);
	return gauss;
}




















