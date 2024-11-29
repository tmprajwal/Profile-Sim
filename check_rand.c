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
#define tnoz 1.3;		//Nozzle temperature, K
#define ntot 1000000;		//Number of particles
//#define ntot 50;
#define stept 3e-5;		//Time steps, in seconds
#define gradq 112.5;		//Gradient in quads, T/m (9 kG on the tip)
//Definitions are available to all functions. However, they generally cannot/should not be altered within the code. This makes them a good substitute for global variables, but only for variables that are expected to remain constant throughout the code.



void check_rand()	/* generates random number between min and max */
{
	double xstore[100000];
	double xbox[10000], errxbox[10000], errxgram[10000];
	double xgram[10000]={};
	double max=pi/2;
	double min=0;
	double xbins=100;
	int check=0;
	int polcount=0;
	double xmin, xmax, xblip, x1, x2;
	double r, comp, res, theta;


	for(int i=0;i<100000;i++)
	{
		res=min+(double(max)-min)*rand()/RAND_MAX;
//		Uses predefined rand() function and picks a random number between min and max. "res" is set to this value. "RAND_MAX" is the maximum value of the rand() function. MB

		do
		{
			r=min+(double(max)-min)*rand()/RAND_MAX;
			comp=min+(double(max)-min)*rand()/RAND_MAX;
//			The two lines above each generate a random number according to the randgen function (defined elsewhere in code) with a minimum value of 0 and a maximum value of r0. MB
		} while(comp>r);
//	The "do" loop is performed until r is equal to or greater than comp. MB

		do
		{
			theta=min+(double(max)-min)*rand()/RAND_MAX;
//			Obtain random value for theta using randgen(). MB
			comp=min+(double(0.5)-min)*rand()/RAND_MAX;
//			Obtain random value for comp using randgen(). MB
		} while(comp>sin(theta)*cos(theta));

		xstore[i]=r;
		polcount++;
	}

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

	for(int p=0;p<xbins;p++)
	{
		xbox[p]=(x1+x2)/2;
		errxbox[p]=xblip/1000;
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

	for(int p=0;p<xbins;p++)
	{errxgram[p]=sqrt(xgram[p]);}





	TCanvas *ABS= new TCanvas("ABS","Rand Dis",10,10,1280,620);
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
//		gr0=new TGraphErrors(xbins, xbox, xgram, 0, 0);
		gr0=new TGraphErrors(xbins, xbox, xgram, errxbox, errxgram);		

		fitp=new TF1("poly","pol1",0,10);
/*
		fitG=new TF1("fitg","gaus",-60,60);
		fitG->SetParameter(0,1800);
		fitG->SetParameter(1,0);
		fitG->SetParameter(2,2);
//*/
///*
		fitG=new TF1("fitg",fitg,-30,30,4);
//		fitG->SetParameter(0,200);
		fitG->SetParameter(0,1800);
		fitG->SetParameter(1,0.5);
		fitG->SetParameter(2,4);
		fitG->SetParameter(3,0);
		fitG->SetParName(0,"Peak");
		fitG->SetParName(1,"Center");
		fitG->SetParName(2,"Sigma");
		fitG->SetParName(3,"Background");
//*/
		fitGC=new TF1("fitgc",fitgc,-110,130,6);
		fitGC->SetParameter(0,450);
		fitGC->SetParameter(1,0);
		fitGC->SetParameter(2,2);
		fitGC->SetParameter(3,50);
		fitGC->SetParameter(4,1/28);
		fitGC->SetParameter(5,0);
		fitGC->SetParName(0,"Peak");
		fitGC->SetParName(1,"Center");
		fitGC->SetParName(2,"Sigma");
		fitGC->SetParName(3,"Amplitude");
		fitGC->SetParName(4,"Wavelength");
		fitGC->SetParName(5,"Phase");

		Cauchy=new TF1("cauchy",cauchy,-7,7,4);
//		Cauchy->SetParameter(0,200);
		Cauchy->SetParameter(0,1800);
		Cauchy->SetParameter(1,0.1);
		Cauchy->SetParameter(2,4.4);
		Cauchy->SetParameter(3,0);
		Cauchy->SetParName(0,"Peak");
		Cauchy->SetParName(1,"Center");
		Cauchy->SetParName(2,"Sigma");
		Cauchy->SetParName(3,"Background");

		fitE=new TF1("fite",fite,-7,0.1,3);
//		fitE->SetParameter(0,200);
		fitE->SetParameter(0,1000);
		fitE->SetParameter(1,1);
		fitE->SetParameter(2,0);
		fitE->SetParameter(3,0);
		fitE->SetParName(0,"Coefficient");
		fitE->SetParName(1,"Half");
		fitE->SetParName(2,"Offset");
		fitE->SetParName(3,"Background");

		gr0->SetMarkerStyle(20);
		gr0->SetMarkerColor(1);
		gr0->SetMarkerSize(1.0);
		gr0->SetLineWidth(0.5);
		gr0->SetLineColor(1);

		gr0->SetTitle("Rand Max Dis");
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
//		gr0->Fit("fitg","r");
//		gr0->Fit("cauchy","r");
//		gr0->Fit("fite","r");
		gr0->Fit("poly","r");

		double peak, errpeak, center, errcenter, sigma, errsigma, back, errback;
		peak=fitG->GetParameter(0);
		center=fitG->GetParameter(1);
		sigma=fitG->GetParameter(2);
		back=fitG->GetParameter(3);
		errpeak=fitG->GetParError(0);
		errcenter=fitG->GetParError(1);
		errsigma=fitG->GetParError(2);
		errback=fitG->GetParError(3);
/*
		double integral=fitG->Integral(35,90);
		double errintegral=fitG->IntegralError(35,90);




	cout<<endl<<"Area: "<<integral<<" : "<<errintegral<<" (torr*mm^2)"<<endl<<endl;

	cout<<"peak: "<<peak<<" : "<<errpeak<<" (torr)"<<endl;
	cout<<"center: "<<center<<" : "<<errcenter<<" (mm)"<<endl;
	cout<<"sigma: "<<sigma<<" : "<<errsigma<<" (mm)"<<endl;
	cout<<"background: "<<back<<" : "<<errback<<" (torr)"<<endl;
//*/

}


double fitg(double *xf, double *par){

	double zf=xf[0];

	double gauss=par[0]*exp(-0.5*((zf-par[1])/par[2])**2)+par[3];
	return gauss;
}

double fitgc(double *xf, double *par){

	double zf=xf[0];

	double gauss=par[0]*exp(-0.5*((zf-par[1])/par[2])**2)+par[3]*cos(par[4]*zf-par[5]);
	return gauss;
}


double fitgg(double *xf, double *par){

	double zf=xf[0];

	double gauss=par[0]*exp(-0.5*((zf-par[1])/par[2])**2)+par[3]*exp(-0.5*((zf-par[4])/par[5])**2);
	return gauss;
}

double cauchy(double *xf, double *par)
{
	double zf=xf[0];

	double cauchy=par[3]+(par[0]/(pi*par[2]))*(par[2]*par[2]/((zf-par[1])*(zf-par[1])+par[2]*par[2]));
	return cauchy;
}

double sech(double *xf, double *par)
{
	double zf=xf[0];

	double sech=par[2]+par[0]*2/(exp(par[1]*zf)+exp(par[1]*zf));
	return sech;
}

double fite(double *xf, double *par)
{
	double zf=xf[0];
	double dexp=par[3]+par[0]*exp(par[1]*zf+par[2]);
	return dexp;
}





