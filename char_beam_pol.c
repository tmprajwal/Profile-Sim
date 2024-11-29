//Maximum array is 1045769
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<algorithm>

#define pi 3.14159
/*
double fitg(double *xf, double *par);

double fitgc(double *xf, double *par);

double fitgg(double *xf, double *par);

double cauchy(double *xf, double *par);

double sech(double *xf, double *par);

double fite(double *xf, double *par);

double fitq(double *xf, double *par);

double poly(double *xf, double *par);
*/
void char_beam_pol()
{
	Double_t *xf= new Double_t[10000000];
	Double_t *yf= new Double_t[10000000];
	Double_t *zf= new Double_t[10000000];
	Double_t *xqstore= new Double_t[10000000];
	Double_t *yqstore= new Double_t[10000000];
	Double_t *zqstore= new Double_t[10000000];
	Double_t *xustore= new Double_t[10000000];
	Double_t *yustore= new Double_t[10000000];
	Double_t *zustore= new Double_t[10000000];
	Double_t *xdstore= new Double_t[10000000];
	Double_t *ydstore= new Double_t[10000000];
	Double_t *zdstore= new Double_t[10000000];
	Double_t *rho= new Double_t[10000000];
	Double_t *pol= new Double_t[10000000];
	Double_t *r= new Double_t[10000000];
	int polcount=0;
	
	Double_t *vx= new Double_t[10000000];
	Double_t *vy= new Double_t[10000000];
	Double_t *vz= new Double_t[10000000];

	Double_t *xstore= new Double_t[10000000];
	Double_t *ystore= new Double_t[10000000];
	Double_t *zstore= new Double_t[10000000];
	Double_t *vratio= new Double_t[10000000];
	Double_t *vR= new Double_t[10000000];
	Double_t *vrho= new Double_t[10000000];
	Double_t *theta= new Double_t[10000000];
	Double_t *vtheta= new Double_t[10000000];
	Double_t *phi= new Double_t[10000000];
	Double_t *bounce= new Double_t[10000000];
	
	Double_t *R= new Double_t[10000000];
	int *number= new int[10000000];
	
	double time=0;
	double zcharf=57.41;
	double rchar=1.651;
	int survive=0;
	double polarization=0;

	ifstream input;
//	input.open("profile_quad_outlet_pos_vel_pol_1K.txt");
//	input.open("profile_quad_outlet_pos_vel_pol_1K_proper.txt");
//	input.open("profile_quad_outlet_pos_vel_pol_1p5K.txt");
	input.open("profile_quad_outlet_pos_vel_pol_2K.txt");

	while(!input.eof())
	{	
		input>>xstore[polcount]>>ystore[polcount]>>zstore[polcount]>>vx[polcount]>>vy[polcount]>>vz[polcount]>>pol[polcount];
		vrho[polcount]=sqrt(vx[polcount]*vx[polcount]+vy[polcount]*vy[polcount]);
		theta[polcount]=(atan2(vrho[polcount],vz[polcount]))*180/pi;
		
		if(vz[polcount]>0)
		{
			time=zcharf/vz[polcount];
			xf[polcount]=xstore[polcount]+vx[polcount]*time;
			yf[polcount]=ystore[polcount]+vy[polcount]*time;
			zf[polcount]=zcharf;
			rho[polcount]=sqrt(xf[polcount]*xf[polcount]+yf[polcount]*yf[polcount]);
		}

		if(vz[polcount]<=0)
		{rho[polcount]=100000;}
		
		if(rho[polcount]<rchar)
		{
			polarization+=pol[polcount];
			survive++;
		}

		polcount++;
	}
	
	polarization=polarization/survive;
	cout<<polarization*(0-100)<<"%"<<endl;
/*	
	xmax=theta[0];
	xmin=theta[0];

	for(int p=0;p<polcount;p++)
	{		
		if(theta[p]>xmax)
		{xmax=theta[p];}
		if(theta[p]<xmin)
		{xmin=theta[p];}
	}
	
	cout<<xmin<<" "<<xmax<<endl;
	
//	xmin=-90; xmax=90;
	xmin=0; xmax=5;

	xblip=(xmax-xmin)/xbins;
//	xblip=3;
	xbins=(xmax-xmin)/xblip;
	x1=xmin;
	x2=xmin+xblip;

	for(int p=0;p<xbins;p++)
	{
		xbox[p]=(x1+x2)/2;
		errxbox[p]=xblip/1000;
		x1+=xblip;
		x2+=xblip;
	}

	check=0;

	for(int p=0;p<polcount;p++)
	{
		x1=xmin;
		x2=xmin+xblip;

		for(int px=0;px<xbins;px++)
		{

//			if(ystore[p]>=(0-3.175) && ystore[p]<=3.175)
//			if(ystore[p]>=(0-1) && ystore[p]<=1)
//			{
//				center1=(x1+x2)/2;
//				if(xstore[p]>=center1-3.175 && xstore[p]<=center1+3.175)
				if(theta[p]>=x1 && theta[p]<=x2 && pol[p]>0)
				{
//					radius=sqrt((xstore[p]-center1)*(xstore[p]-center1)+ystore[p]*ystore[p]);
//					if(radius<=3.175)
//					{
						xgram[px]++;
						if(theta[p]<1.59/2)			//degrees
						{check++;}
//					}

				}

//			}

			x1+=xblip;
			x2+=xblip;
		}

	}

	cout<<polcount<<" "<<check<<endl<<double(check)/double(polcount)*100<<"%"<<endl;

	double grammax=0;
	
	for(int p=0;p<xbins;p++)
	{
		xgram[p]=xgram[p]/abs(sin(xbox[p]*pi/180));
		if(xgram[p]>grammax)
		{grammax=xgram[p];}	
	}

	for(int p=0;p<xbins;p++)
	{
		xgram[p]=xgram[p]/grammax;
		errxgram[p]=0.0001;
//		cout<<xgram[p]<<" "<<xbox[p]<<endl;
	}
//*/

/*
	TCanvas *ABS= new TCanvas("ABS","Predicted Beam Profile",10,10,1280,620);
	ABS->SetGridx(0);
	ABS->SetGridy(0);
	ABS->SetFillColor(10);
	ABS->GetFrame()->SetFillColor(-10);
	ABS->SetBorderMode(0);
	ABS->SetBorderSize(2);

	gStyle->SetOptFit();
	gStyle->SetTitleSize(.08,"t");



		ABS->cd(1);	
		auto gr0=new TGraphErrors(xbins, xbox, xgram, 0, 0);
//*/	
/*
		TF1 *fitG=new TF1("fitg",fitg,-20,20,4);
		fitG->SetParameter(0,5000);
		fitG->SetParameter(1,0);
		fitG->SetParameter(2,1);
		fitG->SetParameter(3,0);
		fitG->SetParName(0,"Peak");
		fitG->SetParName(1,"Center");
		fitG->SetParName(2,"Sigma");
		fitG->SetParName(3,"Background");
*/
/*
		gr0->SetMarkerStyle(20);
		gr0->SetMarkerColor(1);
		gr0->SetMarkerSize(1.0);
		gr0->SetLineWidth(1);
		gr0->SetLineColor(1);

		gr0->SetTitle("Simulated Quad Exit Profile (Rejected Atoms)");
		gr0->GetXaxis()->SetTitle("Angle (degrees)");
		gr0->GetYaxis()->SetTitle("Normalized Pressure");
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
//*/


	input.close();
	
	
}
/*
double fitg(double *xf, double *par){

	double zf=xf[0];

	double gauss=par[0]*exp(-0.5*pow(((zf-par[1])/par[2]),2))+par[3];
	return gauss;
}

double fitgc(double *xf, double *par){

	double zf=xf[0];

//	double gauss=par[0]*exp(-0.5*((zf-par[1])/par[2])**2)+par[3]*cos(par[4]*zf-par[5]);
	double gauss=par[3]+par[0]*cos(par[1]*zf-par[2]);
	return gauss;
}


double fitgg(double *xf, double *par){

	double zf=xf[0];

	double gauss=par[0]*exp(-0.5*pow(((zf-par[1])/par[2]),2))+par[3]*exp(-0.5*pow(((zf-par[4])/par[5]),2));
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

double fitq(double *xf, double *par)
{
	double zf=xf[0];
	double dexp=par[0]*(zf-par[1])*(zf-par[1]);
	return dexp;
}

double poly(double *xf, double *par)
{
	double zf=xf[0];
	double polynomial=par[0]+par[1]*zf+par[2]*zf*zf+par[3]*zf*zf*zf;
	return polynomial;
}
*/
