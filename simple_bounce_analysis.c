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
void simple_bounce_analysis()
{
	Double_t *num= new Double_t[10000000];
	Double_t *vtheta= new Double_t[10000000];
	Double_t *phi= new Double_t[10000000];
	Double_t *bounce= new Double_t[10000000];
	int *state= new int[10000000];
	
	Double_t *yqstore= new Double_t[10000000];
	Double_t *zqstore= new Double_t[10000000];
	Double_t *xustore= new Double_t[10000000];
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
	
	Double_t *R= new Double_t[10000000];
	
//	char *pass= new char[10000000];
	string pass;

	double xbins=1000;
	double ybins=20;
	int pnum=4;
	double xmax, ymax, xmin, ymin, xblip, yblip, x1, y1, x2, y2;
	double xmux, xmun, xblup, xu1, xu2;
	double xgram[10000]={};
	double ygram[10000]={};
	double rgram[10000]={};
	double errxgram[10000], errygram[10000], errrgram[10000];
	double xbox[10000], ybox[10000], errxbox[10000], errybox[10000];
	double rbox[10000], errrbox[10000];
	double xslot[10000]={};
	double errxslot[10000];
	double xhist[10000], errxhist[10000];
	double check=0;
	double dcheck=0;
	int count=0;
	double radius, center1;
	double bounce_avg=0;

	ifstream input;
//	input.open("qt_bounce.txt");
	input.open("blocked_bounce.txt");
	
	cout<<"Data Loaded"<<endl;

	while(!input.eof())
	{
		input>>pass>>state[polcount]>>bounce[polcount];
		bounce_avg+=bounce[polcount];
		polcount++;
	}
	
	cout<<"Data Acquired"<<endl;
	
	bounce_avg=bounce_avg/polcount;
	
	cout<<"Average # of bounces: "<<bounce_avg<<endl;
	
	polcount=polcount-1;

	xmax=xstore[0];
	xmin=xstore[0];
	
	xmax=theta[0];
	xmin=theta[0];
//*
	for(int p=0;p<polcount;p++)
	{
		if(bounce[p]>xmax)
		{xmax=bounce[p];}
		if(bounce[p]<xmin)
		{xmin=bounce[p];}
	}
//*/

//*	
	cout<<xmin<<" "<<xmax<<endl;
	
//	xmin=-90; xmax=90;
//	xmin=-40; xmax=40;

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
			vtheta[p]=bounce[p];

			if(vtheta[p]>=x1 && vtheta[p]<=x2)
			{					
				xgram[px]++;
				check++;
			}
			
			x1+=xblip;
			x2+=xblip;

		}


	}

//	cout<<polcount<<" "<<check<<endl;
/*
	for(int p=0;p<xbins;p++)
	{
		xgram[p]=xgram[p]/(40*abs(sin(xbox[p]*pi/180)));
//		xgram[p]=xgram[p]/abs(xbox[p]);
//		errxgram[p]=400;
//		cout<<xgram[p]<<" "<<xbox[p]<<endl;
	}
//*/
//*/

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
//		auto gr0=new TGraphErrors(polcount, theta, bounce,0,0);
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

		gr0->SetMarkerStyle(20);
		gr0->SetMarkerColor(1);
		gr0->SetMarkerSize(1.0);
		gr0->SetLineWidth(1);
		gr0->SetLineColor(1);
/*
		gr0->SetTitle("Angular Velocity Distribution (Single Channel)");
		gr0->GetXaxis()->SetTitle("Angle with cylindrical axis (degrees)");
*/
		gr0->SetTitle("Bounce Distribution (Single Channel)");
		gr0->GetXaxis()->SetTitle("Number of Bounces");
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
