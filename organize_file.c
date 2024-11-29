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
void organize_file()
{
	Double_t *num= new Double_t[10000000];
	Double_t *vtheta= new Double_t[10000000];
	Double_t *phi= new Double_t[10000000];
	Double_t *bounce= new Double_t[10000000];
	
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

	double xbins=100;
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
	int skip;

	ifstream input;
	input.open("profile_1K_solidworks_dims.txt");
	
	ofstream output;
	output.open("quad_outlet_profile_1K_solidworks_dims.txt");

	while(!input.eof())
	{
		input>>xstore[polcount]>>ystore[polcount]>>zstore[polcount]>>vx[polcount]>>vy[polcount]>>vz[polcount]>>skip;

		output<<xstore[polcount]<<" "<<ystore[polcount]<<" "<<zstore[polcount]<<" "<<vx[polcount]<<" "<<vy[polcount]<<" "<<vz[polcount]<<endl;
		polcount++;
	}
	
	input.close();
	output.close();
	
	
}

