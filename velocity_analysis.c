//Maximum array is 1045769
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<algorithm>

#define pi 3.14159

double fitg(double *xf, double *par);

double fitgc(double *xf, double *par);

double fitgg(double *xf, double *par);

double cauchy(double *xf, double *par);

double sech(double *xf, double *par);

double fite(double *xf, double *par);

double fitq(double *xf, double *par);

double poly(double *xf, double *par);

void velocity_analysis()
{
	Double_t *xi= new Double_t[10000000];
	Double_t *yi= new Double_t[10000000];
	Double_t *zi= new Double_t[10000000];
	Double_t *xqstore= new Double_t[10000000];
	Double_t *yqstore= new Double_t[10000000];
	Double_t *zqstore= new Double_t[10000000];
	Double_t *vmag= new Double_t[10000000];
	Double_t *vqx= new Double_t[10000000];
	Double_t *vqy= new Double_t[10000000];
	Double_t *vqz= new Double_t[10000000];

	Double_t *theta= new Double_t[10000000];
	Double_t *phi= new Double_t[10000000];
	Double_t *hittime= new Double_t[10000000];
/*
	Double_t *pol= new Double_t[10000000];
	Double_t *r= new Double_t[10000000];
*/
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
	
	Double_t *R= new Double_t[10000000];
	Double_t *Rq= new Double_t[10000000];
	Double_t *Ru= new Double_t[10000000];
	Double_t *Rd= new Double_t[10000000];
	
	Double_t *pol= new Double_t[10000000];
	
	int tempcount=0;

	double vbins=100;
	double ybins=20;
	int pnum=4;
	double vmax, ymax, vmin, ymin, vblip, yblip, v1, y1, v2, y2;
	double xmux, xmun, xblup, xu1, xu2;
	double vgram[10000]={};
	double ygram[10000]={};
	double rgram[10000]={};
	double errvgram[10000], errygram[10000], errrgram[10000];
	double vbox[10000], ybox[10000], errvbox[10000], errybox[10000];
	double rbox[10000], errrbox[10000];
	double xslot[10000]={};
	double errxslot[10000];
	double xhist[10000], errxhist[10000];
	double check=0;
	double dcheck=0;
	double polfull=0;
	double polq=0;
	double polus=0;
	double polds=0;
	int count=0;
	int qn=0;
	int un=0;
	int dn=0;
	double radius, center1;
	
	double skip;

	ifstream input;
//	input.open("profile_quad_outlet_pos_vel_pol_1K.txt");
	input.open("profile_quad_outlet_pos_vel_pol_1p5K.txt");
	cout<<"data_loaded"<<endl;

	while(!input.eof())
	{
		input>>xqstore[polcount]>>yqstore[polcount]>>zqstore[polcount]>>vqx[polcount]>>vqy[polcount]>>vqz[polcount]>>pol[polcount];
		vmag[polcount]=sqrt(vqx[polcount]*vqx[polcount]+vqy[polcount]*vqy[polcount]+vqz[polcount]*vqz[polcount]);
		polcount++;
	}

	cout<<"data_read"<<endl;
	cout<<polcount<<endl;

	vmax=vmag[0];
	vmin=vmag[0];

	for(int p=0;p<polcount;p++)
	{
		if(vmag[p]>vmax)
		{vmax=vmag[p];}
		if(vmag[p]<vmin)
		{vmin=vmag[p];}
	}
	
	cout<<vmin<<" "<<vmax<<endl;
	
//	vmin=-40; vmax=40;

	vblip=(vmax-vmin)/vbins;
//	vblip=3;
	vbins=(vmax-vmin)/vblip;
	v1=vmin;
	v2=vmin+vblip;

	for(int p=0;p<vbins;p++)
	{
		vbox[p]=(v1+v2)/2;
		errvbox[p]=vblip/1000;
		v1+=vblip;
		v2+=vblip;
	}

	check=0;

	for(int p=0;p<polcount;p++)
	{
		v1=vmin;
		v2=vmin+vblip;

		for(int pv=0;pv<vbins;pv++)
		{

			if(vmag[p]>=v1 && vmag[p]<v2 && pol[p]>0)
			{
				vgram[pv]++;
				check++;
			}

			v1+=vblip;
			v2+=vblip;
		}

	}

	cout<<polcount<<" "<<check<<endl;
/*
	for(int p=0;p<vbins;p++)
	{
		vgram[p]=vgram[p]/2.72;
		errvgram[p]=400;
//		cout<<xgram[p]<<" "<<xbox[p]<<endl;
	}
*/


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
		auto gr0=new TGraphErrors(vbins, vbox, vgram, 0, 0);	

		TF1 *fitG=new TF1("fitg",fitg,-40,40,3);
		fitG->SetParameter(0,900);
		fitG->SetParameter(1,0);
		fitG->SetParameter(2,4.4);
//		fitG->SetParameter(3,0);
		fitG->SetParName(0,"Peak");
		fitG->SetParName(1,"Center");
		fitG->SetParName(2,"Sigma");
//		fitG->SetParName(3,"Background");


		gr0->SetMarkerStyle(20);
		gr0->SetMarkerColor(1);
		gr0->SetMarkerSize(1.0);
		gr0->SetLineWidth(1);
		gr0->SetLineColor(1);

		gr0->SetTitle("Uniform 1K Temperature");
		gr0->GetXaxis()->SetTitle("velocity (m/s)");
		gr0->GetYaxis()->SetTitle("Number of Particles");
		gr0->GetXaxis()->SetTitleSize(.05);
		gr0->GetYaxis()->SetTitleSize(.05);
		gr0->GetXaxis()->CenterTitle(1);
		gr0->GetYaxis()->CenterTitle(1);
		gr0->GetXaxis()->SetTitleOffset(.91);
		gr0->GetYaxis()->SetTitleOffset(1);
		gr0->GetXaxis()->SetLabelSize(.05);
		gr0->GetYaxis()->SetLabelSize(.05);
//		gr0->GetYaxis()->SetRangeUser(0,3e3);
		gr0->Draw("ap");
//		gr0->Fit("fitg","r");



	input.close();
	
	
}

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
