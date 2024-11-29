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

void wtf()
{
	Double_t *xi= new Double_t[10000000];
	Double_t *yi= new Double_t[10000000];
	Double_t *zi= new Double_t[10000000];
	Double_t *xqstore= new Double_t[10000000];
	Double_t *yqstore= new Double_t[10000000];
	Double_t *zqstore= new Double_t[10000000];
	Double_t *xustore= new Double_t[10000000];
	Double_t *yustore= new Double_t[10000000];
	Double_t *zustore= new Double_t[10000000];
	Double_t *xdstore= new Double_t[10000000];
	Double_t *ydstore= new Double_t[10000000];
	Double_t *zdstore= new Double_t[10000000];
	Double_t *vqx= new Double_t[10000000];
	Double_t *vqy= new Double_t[10000000];
	Double_t *vqz= new Double_t[10000000];
	Double_t *vux= new Double_t[10000000];
	Double_t *vuy= new Double_t[10000000];
	Double_t *vuz= new Double_t[10000000];
	Double_t *vdx= new Double_t[10000000];
	Double_t *vdy= new Double_t[10000000];
	Double_t *vdz= new Double_t[10000000];
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
	
	Double_t *R= new Double_t[10000000];
	Double_t *Rq= new Double_t[10000000];
	Double_t *Ru= new Double_t[10000000];
	Double_t *Rd= new Double_t[10000000];

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
	double polfull=0;
	double polq=0;
	double polus=0;
	double polds=0;
	int count=0;
	int qn=0;
	int un=0;
	int dn=0;
	double radius, center1;

	ifstream input;
//	input.open("profile_MCplate_ideal.txt");
//	input.open("profile_6_23_2020.txt");
//	input.open("profile_old_nozzle_real.txt");
//	input.open("profile_mcplate_ynoz_9mm.txt");
	input.open("profile_temp_2k.txt");
	cout<<"data_loaded"<<endl;

	while(!input.eof())
	{
//		input>>xstore[polcount]>>ystore[polcount]>>zstore[polcount];
//		R[polcount]=sqrt(xstore[polcount]*xstore[polcount]+ystore[polcount]*ystore[polcount]);

		input>>xqstore[polcount]>>yqstore[polcount]>>zqstore[polcount]>>vqx[polcount]>>vqy[polcount]>>vqz[polcount]>>pol[polcount]>>xustore[polcount]>>yustore[polcount]>>zustore[polcount]>>vux[polcount]>>vuy[polcount]>>vuz[polcount]>>pol[polcount]>>xdstore[polcount]>>ydstore[polcount]>>zdstore[polcount]>>vdx[polcount]>>vdy[polcount]>>vdz[polcount]>>pol[polcount];

		xstore[polcount]=xdstore[polcount];
		ystore[polcount]=ydstore[polcount];
		zstore[polcount]=zdstore[polcount];
		R[polcount]=sqrt(xstore[polcount]*xstore[polcount]+ystore[polcount]*ystore[polcount]);
		
		Rq[polcount]=sqrt(xqstore[polcount]*xqstore[polcount]+yqstore[polcount]*yqstore[polcount]);
		Ru[polcount]=sqrt(xustore[polcount]*xustore[polcount]+yustore[polcount]*yustore[polcount]);
		Rd[polcount]=sqrt(xdstore[polcount]*xdstore[polcount]+ydstore[polcount]*ydstore[polcount]);
		
		vratio[polcount]=vdz[polcount]/sqrt(vdx[polcount]*vdx[polcount]+vdy[polcount]*vdy[polcount]+vdz[polcount]*vdz[polcount]);
		vR[polcount]=abs(vdx[polcount])/sqrt(vdx[polcount]*vdx[polcount]+vdy[polcount]*vdy[polcount]+vdz[polcount]*vdz[polcount]);
		vrho[polcount]=sqrt(vdx[polcount]*vdx[polcount]+vdy[polcount]*vdy[polcount])/vdz[polcount];

		polcount++;
	}

	cout<<"data_read"<<endl;
	
	for(int n=0;n<polcount;n++)
	{
		polfull=polfull+pol[n];
		
		if(Rq[n]<=8.255)
		{
			polq=polq+pol[n];
			qn++;
		}
		
		if(Ru[n]<=8.225)
		{
			polus=polus+pol[n];
			un++;
		}
		
		if(Rd[n]<=8.225)
		{
			polds=polds+pol[n];
			dn++;
		}
		
	}
	
	polfull=polfull/polcount;
	polq=polq/qn;
	polus=polus/un;
	polds=polds/dn;
	cout<<"Full Pol	Quad Pol	US Pol		DS Pol"<<endl<<polfull<<"	"<<polq<<"	"<<polus<<"	"<<polds<<endl;

	xmax=xstore[0];
	xmin=xstore[0];

	for(int p=0;p<polcount;p++)
	{
		if(xstore[p]>xmax)
		{xmax=xstore[p];}
		if(xstore[p]<xmin)
		{xmin=xstore[p];}
	}
	
	cout<<xmin<<" "<<xmax<<endl;
	
	xmin=-40; xmax=40;
	
//	xmin=-4.3; xmax=4.3;

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

			if(ystore[p]>=(0-3.175) && ystore[p]<=3.175)
			{
				center1=(x1+x2)/2;
				if(xstore[p]>=center1-3.175 && xstore[p]<=center1+3.175)
				{
					radius=sqrt((xstore[p]-center1)*(xstore[p]-center1)+ystore[p]*ystore[p]);
					if(radius<=3.175)
					{
						xgram[px]++;
						check++;
					}

				}

			}

			x1+=xblip;
			x2+=xblip;
		}

	}

	cout<<polcount<<" "<<check<<endl;

	for(int p=0;p<xbins;p++)
	{
		xgram[p]=xgram[p]/2.72;
		errxgram[p]=400;
//		cout<<xgram[p]<<" "<<xbox[p]<<endl;
	}



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

		gr0->SetTitle("Simulated Profile DS RGA");
		gr0->GetXaxis()->SetTitle("Off-axis position (mm)");
		gr0->GetYaxis()->SetTitle("Number of Particles");
		gr0->GetXaxis()->SetTitleSize(.05);
		gr0->GetYaxis()->SetTitleSize(.05);
		gr0->GetXaxis()->CenterTitle(1);
		gr0->GetYaxis()->CenterTitle(1);
		gr0->GetXaxis()->SetTitleOffset(.91);
		gr0->GetYaxis()->SetTitleOffset(.9);
		gr0->GetXaxis()->SetLabelSize(.05);
		gr0->GetYaxis()->SetLabelSize(.05);
//		gr0->GetYaxis()->SetRangeUser(0,3e3);
		gr0->Draw("ap");
		gr0->Fit("fitg","r");



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
