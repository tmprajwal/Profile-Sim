//This code is for fitting beam profiles saved from the Genya (rewritten) simulation.

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<algorithm>

#define pi 3.14159

double fitg(double *xf, double *par);

void test()
{
	cout<<"hello"<<endl;
	
	double xi[100000000], yi[100000000], zi[100000000], xqstore[100000000];
/*	double yi[100000000], zi[100000000], xqstore[10000000], yqstore[10000000], zqstore[10000000], xustore[10000000], yustore[10000000], zustore[10000000], xdstore[10000000], ydstore[10000000], zdstore[10000000];
	double vqx[10000000], vqy[10000000], vqz[10000000], vux[10000000], vuy[10000000], vuz[10000000], vdx[10000000], vdy[10000000], vdz[10000000], pol[10000000], r[10000000];
	int polcount=0;

	double vx[10000000],vy[10000000],vz[10000000];

	double xstore[10000000], ystore[10000000], zstore[10000000], vratio[10000000], vR[10000000], vrho[10000000];

	double xbins=12;
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

	ifstream input;
	input.open("profiles_small_stept_largeN.txt");

	while(!input.eof())
	{
		input>>xqstore[polcount]>>yqstore[polcount]>>zqstore[polcount]>>vqx[polcount]>>vqy[polcount]>>vqz[polcount]>>pol[polcount]>>xustore[polcount]>>yustore[polcount]>>zustore[polcount]>>vux[polcount]>>vuy[polcount]>>vuz[polcount]>>pol[polcount]>>xdstore[polcount]>>ydstore[polcount]>>zdstore[polcount]>>vdx[polcount]>>vdy[polcount]>>vdz[polcount]>>pol[polcount];

		xstore[polcount]=xdstore[polcount];
		ystore[polcount]=ydstore[polcount];
		zstore[polcount]=zdstore[polcount];
		vratio[polcount]=vdz[polcount]/sqrt(vdx[polcount]*vdx[polcount]+vdy[polcount]*vdy[polcount]+vdz[polcount]*vdz[polcount]);
		vR[polcount]=abs(vdx[polcount])/sqrt(vdx[polcount]*vdx[polcount]+vdy[polcount]*vdy[polcount]+vdz[polcount]*vdz[polcount]);
		vrho[polcount]=sqrt(vdx[polcount]*vdx[polcount]+vdy[polcount]*vdy[polcount])/vdz[polcount];

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
	xblip=3;
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
		auto gr0=new TGraphErrors(polcount, xstore, ystore, 0, 0);	

		TF1 *fitG=new TF1("fitg",fitg,-20,20,4);
		fitG->SetParameter(0,5000);
		fitG->SetParameter(1,0);
		fitG->SetParameter(2,1);
		fitG->SetParameter(3,0);
		fitG->SetParName(0,"Peak");
		fitG->SetParName(1,"Center");
		fitG->SetParName(2,"Sigma");
		fitG->SetParName(3,"Background");


		gr0->SetMarkerStyle(20);
		gr0->SetMarkerColor(1);
		gr0->SetMarkerSize(1.0);
		gr0->SetLineWidth(0.5);
		gr0->SetLineColor(1);

		gr0->SetTitle("Simulated Profile At DS RGA");
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
		gr0->Draw("ap");
		gr0->Fit("fitg","r");



	input.close();

}


double fitg(double *xf, double *par)
{

	double zf=xf[0];

	double gauss=par[0]*exp(0.5*pow(((zf-par[1])/par[2]),2))+par[3];
	return gauss;
}

