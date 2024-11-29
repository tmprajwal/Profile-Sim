#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<algorithm>

#define pi 3.14159

void countour_plot()
{
	double xi[100000000], yi[100000000], zi[100000000], xqstore[10000000], yqstore[10000000], zqstore[10000000], xustore[10000000], yustore[10000000], zustore[10000000], xdstore[10000000], ydstore[10000000], zdstore[10000000];
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
	double radius, center;

	ifstream input;
	input.open("profiles_small_stept_largeN.txt");
//	input.open("profiles_PB_small_step_VLN_theta0p1.txt");
//	input.open("profiles_half_grad.txt");

	while(!input.eof())
	{
//		input>>xstore[polcount]>>ystore[polcount]>>zstore[polcount]>>vx[polcount]>>vy[polcount]>>vz[polcount]>>pol[polcount];

		input>>xqstore[polcount]>>yqstore[polcount]>>zqstore[polcount]>>vqx[polcount]>>vqy[polcount]>>vqz[polcount]>>pol[polcount]>>xustore[polcount]>>yustore[polcount]>>zustore[polcount]>>vux[polcount]>>vuy[polcount]>>vuz[polcount]>>pol[polcount]>>xdstore[polcount]>>ydstore[polcount]>>zdstore[polcount]>>vdx[polcount]>>vdy[polcount]>>vdz[polcount]>>pol[polcount];

//		input>>xqstore[polcount]>>yqstore[polcount]>>zqstore[polcount]>>xustore[polcount]>>yustore[polcount]>>zustore[polcount]>>xdstore[polcount]>>ydstore[polcount]>>zdstore[polcount];
/*
		xstore[polcount]=xqstore[polcount];
		ystore[polcount]=yqstore[polcount];
		zstore[polcount]=zqstore[polcount];
//*/
/*
		xstore[polcount]=xustore[polcount];
		ystore[polcount]=yustore[polcount];
		zstore[polcount]=zustore[polcount];
//*/
//*
		xstore[polcount]=xdstore[polcount];
		ystore[polcount]=ydstore[polcount];
		zstore[polcount]=zdstore[polcount];
		vratio[polcount]=vdz[polcount]/sqrt(vdx[polcount]*vdx[polcount]+vdy[polcount]*vdy[polcount]+vdz[polcount]*vdz[polcount]);
		vR[polcount]=abs(vdx[polcount])/sqrt(vdx[polcount]*vdx[polcount]+vdy[polcount]*vdy[polcount]+vdz[polcount]*vdz[polcount]);
		vrho[polcount]=sqrt(vdx[polcount]*vdx[polcount]+vdy[polcount]*vdy[polcount])/vz;
//*/
		polcount++;
	}


	TCanvas *ABS= new TCanvas("ABS","Predicted Beam Profile",10,10,1280,620);
//	TCanvas *ABS= new TCanvas("ABS","Predicted Beam Profile",10,10,1000,1000);
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

		gr0=new TGraphErrors(polcount, xstore, ystore, 0, 0);

		gr0->SetMarkerStyle(20);
		gr0->SetMarkerColor(1);
		gr0->SetMarkerSize(1.0);
		gr0->SetLineWidth(0.5);
		gr0->SetLineColor(1);

		gr0->SetTitle("Simulated Profile At DS RGA");
//		gr0->GetTitle()->SetTitleSize(.08);
		gr0->GetXaxis()->SetTitle("Off-axis position (mm)");
		gr0->GetYaxis()->SetTitle("Number of Particles");
//		gr0->GetYaxis()->SetTitle("velocity (m/s)");
		gr0->GetXaxis()->SetTitleSize(.05);
		gr0->GetYaxis()->SetTitleSize(.05);
		gr0->GetXaxis()->CenterTitle(1);
		gr0->GetYaxis()->CenterTitle(1);
		gr0->GetXaxis()->SetTitleOffset(.91);
		gr0->GetYaxis()->SetTitleOffset(.9);
		gr0->GetXaxis()->SetLabelSize(.05);
		gr0->GetYaxis()->SetLabelSize(.05);
		gr0->Draw("ap");




	input.close();


}
