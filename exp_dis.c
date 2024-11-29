#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<algorithm>

#define ntot 10000000

double randgen(double min, double max)	/* generates random number between min and max */
{
double res;	//Variable just for use in this function

	res=min+(double(max)-min)*rand()/RAND_MAX;
//	Uses predefined rand() function and picks a random number between min and max. "res" is set to this value. "RAND_MAX" is the maximum value of the rand() function. MB
	return res;	//function returns the value "res" to the main program. MB
}

double ranz(double z0)
{
	double z, comp;
//	double A=0-6.42e-3;
	double A=0-6.6e-3;
	double B=0.56;
	double cmax=exp(B);
	do
	{
		z=randgen(0.0,z0);
		comp=randgen(0.0,cmax);
	}while(comp>(exp(A*z+B)));
	
	return z;
}

void exp_dis()
{
	double *z= new double[ntot];
	double zgram[10000]={};
	double zbox[10000], errzbox[10000];
	
	double zmin=0;
	double zmax=1295/2;
	double nbins=100;
	double zblip, z1, z2, check;
	
	for(int i=0;i<ntot;i++)
	{
		z[i]=ranz(zmax);
	}
	
	zblip=(zmax-zmin)/nbins;
	z1=zmin;
	z2=zmin+zblip;
	
	for(int p=0;p<nbins;p++)
	{
		zbox[p]=(z1+z2)/2;
		errzbox[p]=zblip/1000;
		z1+=zblip;
		z2+=zblip;
	}
	
	check=0;

	for(int p=0;p<ntot;p++)
	{
		z1=zmin;
		z2=zmin+zblip;

		for(int px=0;px<nbins;px++)
		{

			if(z[p]>=z1 && z[p]<=z2)
			{
				zgram[px]++;
				check++;
			}

			z1+=zblip;
			z2+=zblip;
		}

	}

	cout<<ntot<<" "<<check<<endl;
//*	
	double grammax=0;
	
	for(int p=0;p<nbins;p++)
	{
		if(zgram[p]>grammax)
		{grammax=zgram[p];}	
	}

	for(int p=0;p<nbins;p++)
	{
		zgram[p]=zgram[p]/grammax;
//		errzgram[p]=0.0001;
//		cout<<xgram[p]<<" "<<xbox[p]<<endl;
	}
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
		auto gr0=new TGraphErrors(nbins, zbox, zgram, 0, 0);	
		gr0->SetMarkerStyle(20);
		gr0->SetMarkerColor(1);
		gr0->SetMarkerSize(1.0);
		gr0->SetLineWidth(1);
		gr0->SetLineColor(1);

		gr0->SetTitle("Simulated MC Plate Profile Near Skimmer");
		gr0->GetXaxis()->SetTitle("z (mm)");
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
	
	
	
	
	
	
	
	
	
}
