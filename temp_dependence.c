/*
#include <iostream>
#include <TF1>
#include <TList>
#include <TH1>
#include <TCanvas>
#include <TStyle>
#include <iomanip>
#include <fstream>
#include <TH1F>
#include <TRandom>
#include <TString>
*/
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<algorithm>

double negexp(double *xf, double *par);

double exp(double *xc, double *par);

double simple(double *xs, double *par);

void temp_dependence()
{

	double temp[100], pol[100], trans[100], ratio[100], normratio[100], polq[100], polus[100], polds[100], delpol[100];
	double errpol[100];
	int num=0;
	
	cout<<"initialized"<<endl;
	
	ifstream input;
	input.open("nozzle_temps_1to2K_annulus.txt");
	
	cout<<"File Found"<<endl;
	
	while(!input.eof())
	{
/*		input>>temp[num]>>pol[num]>>errpol[num]>>trans[num]>>ratio[num];
		pol[num]=pol[num]*100;
		errpol[num]=errpol[num]*100;
		trans[num]=trans[num]*100;
		normratio[num]=ratio[num]*100/ratio[0];
//*/		
		input>>temp[num]>>pol[num]>>polq[num]>>polus[num]>>polds[num]>>delpol[num];
		pol[num]=pol[num]*100;
		polq[num]=polq[num]*100;
		polus[num]=polus[num]*100;
		polds[num]=polds[num]*100;
		delpol[num]=delpol[num]*100;
		
		cout<<num<<" "<<temp[num]<<" "<<pol[num]<<" "<<polq[num]<<" "<<polus[num]<<" "<<polds[num]<<" "<<delpol[num]<<endl;
		
		num++;
	}
	
	num--;
	
	cout<<temp[num]<<" "<<temp[num-1]<<endl;




	TCanvas *ABS= new TCanvas("ABS1","Simple Plot",10,10,1280,620);
	ABS->SetGridx(0);
	ABS->SetGridy(0);
	ABS->SetFillColor(10);
	ABS->GetFrame()->SetFillColor(-10);
	ABS->SetBorderMode(0);
	ABS->SetBorderSize(2);
	ABS->Divide(2,2);
//	ABS->Divide(1,2);

	gStyle->SetOptFit();
	gStyle->SetTitleSize(.08,"t");

		ABS->cd(1);
//		auto gr0=new TGraphErrors(num,temp,pol,0,errpol);
		auto gr0=new TGraphErrors(num,temp,pol,0,0);
		gr0->SetMarkerStyle(20);
		gr0->SetMarkerColor(1);
		gr0->SetMarkerSize(1.0);
		gr0->SetLineWidth(1);
		gr0->SetLineColor(1);
		
		auto fit1=new TF1("fit1",negexp,1,5,3);
		fit1->SetParName(0,"Function Coefficient");
		fit1->SetParName(1,"Argument Coefficient");
		fit1->SetParName(2,"Offset");
		fit1->SetLineColor(2);
		
		gr0->SetTitle("Polarization");
//		gr0->GetTitle()->SetTitleSize(.08);
		gr0->GetXaxis()->SetTitle("Nozzle Temperature (K)");
		gr0->GetYaxis()->SetTitle("Polarization (%)");
		gr0->GetXaxis()->SetTitleSize(.05);
		gr0->GetYaxis()->SetTitleSize(.05);
		gr0->GetXaxis()->CenterTitle(1);
		gr0->GetYaxis()->CenterTitle(1);
		gr0->GetXaxis()->SetTitleOffset(.91);
		gr0->GetYaxis()->SetTitleOffset(.9);
		gr0->GetXaxis()->SetLabelSize(.05);
		gr0->GetYaxis()->SetLabelSize(.05);
		gr0->Draw("ap");
		gr0->Fit("fit1","r");
		
		
		ABS->cd(2);
//		auto gr1=new TGraphErrors(num,temp,trans,0,0);
		auto gr1=new TGraphErrors(num,temp,polus,0,0);
		gr1->SetMarkerStyle(20);
		gr1->SetMarkerColor(2);
		gr1->SetMarkerSize(1.0);
		gr1->SetLineWidth(1);
		gr1->SetLineColor(2);
/*		
		auto fit2=new TF1("fit2",exp,1,5,3);
		fit2->SetParName(0,"Function Coefficient");
		fit2->SetParName(1,"Argument Coefficient");
		fit2->SetParName(2,"Offset");
		fit2->SetLineColor(1);
//*/		
		auto fit2=new TF1("fit2",negexp,1,5,3);
		fit2->SetParName(0,"Function Coefficient");
		fit2->SetParName(1,"Argument Coefficient");
		fit2->SetParName(2,"Offset");
		fit2->SetLineColor(1);
		
//		gr1->SetTitle("Transmission");
		gr1->SetTitle("US Annulus");
//		gr1->GetTitle()->SetTitleSize(.08);
		gr1->GetXaxis()->SetTitle("Nozzle Temperature (K)");
//		gr1->GetYaxis()->SetTitle("Transmission (%)");
		gr1->GetYaxis()->SetTitle("Polarization (%)");
		gr1->GetXaxis()->SetTitleSize(.05);
		gr1->GetYaxis()->SetTitleSize(.05);
		gr1->GetXaxis()->CenterTitle(1);
		gr1->GetYaxis()->CenterTitle(1);
		gr1->GetXaxis()->SetTitleOffset(.91);
		gr1->GetYaxis()->SetTitleOffset(.9);
		gr1->GetXaxis()->SetLabelSize(.05);
		gr1->GetYaxis()->SetLabelSize(.05);
		gr1->Draw("ap");
		gr1->Fit("fit2","r");
		
		
		ABS->cd(3);
//		auto gr2=new TGraphErrors(num,temp,ratio,0,0);
		auto gr2=new TGraphErrors(num,temp,polds,0,0);
		gr2->SetMarkerStyle(20);
		gr2->SetMarkerColor(3);
		gr2->SetMarkerSize(1.0);
		gr2->SetLineWidth(1);
		gr2->SetLineColor(3);
		
		fit2->SetLineColor(1);
		
//		gr2->SetTitle("T/R Ratio");
		gr2->SetTitle("DS Annulus");
//		gr2->GetTitle()->SetTitleSize(.08);
		gr2->GetXaxis()->SetTitle("Nozzle Temperature (K)");
//		gr2->GetYaxis()->SetTitle("T/R Ratio (Unitless)");
		gr2->GetYaxis()->SetTitle("Polarization (%)");
		gr2->GetXaxis()->SetTitleSize(.05);
		gr2->GetYaxis()->SetTitleSize(.05);
		gr2->GetXaxis()->CenterTitle(1);
		gr2->GetYaxis()->CenterTitle(1);
		gr2->GetXaxis()->SetTitleOffset(.91);
		gr2->GetYaxis()->SetTitleOffset(.9);
		gr2->GetXaxis()->SetLabelSize(.05);
		gr2->GetYaxis()->SetLabelSize(.05);
		gr2->Draw("ap");
		gr2->Fit("fit2","r");
		
		
		ABS->cd(4);
//		auto gr3=new TGraphErrors(num,temp,normratio,0,0);
		auto gr3=new TGraphErrors(num,temp,delpol,0,0);
		gr3->SetMarkerStyle(20);
		gr3->SetMarkerColor(4);
		gr3->SetMarkerSize(1.0);
		gr3->SetLineWidth(1);
		gr3->SetLineColor(4);
/*		
		fit2->SetParameter(0,146);
		fit2->SetParameter(1,-0.791);
		fit2->SetParameter(2,31.2);
		fit2->SetLineColor(1);
//*/		
		auto fit4=new TF1("fit4","pol1",1,2);
		
//		gr3->SetTitle("Normed T/R Ratio");
		gr3->SetTitle("Polarization Difference");
//		gr3->GetTitle()->SetTitleSize(.08);
		gr3->GetXaxis()->SetTitle("Nozzle Temperature (K)");
//		gr3->GetYaxis()->SetTitle("T/R Ratio (Normalized at T=1 K)");
		gr3->GetYaxis()->SetTitle("Difference (%)");
		gr3->GetXaxis()->SetTitleSize(.05);
		gr3->GetYaxis()->SetTitleSize(.05);
		gr3->GetXaxis()->CenterTitle(1);
		gr3->GetYaxis()->CenterTitle(1);
		gr3->GetXaxis()->SetTitleOffset(.91);
		gr3->GetYaxis()->SetTitleOffset(.9);
		gr3->GetXaxis()->SetLabelSize(.05);
		gr3->GetYaxis()->SetLabelSize(.05);
		gr3->Draw("ap");
//		gr3->Fit("fit2","r");
		gr3->Fit("fit4","r");
		
		
//		auto fit1=new TF1("fit1","pol0",0,50);
//		fit1->SetLineColor(2);

/*
		double peak, errpeak, center, errcenter, sigma, errsigma, back, errback;
		peak=fitG->GetParameter(0);
		center=fitG->GetParameter(1);
		sigma=fitG->GetParameter(2);
		back=fitG->GetParameter(3);
		errpeak=fitG->GetParError(0);
		errcenter=fitG->GetParError(1);
		errsigma=fitG->GetParError(2);
		errback=fitG->GetParError(3);
*/
//		double integral=fitG->Integral(35,90);
//		double errintegral=fitG->IntegralError(35,90);

	input.close();

}



double negexp(double *xf, double *par){

	double zf=xf[0];

	double gauss=100-par[0]*exp(0-par[1]/zf)+par[2];
	return gauss;

//	double cosine=7e-9+(par[0]*cos((28.648*zf)-25));
//	return cosine;


}

double exp(double *xc, double *par){

	double zc=xc[0];

	double cosin=par[0]*exp(par[1]*zc)+par[2];
	return cosin;

//	double cosine=7e-9+(par[0]*cos((28.648*zf)-25));
//	return cosine;


}

double simple(double *xs, double *par){

	double zs=xs[0];

	double simplefit=par[3]+par[0]*sqrt(par[1]*par[1]-((zs-par[2])*(zs-par[2])))/par[1];
//	double simplefit=par[0]*(1-0.5*(zs-par[2])*(zs-par[2])/(par[1]*par[1]));
	return simplefit;




}
