#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<algorithm>

void ratio_optimization()
{
	double *radius= new double[1000000];
	double *trans= new double[1000000];
	double *pol= new double[1000000];
	double *polerr= new double[1000000];
	double *ratio= new double[1000000];

	
	int num=0;
	
	ifstream input;
	input.open("nozzle_6mm_ratio_maxing.txt");
	
	while(!input.eof())
	{
		input>>radius[num]>>trans[num]>>pol[num]>>polerr[num]>>ratio[num];
		trans[num]=100*trans[num];
		pol[num]=100*pol[num];
		polerr[num]=100*polerr[num];
		num++;
	}
	
	num=num-2;
	cout<<num<<" "<<radius[num]<<endl;
	
	TCanvas *ABS= new TCanvas("ABS","Nozzle Radius",10,10,1280,620);
	ABS->SetGridx(0);
	ABS->SetGridy(0);
	ABS->SetFillColor(10);
	ABS->GetFrame()->SetFillColor(-10);
	ABS->SetBorderMode(0);
	ABS->SetBorderSize(2);
	ABS->Divide(2,2);
	
	gStyle->SetOptFit();
	gStyle->SetTitleSize(.08,"t");
	
		ABS->cd(1);
		auto gr0=new TGraphErrors(num,radius,trans,0,0);
		gr0->SetMarkerStyle(20);
		gr0->SetMarkerColor(1);
		gr0->SetMarkerSize(1.0);
		gr0->SetLineWidth(1);
		gr0->SetLineColor(1);
		gr0->SetTitle("Transmission");
//		gr0->GetTitle()->SetTitleSize(.08);
		gr0->GetXaxis()->SetTitle("Aperture Radius (mm");
		gr0->GetYaxis()->SetTitle("Transmission (%)");
		gr0->GetXaxis()->SetTitleSize(.05);
		gr0->GetYaxis()->SetTitleSize(.05);
		gr0->GetXaxis()->CenterTitle(1);
		gr0->GetYaxis()->CenterTitle(1);
		gr0->GetXaxis()->SetTitleOffset(.91);
		gr0->GetYaxis()->SetTitleOffset(.9);
		gr0->GetXaxis()->SetLabelSize(.05);
		gr0->GetYaxis()->SetLabelSize(.05);
		gr0->Draw("ap");
		
		ABS->cd(2);
		auto gr1=new TGraphErrors(num,radius,pol,0,polerr);
		gr1->SetMarkerStyle(20);
		gr1->SetMarkerColor(2);
		gr1->SetMarkerSize(1.0);
		gr1->SetLineWidth(1);
		gr1->SetLineColor(2);
		gr1->SetTitle("Polarization");
//		gr1->GetTitle()->SetTitleSize(.08);
		gr1->GetXaxis()->SetTitle("Aperture Radius (mm)");
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
		
		ABS->cd(3);
		auto gr2=new TGraphErrors(num,radius,ratio,0,0);
		gr2->SetMarkerStyle(20);
		gr2->SetMarkerColor(3);
		gr2->SetMarkerSize(1.0);
		gr2->SetLineWidth(1);
		gr2->SetLineColor(1);
		gr2->SetTitle("T/R ratio");
//		gr2->GetTitle()->SetTitleSize(.08);
		gr2->GetXaxis()->SetTitle("Aperture Radius (mm)");
		gr2->GetYaxis()->SetTitle("T/R Ratio (unitless)");
		gr2->GetXaxis()->SetTitleSize(.05);
		gr2->GetYaxis()->SetTitleSize(.05);
		gr2->GetXaxis()->CenterTitle(1);
		gr2->GetYaxis()->CenterTitle(1);
		gr2->GetXaxis()->SetTitleOffset(.91);
		gr2->GetYaxis()->SetTitleOffset(.9);
		gr2->GetXaxis()->SetLabelSize(.05);
		gr2->GetYaxis()->SetLabelSize(.05);
		gr2->Draw("ap");
		
		ABS->cd(4);
		auto gr3=new TGraphErrors(num,trans,ratio,0,0);
		gr3->SetMarkerStyle(20);
		gr3->SetMarkerColor(4);
		gr3->SetMarkerSize(1.0);
		gr3->SetLineWidth(1);
		gr3->SetLineColor(1);
		
		auto fit1=new TF1("fit1","pol1",0,0.5);
		
		gr3->SetTitle("T/R vs Transmission");
//		gr3->GetTitle()->SetTitleSize(.08);
		gr3->GetXaxis()->SetTitle("Transmission (%)");
		gr3->GetYaxis()->SetTitle("T/R Ratio (unitless)");
		gr3->GetXaxis()->SetTitleSize(.05);
		gr3->GetYaxis()->SetTitleSize(.05);
		gr3->GetXaxis()->CenterTitle(1);
		gr3->GetYaxis()->CenterTitle(1);
		gr3->GetXaxis()->SetTitleOffset(.91);
		gr3->GetYaxis()->SetTitleOffset(.9);
		gr3->GetXaxis()->SetLabelSize(.05);
		gr3->GetYaxis()->SetLabelSize(.05);
		gr3->Draw("ap");
		gr3->Fit("fit1","r");
/*		
		ABS->cd(5);
		auto gr4=new TGraphErrors(num,radius,ratio,0,0);
		gr4->SetMarkerStyle(20);
		gr4->SetMarkerColor(9);
		gr4->SetMarkerSize(1.0);
		gr4->SetLineWidth(1);
		gr4->SetLineColor(1);
		gr4->SetTitle("T/R Ratio");
//		gr4->GetTitle()->SetTitleSize(.08);
		gr4->GetXaxis()->SetTitle("Nozzle Radius (mm)");
		gr4->GetYaxis()->SetTitle("T/R Ratio (unitless)");
		gr4->GetXaxis()->SetTitleSize(.05);
		gr4->GetYaxis()->SetTitleSize(.05);
		gr4->GetXaxis()->CenterTitle(1);
		gr4->GetYaxis()->CenterTitle(1);
		gr4->GetXaxis()->SetTitleOffset(.91);
//		gr4->GetYaxis()->SetTitleOffset(.9);
		gr4->GetXaxis()->SetLabelSize(.05);
		gr4->GetYaxis()->SetLabelSize(.05);
		gr4->Draw("ap");
		
		ABS->cd(6);
		auto gr5=new TGraphErrors(num,trans,pol,0,0);
		gr5->SetMarkerStyle(20);
		gr5->SetMarkerColor(6);
		gr5->SetMarkerSize(1.0);
		gr5->SetLineWidth(1);
		gr5->SetLineColor(1);
		gr5->SetTitle("Polarization Vs Transmission");
//		gr5->GetTitle()->SetTitleSize(.08);
		gr5->GetXaxis()->SetTitle("Transmission (%)");
		gr5->GetYaxis()->SetTitle("Polarization (%)");
		gr5->GetXaxis()->SetTitleSize(.05);
		gr5->GetYaxis()->SetTitleSize(.05);
		gr5->GetXaxis()->CenterTitle(1);
		gr5->GetYaxis()->CenterTitle(1);
		gr5->GetXaxis()->SetTitleOffset(.91);
//		gr5->GetYaxis()->SetTitleOffset(.9);
		gr5->GetXaxis()->SetLabelSize(.05);
		gr5->GetYaxis()->SetLabelSize(.05);
		gr5->Draw("ap");
//*/
/*		
		auto gr6=new TGraphErrors(num,trans,ratio,0,0);
		gr6->SetMarkerStyle(20);
		gr6->SetMarkerColor(7);
		gr6->SetMarkerSize(1.0);
		gr6->SetLineWidth(1);
		gr6->SetLineColor(1);
		gr6->SetTitle("T/R Ratio Vs Transmission");
//		gr6->GetTitle()->SetTitleSize(.08);
		gr6->GetXaxis()->SetTitle("Transmission (%)");
		gr6->GetYaxis()->SetTitle("T/R Ratio (unitless)");
		gr6->GetXaxis()->SetTitleSize(.05);
		gr6->GetYaxis()->SetTitleSize(.05);
		gr6->GetXaxis()->CenterTitle(1);
		gr6->GetYaxis()->CenterTitle(1);
		gr6->GetXaxis()->SetTitleOffset(.91);
//		gr6->GetYaxis()->SetTitleOffset(.9);
		gr6->GetXaxis()->SetLabelSize(.05);
		gr6->GetYaxis()->SetLabelSize(.05);
		gr6->Draw("ap");
//*/
/*		
		auto gr7=new TGraphErrors(num,time,line4,0,0);
		gr7->SetMarkerStyle(20);
		gr7->SetMarkerColor(8);
		gr7->SetMarkerSize(1.0);
		gr7->SetLineWidth(1);
		gr7->SetLineColor(1);
//*/
/*				
		auto gr8=new TGraphErrors(num,time,Nozzle,0,0);
		gr8->SetMarkerStyle(20);
		gr8->SetMarkerColorAlpha(6, 0.01);
		gr8->SetMarkerSize(1.0);
		gr8->SetLineWidth(1);
		gr8->SetLineColor(1);
//*/
	
	
	
	
	input.close();
	
}
