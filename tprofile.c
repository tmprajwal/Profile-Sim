#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<algorithm>

#define pi 3.14159

double fitg(double *xf, double *par);

void tprofile()
{
	
//	cout<<"hello"<<endl;
	double x[100], y[100], integral[100], func[100], errint[100];
	double dum[1000];
	int num=0;
	int count=0;
	int space=0;

	ifstream input;
//	input.open("dum1_expB_1B.csv");
//	input.open("dum1_vrovzB_simlargeN.csv");
//	input.open("dum1_gaussB_simlargeN.csv");
//	input.open("dum1_simulation_pure_circ.csv");
//	input.open("dum1_simresults_3mm.csv");
//	input.open("dum1_simx0O980mm_m5.csv");
	input.open("dum1_sim_vrhoOvz.csv");
//	input.open("dum1_simcirc_1mmbins.csv");


	space=3;
	while(!input.eof())
	{
		if(count%space==0)
		{
			input>>x[num]>>y[num]>>integral[num]>>func[num];

	//		integral[num]=100000*integral[num];
	//		errint[num]=sqrt(integral[num])*2;
			errint[num]=200;
	//		errint[num]=10;
			num++;
		}

		if(count%space!=0)
		{
			input>>dum[num]>>dum[num]>>dum[num]>>dum[num];
		}

		count++;
	}

	num--;




	TCanvas *ABS= new TCanvas("ABS","Predicted Beam Profile",10,10,1280,620);
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

		ABS->cd(1);
		auto gr1=new TGraphErrors(num, x, integral, 0, errint);
		gr1->SetMarkerStyle(20);
		gr1->SetMarkerColor(1);
		gr1->SetMarkerSize(1.0);
		gr1->SetLineWidth(0.5);
		gr1->SetLineColor(1);

		TF1 *fitG=new TF1("fitg",fitg,-20,20,4);
//		fitG->SetParameter(0,1e4);
		fitG->SetParameter(0,450);
		fitG->SetParameter(1,0);
		fitG->SetParameter(2,3);
		fitG->SetParameter(3,0);
		fitG->SetParName(0,"Peak");
		fitG->SetParName(1,"Center");
		fitG->SetParName(2,"Sigma");
		fitG->SetParName(3,"Background");

		gr1->SetTitle("Downstream RGA");
//		gr1->GetTitle()->SetTitleSize(.08);
		gr1->GetXaxis()->SetTitle("radial position (mm)");
		gr1->GetYaxis()->SetTitle("Number of Particles");
		gr1->GetXaxis()->SetTitleSize(.05);
		gr1->GetYaxis()->SetTitleSize(.05);
		gr1->GetXaxis()->CenterTitle(1);
		gr1->GetYaxis()->CenterTitle(1);
		gr1->GetXaxis()->SetTitleOffset(.91);
		gr1->GetYaxis()->SetTitleOffset(.9);
		gr1->GetXaxis()->SetLabelSize(.05);
		gr1->GetYaxis()->SetLabelSize(.05);
		gr1->Draw("ap");
		gr1->Fit("fitg","r");





	input.close();



}

double fitg(double *xf, double *par){

	double zf=xf[0];

	double gauss=par[0]*exp(-0.5*pow(((zf-par[1])/par[2]),2))+par[3];
	return gauss;
}
