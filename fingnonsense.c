#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<algorithm>

void fingnonsense()
{
//out2<<x<<" "<<y<<" "<<z<<" "<<vx<<" "<<vy<<" "<<vz<<endl;
	double *x= new double[1000000];
	double *y= new double[1000000];
	double *z= new double[1000000];
	double *vx= new double[1000000];
	double *vy= new double[1000000];
	double *vz= new double[1000000];
	double *length= new double[100000];
	double *npol= new double[100000];
	
	double **xish= new double*[10000];
	double **yish= new double*[10000];
	double **zish= new double*[10000];
	
	for(int i=0;i<10000;++i)
	{
		xish[i]= new double[1000];
		yish[i]= new double[1000];
		zish[i]= new double[1000];
	}
	
	int num=0;
	int n=0;
	int spot=0;
	int check=0;

	cout<<"hello"<<endl;

	ifstream input;
	input.open("position_1000parts_4.txt");
//	input.open("position_test.txt");
	
	cout<<"goodbye"<<endl;

	while(!input.eof())
	{
		input>>x[num]>>y[num]>>z[num]>>vx[num]>>vy[num]>>vz[num]>>npol[num];

		if(num>0)
		{
//			if(z[num]==0 && npol[num]==-1)
			if(z[num]==0)
			{
				length[n]=num-spot;
				n++;
				spot=num;
				
//				cout<<n<<"	"<<spot<<endl;
			}

		}

//		if(npol[num]==-1)
//		{
			xish[n][num-spot]=x[num];
			yish[n][num-spot]=y[num];
			zish[n][num-spot]=z[num];
//		}

//		if(npol[num]==1)
//		{check++;}

		num++;
	}

	n--;
	num--;

//	cout<<check<<endl;
	cout<<"data acquired"<<endl;


	auto gr0=new TMultiGraph();
	double gx[1000], gy[1000], gz[1000];
	
	cout<<"multigraph initialized"<<endl;
	
	auto gr1=new TGraphErrors(num,x,y,0,0);
	cout<<"gr1 initialized"<<endl;	

	for(int i=3;i<n;i++)
	{
		for(int j=0;j<length[i]-1;j++)
		{
			gx[j]=xish[i][j];
			gy[j]=yish[i][j];
			gz[j]=zish[i][j];
		}

		gr1=new TGraphErrors(length[i]-1,gz,gx,0,0);
		gr1->SetMarkerStyle(20);
		gr1->SetMarkerColor(1);
		gr1->SetMarkerSize(2.0);
		gr1->SetLineWidth(2);
		gr1->SetLineColor(1);

		gr0->Add(gr1);

		for(int j=0;j<length[i]-1;j++)
		{
			gx[j]=0;
			gy[j]=0;
			gz[j]=0;
		}
		
	}
	
	cout<<"multigraph complete"<<endl;




	TCanvas *ABS= new TCanvas("ABS","Particle Tracker",10,10,1280,620);
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
//		gr0=new TGraphErrors(num, z, x, 0, 0);

		TF1 *fit1=new TF1("fit1","pol1",0,187);

		TF1 *fit2=new TF1("fit2","pol2",188,1000);

		TF1 *fits=new TF1("fits","[0]*sin([1]*x+[2])",188,1600);
		fits->SetParameter(0,4);
		fits->SetParameter(1,-0.007507/2);
		fits->SetParameter(2,402.7);
		fits->SetParName(0,"Amplitude");
		fits->SetParName(1,"Wavelength");
		fits->SetParName(2,"Phase");
/*
		gr0->SetMarkerStyle(20);
		gr0->SetMarkerColor(1);
		gr0->SetMarkerSize(1.0);
		gr0->SetLineWidth(0.5);
		gr0->SetLineColor(1);
*/
		gr0->SetTitle("Wanted Particles;Direction of Quad Axis (mm); Radial displacement from Quad Axis (mm)");
//		gr0->GetTitle()->SetTitleSize(.08);
/*
		gr0->GetXaxis()->SetTitle("on axis (mm)");
		gr0->GetYaxis()->SetTitle("off axis (mm)");
		gr0->GetXaxis()->SetTitleSize(.05);
		gr0->GetYaxis()->SetTitleSize(.05);
		gr0->GetXaxis()->CenterTitle(1);
		gr0->GetYaxis()->CenterTitle(1);
		gr0->GetXaxis()->SetTitleOffset(.91);
		gr0->GetYaxis()->SetTitleOffset(.9);
		gr0->GetXaxis()->SetLabelSize(.05);
		gr0->GetYaxis()->SetLabelSize(.05);
//*/
		gr0->Draw("ac");
//		gr0->Fit("fits","r");



/*
		ABS->cd(2);
		gr1=new TGraphErrors(num, x, y, 0, 0);

		fit1=new TF1("fit1","pol1",0,187);

		fit2=new TF1("fit2","pol2",188,1000);

		fits=new TF1("fits","[0]*sin([1]*x+[2])",188,1600);
		fits->SetParameter(0,4);
		fits->SetParameter(1,-0.007507/2);
		fits->SetParameter(2,402.7);
		fits->SetParName(0,"Amplitude");
		fits->SetParName(1,"Wavelength");
		fits->SetParName(2,"Phase");

		gr1->SetMarkerStyle(20);
		gr1->SetMarkerColor(1);
		gr1->SetMarkerSize(1.0);
		gr1->SetLineWidth(0.5);
		gr1->SetLineColor(1);

		gr1->SetTitle("Particle Position");
//		gr1->GetTitle()->SetTitleSize(.08);
		gr1->GetXaxis()->SetTitle("on axis (mm)");
		gr1->GetYaxis()->SetTitle("off axis (mm)");
		gr1->GetXaxis()->SetTitleSize(.05);
		gr1->GetYaxis()->SetTitleSize(.05);
		gr1->GetXaxis()->CenterTitle(1);
		gr1->GetYaxis()->CenterTitle(1);
		gr1->GetXaxis()->SetTitleOffset(.91);
		gr1->GetYaxis()->SetTitleOffset(.9);
		gr1->GetXaxis()->SetLabelSize(.05);
		gr1->GetYaxis()->SetLabelSize(.05);
		gr1->Draw("ap");
		gr1->Fit("fits","r");
//*/

		input.close();




}





