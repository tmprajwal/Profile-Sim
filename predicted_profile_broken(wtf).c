//This code is for fitting beam profiles saved from the Genya (rewritten) simulation.

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

void predicted_profile()
{
	cout<<"hello"<<endl;
	
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
	double R[100000000];
	double radius, center1;

	ifstream input;
	input.open("profiles_small_stept_largeN.txt");
//	input.open("profiles_PB_small_step_VLN_theta0p1.txt");
//	input.open("profiles_half_grad.txt");
/*
	ifstream in2;
	in2.open("initial_profile_1p5K_1.txt");

	while(!in2.eof())
	{
		in2>>xi[count]>>yi[count]>>zi[count];
		count++;
	}

	xmux=xi[0];
	xmun=xi[0];

	for(int p=0;p<count/100;p++)
	{
		if(xi[p]>xmux)
		{xmux=xi[p];}
		if(xi[p]<xmun)
		{xmun=xi[p];}
	}

	xblup=(xmux-xmun)/xbins;
	xu1=xmun;
	xu2=xmun+xblup;
	cout<<xmun<<" "<<xmux<<endl;

	for(int p=0;p<xbins;p++)
	{
		xslot[p]=(xu1+xu2)/2;
		errxslot[p]=xblup/1000;
		xu1+=xblup;
		xu2+=xblup;
//		cout<<xslot[p]<<endl;
//		cout<<xblip<<" "<<xmin<<" "<<xmax<<" "<<x1<<" "<<x2<<endl;
	}

	for(int p=0;p<count/100;p++)
	{
		xu1=xmun;
		xu2=xmun+xblup;

		for(int px=0;px<xbins;px++)
		{
//			cout<<xu1<<" "<<xu2<<" "<<xi[p]<<endl;

			if(xi[p]>=xu1 && xi[p]<=xu2)
			{
				xhist[px]++;
				check++;
				break;
			}

//			if(px==0){dcheck++;}
			xu1+=xblup;
			xu2+=xblup;
		}

	}

	cout<<count<<" "<<check<<endl;

	for(int p=0;p<xbins;p++)
	{errxhist[p]=sqrt(xhist[p]);}
//*/
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
		R[polcount]=sqrt(xstore[polcount]*xstore[polcount]+ystore[polcount]*ystore[polcount]);
		
		vratio[polcount]=vdz[polcount]/sqrt(vdx[polcount]*vdx[polcount]+vdy[polcount]*vdy[polcount]+vdz[polcount]*vdz[polcount]);
		vR[polcount]=abs(vdx[polcount])/sqrt(vdx[polcount]*vdx[polcount]+vdy[polcount]*vdy[polcount]+vdz[polcount]*vdz[polcount]);
		vrho[polcount]=sqrt(vdx[polcount]*vdx[polcount]+vdy[polcount]*vdy[polcount])/vdz[polcount];
//*/
//		r[polcount]=sqrt(xstore[polcount]*xstore[polcount]+ystore[polcount]*ystore[polcount]);
/*
		if(abs(xstore[polcount])>25)
		{
			xstore[polcount]=0-25+(double(25)+25)*rand()/RAND_MAX;
		}

		if(abs(ystore[polcount])>25)
		{
			ystore[polcount]=0-25+(double(25)+25)*rand()/RAND_MAX;
		}
*/
		polcount++;
	}

//	double rmax=r[0];
//	double rmin=r[0];

	xmax=xstore[0];
	xmin=xstore[0];

	for(int p=0;p<polcount;p++)
	{
		if(xstore[p]>xmax/* && xstore[p]<=20*/)
		{xmax=xstore[p];}
		if(xstore[p]<xmin/* && xstore[p]>=0-20*/)
		{xmin=xstore[p];}
/*
		if(r[p]>rmax)
		{rmax=r[p];}
		if(r[p]<rmin)
		{rmin=r[p];}
*/
	}

	xblip=(xmax-xmin)/xbins;
	xblip=3;
	xbins=(xmax-xmin)/xblip;
	x1=xmin;
	x2=xmin+xblip;
/*
	double rblip=(rmax-rmin)/xbins;
	double r1=rmin;
	double r2=rmin+rblip;
*/
	for(int p=0;p<xbins;p++)
	{
		xbox[p]=(x1+x2)/2;
		errxbox[p]=xblip/1000;
		x1+=xblip;
		x2+=xblip;
//		cout<<xblip<<" "<<xmin<<" "<<xmax<<" "<<x1<<" "<<x2<<endl;
/*
		rbox[p]=(r1+r2)/2;
		errrbox[p]=rblip/1000;
		r1+=rblip;
		r2+=rblip;
*/
	}

	check=0;

	for(int p=0;p<polcount;p++)
	{
		x1=xmin;
		x2=xmin+xblip;
/*
		r1=rmin;
		r2=rmin+rblip;
*/
		for(int px=0;px<xbins;px++)
		{
//*
			if(ystore[p]>=(0-3.175) && ystore[p]<=3.175)
			{
				center1=(x1+x2)/2;
//				if(xstore[p]>=x1 && xstore[p]<=x2)
				if(xstore[p]>=center1-3.175 && xstore[p]<=center1+3.175)
				{
					radius=sqrt((xstore[p]-center1)*(xstore[p]-center1)+ystore[p]*ystore[p]);
//					cout<<radius<<endl;
					if(radius<=3.175)
					{
						xgram[px]++;
						check++;
						//break;
					}

				}

			}
//*/
/*
			if(r[p]>=r1 && r[p]<=r2)
			{
				rgram[px]++;
				break;
			}
*/
//			if(px==0){dcheck++;}
			x1+=xblip;
			x2+=xblip;
/*
			r1+=rblip;
			r2+=rblip;
*/		}

	}

	cout<<polcount<<" "<<check<<endl;
/*
	for(int p=0;p<xbins;p++)
	{errxgram[p]=sqrt(xgram[p]);}
//*/
//*
	for(int p=0;p<xbins;p++)
	{
		xgram[p]=xgram[p]/2.72;
		errxgram[p]=400;
	}
//*/



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
/*
		ABS->cd(1);
		gr1=new TGraphErrors(xbins, xslot, xhist, errxslot, errxhist);
		gr1->SetMarkerStyle(20);
		gr1->SetMarkerColor(1);
		gr1->SetMarkerSize(1.0);
		gr1->SetLineWidth(0.5);
		gr1->SetLineColor(1);

		gr1->SetTitle("Initial Profile");
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
//*/




		ABS->cd(1);
		auto gr0=new TGraphErrors(xbins, xbox, xgram, 0, 0);
//		gr0=new TGraphErrors(xbins, xbox, xgram, errxbox, errxgram);	
//		auto gr0=new TGraphErrors(polcount, R, vdz, 0, 0);	
/*
		fitG=new TF1("fitg","gaus",-20,20);
//		fitG->SetParameter(0,35000);
//		fitG->SetParameter(1,0);
//		fitG->SetParameter(2,2);
//*/

		TF1 *fitG=new TF1("fitg",fitg,-20,20,4);
//		fitG->SetParameter(0,800);
		fitG->SetParameter(0,5000);
		fitG->SetParameter(1,0);
		fitG->SetParameter(2,1);
		fitG->SetParameter(3,0);
		fitG->SetParName(0,"Peak");
		fitG->SetParName(1,"Center");
		fitG->SetParName(2,"Sigma");
		fitG->SetParName(3,"Background");
/*
		TF1 *fitGC=new TF1("fitgc",fitgc,-110,130,3);
		fitGC->SetParameter(0,450);
		fitGC->SetParameter(1,0);
		fitGC->SetParameter(2,2);
		fitGC->SetParameter(3,50);
		fitGC->SetParameter(4,1/28);
		fitGC->SetParameter(5,0);
		fitGC->SetParName(0,"Peak");
		fitGC->SetParName(1,"Center");
		fitGC->SetParName(2,"Sigma");
		fitGC->SetParName(3,"Amplitude");
		fitGC->SetParName(4,"Wavelength");
		fitGC->SetParName(5,"Phase");
//		fitGC->SetParameter(0,400);
//		fitGC->SetParameter(1,1/16);
//		fitGC->SetParameter(2,0);

		TF1 *Cauchy=new TF1("cauchy",cauchy,-7,7,4);
//		Cauchy->SetParameter(0,200);
		Cauchy->SetParameter(0,1800);
		Cauchy->SetParameter(1,0.1);
		Cauchy->SetParameter(2,4.4);
		Cauchy->SetParameter(3,0);
		Cauchy->SetParName(0,"Peak");
		Cauchy->SetParName(1,"Center");
		Cauchy->SetParName(2,"Sigma");
		Cauchy->SetParName(3,"Background");

		TF1 *fitE=new TF1("fite",fitq,-5,0,2);
//		fitE->SetParameter(0,32000);
//		fitE->SetParameter(1,0);
//		fitE->SetParameter(2,0);
		fitE->SetParName(0,"Constant");
		fitE->SetParName(1,"Half");
		fitE->SetParName(2,"Offset");
		fitE->SetParName(3,"Background");

		TF1 *fitpoly=new TF1("poly","pol6",0,7);
		fitpoly->SetParameter(0,5134);
		fitpoly->SetParameter(1,5365);
		fitpoly->SetParameter(2,3412);
		fitpoly->SetParameter(3,1268);
		fitpoly->SetParameter(4,260.5);
		fitpoly->SetParameter(5,27.5);
		fitpoly->SetParameter(6,1.2);

//		overlay=new TF1("rsq","1000/(x*x)",-6,0);
*/
		gr0->SetMarkerStyle(20);
		gr0->SetMarkerColor(1);
		gr0->SetMarkerSize(1.0);
		gr0->SetLineWidth(1);
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
		gr0->Fit("fitg","r");
//		gr0->Fit("cauchy","r");
//		gr0->Fit("fite","r");
//		gr0->Fit("poly","r");
//		overlay->Draw();

		double peak, errpeak, center, errcenter, sigma, errsigma, back, errback;
		peak=fitG->GetParameter(0);
		center=fitG->GetParameter(1);
		sigma=fitG->GetParameter(2);
		back=fitG->GetParameter(3);
		errpeak=fitG->GetParError(0);
		errcenter=fitG->GetParError(1);
		errsigma=fitG->GetParError(2);
		errback=fitG->GetParError(3);

		double integral=fitG->Integral(35,90);
		double errintegral=fitG->IntegralError(35,90);




	cout<<endl<<"Area: "<<integral<<" : "<<errintegral<<" (torr*mm^2)"<<endl<<endl;

	cout<<"peak: "<<peak<<" : "<<errpeak<<" (torr)"<<endl;
	cout<<"center: "<<center<<" : "<<errcenter<<" (mm)"<<endl;
	cout<<"sigma: "<<sigma<<" : "<<errsigma<<" (mm)"<<endl;
	cout<<"background: "<<back<<" : "<<errback<<" (torr)"<<endl;
//*/

	input.close();

}


double fitg(double *xf, double *par){

	double zf=xf[0];

	double gauss=par[0]*exp(0.5*pow(((zf-par[1])/par[2]),2))+par[3];
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









