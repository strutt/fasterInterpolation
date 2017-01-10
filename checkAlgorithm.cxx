#include "TFile.h"
#include "TTree.h"
#include "TApplication.h"

#include "UsefulAnitaEvent.h"

#include "RootTools.h"

int main(int argc, char* argv[]){

  TFile* f = TFile::Open("~/Repositories/Install/share/anitaCalib/fakeEventFile.root");
  TTree* t = (TTree*) f->Get("eventTree");


  UsefulAnitaEvent* usefulEvent = 0;

  t->SetBranchAddress("event", &usefulEvent);


  t->GetEntry(0);


  TGraph* gr = usefulEvent->getGraph(0, AnitaPol::kVertical);


  const double deltaT = 1.0/2.6;
  const int n = gr->GetN();
  const int gradOffset = 2;

  // there should be 2 special points at the start, 2 special points at the end, and n-1 gradient points
  // there are n+3 points, with indices from 0 to n+2
  // 0, 1 are special
  // 2 -> n are the calculated gradients (i.e. n-1 points for the gradients)
  // n+1, n+2 are special
  std::vector<double> gradient(n + 3, NAN);

  // do the regular gradient points first
  for(int i=0; i < gr->GetN()-1; i++){
    int j = i + gradOffset;
    gradient.at(j) = (gr->GetY()[i+1] - gr->GetY()[i])/(gr->GetX()[i+1] - gr->GetX()[i]);
  }

  // now do the edges...
  // 145   /* non-periodic boundary conditions */
  // 146   m[-2] = 3.0 * m[0] - 2.0 * m[1];
  // 147   m[-1] = 2.0 * m[0] - m[1];
  // 148   m[size - 1] = 2.0 * m[size - 2] - m[size - 3];
  // 149   m[size] = 3.0 * m[size - 2] - 2.0 * m[size - 3];

  // 0, 1 are blank, fill here
  gradient.at(0) = 3.0 * gradient.at(2) - 2.0 * gradient.at(3);
  gradient.at(1) = 2.0 * gradient.at(2) - 1.0 * gradient.at(3);

  // then fill in the back end, like I do to your mum
  gradient.at(n + 1) = 2.0 * gradient.at(n) - 1.0 * gradient.at(n-1);
  gradient.at(n + 2) = 3.0 * gradient.at(n) - 2.0 * gradient.at(n-1);

  // for(size_t j=0; j < gradient.size(); j++){
  //   int i = j - gradOffset;
  //   double y = i < 0 || i >= n ? -999.9999 : gr->GetY()[i];
  //   std::cout << y << "\t" << gradient.at(j) << std::endl;
  // }

  std::vector<double> b(n-1, 0);
  std::vector<double> c(n-1, 0);
  std::vector<double> d(n-1, 0);

  for (int i = 0; i < (n - 1); i++){
    const double NE = fabs (gradient[i + 1 + 2] - gradient[i + 2]) + fabs (gradient[i - 1 + 2] - gradient[i - 2 + 2]);
    if (NE == 0.0){
	b[i] = gradient[i + 2];
	c[i] = 0.0;
	d[i] = 0.0;
	std::cout << "condition 0:\t" << b[i] << "\t" << c[i] << "\t" << d[i] << std::endl;
    }
    else{
      const double h_i = gr->GetX()[i + 1] - gr->GetX()[i];
      const double NE_next = fabs (gradient[2+i + 2] - gradient[2+i + 1]) + fabs (gradient[2+i] - gradient[2+i - 1]);
      const double alpha_i = fabs (gradient[2+i - 1] - gradient[2+i - 2]) / NE;
      double tL_ip1;
      if (NE_next == 0.0){
	tL_ip1 = gradient[2+i];
      }
      else{
	double alpha_ip1 = fabs (gradient[2+i] - gradient[2+i - 1]) / NE_next;
	tL_ip1 = (1.0 - alpha_ip1) * gradient[2+i] + alpha_ip1 * gradient[2+i + 1];
      }
      b[i] = (1.0 - alpha_i) * gradient[2+i - 1] + alpha_i * gradient[2+i];
      c[i] = (3.0 * gradient[2+i] - 2.0 * b[i] - tL_ip1) / h_i;
      d[i] = (b[i] + tL_ip1 - 2.0 * gradient[2+i]) / (h_i * h_i);

      std::cout << "condition 1:\t" << b[i] << "\t" << c[i] << "\t" << d[i] << std::endl;
    }
  }

  const int newN = 256;
  std::vector<double> newYs(newN, 0);
  std::vector<double> newTs(newN, 0);


  // Double_t thisStartTime = gr->GetX()[0];
  // Double_t lastTime = gr->GetX()[gr->GetN()-1];
  // Quantizes the start and end times so data poInt_ts lie at Int_teger multiples of nominal sampling
  // startTime = correlationDeltaT*TMath::Nint(startTime/correlationDeltaT + 0.5);
  // lastTime = correlationDeltaT*TMath::Nint(lastTime/correlationDeltaT - 0.5);
  Double_t startTime = deltaT*TMath::Nint(gr->GetX()[0]/deltaT + 0.5);
  // lastTime = deltaT*TMath::Nint(lastTime/deltaT - 0.5);
  double newT = startTime;

  int j=0;
  for(int i=0; i < newN; i++){

    while(newT > gr->GetX()[j+1]){
      j++;
    }

    std::cout << i << "\t" << j << "\t" << newT << "\t" << gr->GetX()[j]  << "\t" << gr->GetX()[j+1]
	      << "\t" << (newT >= gr->GetX()[j] && newT < gr->GetX()[j+1]) << std::endl;

    double newY = 0;
    if(j < n){
      double delx = newT - gr->GetX()[j];
      newY = gr->GetY()[j] + delx * (b[j] + delx * (c[j] + d[j] * delx));

      // std::cout << i << "\t" << j << delx << "\t" << gr->GetY()[j] << "\t"
      // 		<<  b[j]  << "\t" << c[j] << "\t" << d[j] << std::endl;
    }

    newYs.at(i) = newY;
    newTs.at(i) = newT;

    newT += deltaT;
  }

  TGraph* grMine = new TGraph(newN, &newTs[0], &newYs[0]);


  // TGraph* grInterp = RootTools::interpolateWithStartTime(gr, grMine->GetX()[0], deltaT, newN);
  TGraph* grInterp = RootTools::interpolateWithStartTime(gr, gr->GetX()[0], deltaT, newN);


  TApplication* theApp = new TApplication("App", &argc, argv);
  auto c1 = new TCanvas();
  grMine->Draw("alp");
  grMine->SetLineColor(kRed);
  grInterp->Draw("lsame");
  c1->Update();
  theApp->Run();

  return 0;


}
