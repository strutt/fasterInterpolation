#include "TFile.h"
#include "TTree.h"

#include "UsefulAnitaEvent.h"

#include "RootTools.h"

int main(){

  TFile* f = TFile::Open("~/Repositories/Install/share/anitaCalib/fakeEventFile.root");
  TTree* t = (TTree*) f->Get("eventTree");


  UsefulAnitaEvent* usefulEvent = 0;

  t->SetBranchAddress("event", &usefulEvent);


  t->GetEntry(0);


  TGraph* gr = usefulEvent->getGraph(0, AnitaPol::kVertical);


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




  TGraph* grInterp = RootTools::interpolateWithStartTime(gr, gr->GetX()[0], 1./2.6, 256);






  return 0;


}
