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

  const int numTest = 100000;
  for(int i=0; i < numTest; i++){

    TGraph* grInterp = RootTools::interpolateWithStartTime(gr, gr->GetX()[0], 1./2.6, 256);
    delete grInterp;

  }

  return 0;


}
