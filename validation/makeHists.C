makeHist(double bExp, double deltaB){
	
	TFile *tf = new TFile("data/histograms.root","recreate");
	TH1F *sigData = new TH1F("signalData","",1,0,1);
	TH1F *controlData = new TH1F("controlData","",1,0,1);
	TH1F *bkgExp = new TH1F("bExp","",1,0,1);
	TH1F *tauHist = new TH1F("tau","",1,0,1);
	TH1F *unit = new TH1F("unit","",1,0,1);

	double tau = bExp/deltaB/deltaB;

	sigData->Fill(0.5,bExp);
	controlData->Fill(0.5,bExp*tau);
	bkgExp->Fill(0.5, bExp);
	tauHist->Fill(0.5, tau);

	unit->Fill(0.5,1);

	//hd_neg->Write("histograms.root");
	tf->Write();
	tf->Close();

  }

  void runAll(double bExp=50, double deltaB=3.){
  	makeHist(bExp,deltaB);
 	gSystem->Exec("hist2workspace config/example_DataDriven.xml");
 	gSystem->Load("$ROOTSYS/tutorials/roostats/StandardHypoTestInvDemo.C");
	StandardHypoTestInvDemo("results/example_DataDriven_combined_GaussExample_model.root", "combined", "ModelConfig","","obsData",2,3,true,10,0,2*sqrt(bExp)+2*deltaB);


  }
