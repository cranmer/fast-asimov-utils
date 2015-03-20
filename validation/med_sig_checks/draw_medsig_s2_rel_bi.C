// a macro to read in xy data from a file and make a plot
// First run root, then
//  .L drawXY.C
//  draw()


#include "Riostream.h"
void draw() {

  gROOT->Reset();

// first make a canvas and a 2D histogram for the axes

  TCanvas* canvas = new TCanvas("canvas", "canvas", 10, 10, 500, 500);
  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);  
  canvas->SetFrameBorderMode(0);   // need this to turn off red hist frame!

  gROOT->SetStyle("Plain");
  canvas->UseCurrentStyle();

  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.07);
  gPad->SetBottomMargin(0.17);

  gStyle->SetOptStat(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleSize(0.04);

  gStyle->SetTextFont(42);
  gStyle->SetTextSize(0.04);
  gStyle->SetTitleFont(42, "hxy");    // for histogram and axis title
  gStyle->SetLabelFont(42, "xyz");    // for axis labels (values)

  gStyle->SetTitleOffset(0.8, "h");        // what does this do?
  gStyle->SetTitleX(0.15);
  gStyle->SetTitleY(0.99);

  gROOT->ForceStyle();

  // can make histogram or alternatively use the histograms automatically
  // connected to the TF1 or TGraph objects

  double xMin = 0.1;
  double xMax = 100.;
  // double yMin = 0.1;
  // double yMax = 10.;
  double yMin = 0.;
  double yMax = 8.;
  TH2F* axhist = new TH2F("axhist", "title", 10, xMin, xMax, 10, yMin, yMax);
  axhist->SetTitle("");
  axhist->SetXTitle("b");
  axhist->SetYTitle("med[Z|s]");
  gPad->SetLogx(1);
  gPad->SetLogy(0);

  double u[20];
  double x[20][500];

  // Read in data from file and insert in TTree

  TString fileName;
  // cout << "Enter file name: ";
  // cin >> fileName;
  fileName = "medsig_s2_rel_bi.txt";

  ifstream inFile;
  inFile.open(fileName);
  if (inFile.fail()) { 
    cout << "Couldn't open file!" << endl;
    exit(1); 
  }

  bool readLine = true;
  int lineNum = 0;
  int ncol;
  while ( readLine ){

    TString line;
    stringstream ss;
    line.ReadLine(inFile);
    readLine = inFile.good();

    if ( readLine ) {

      TString firstChar = line(0,1);
      bool useLine = firstChar != "#";

      if ( useLine ){

        int i = 0;
        stringstream ss;
        ss << line;              // put whole line into ss
        TString token;
        bool getToken = true;

        while ( getToken ) {
          ss >> token;           // extracts one token
          if ( token.Length() > 0 ) {
            u[i] = token.Atof();
            i++;
          } 
          else {
            getToken = false;
          }
        }             // getToken
        ncol = i;  

        for (int i=0; i<ncol; i++){
          x[i][lineNum] = u[i];
        }
        lineNum++;

      }               // useLine

    }         // readLine

  }           // readLine

  int n = lineNum;
  inFile.close();

  //  for (int i=0; i<n; i++){
  //  cout << i << "  " << x[0][i] << "  " << x[3][i] << "  " 
  //	 << x[4][i] << endl;
  //  }

  TGraph* tg1 = new TGraph(n, x[0], x[1]);
  TGraph* tg2 = new TGraph(n, x[0], x[2]);
  TGraph* tg3 = new TGraph(n, x[0], x[3]);
  TGraph* tg4 = new TGraph(n, x[0], x[4]);
  TGraph* tg5 = new TGraph(n, x[0], x[5]);
  TGraph* tg6 = new TGraph(n, x[0], x[6]);
  TGraph* tg7 = new TGraph(n, x[0], x[7]);
  TGraph* tg8 = new TGraph(n, x[0], x[8]);
  TGraph* tg9 = new TGraph(n, x[0], x[9]);
  TGraph* tg10 = new TGraph(n, x[0], x[10]);
  TGraph* tg11 = new TGraph(n, x[0], x[11]);
  TGraph* tg12 = new TGraph(n, x[0], x[12]);
  TGraph* tg13 = new TGraph(n, x[0], x[13]);
  TGraph* tg14 = new TGraph(n, x[0], x[14]);
  TGraph* tg15 = new TGraph(n, x[0], x[15]);

  TAxis* xa = axhist->GetXaxis();
  TAxis* ya = axhist->GetYaxis();

  xa->SetTitleOffset(1.2);    //  factor multiplies default offset
  ya->SetTitleOffset(1.1);

  xa->SetLabelOffset(0.005);
  ya->SetLabelOffset(0.005);

  xa->SetTickLength(0.015);  // default  = 0.03
  ya->SetTickLength(0.015);  // default  = 0.03

  xa->SetTitleSize(0.05);
  ya->SetTitleSize(0.05);

  //  gPad->SetLogx(1);
  //  xa->SetLimits(90., 700.);

  xa->SetNdivisions(-5); // negative value should force number of divisions?
  ya->SetNdivisions(-4);

  xa->SetLabelSize(0.05);
  ya->SetLabelSize(0.05);

  // Draw axes and then add stuff

  // kDot=1, kPlus, kStar, kCircle=4, kMultiply=5,
  // kFullDotSmall=6, kFullDotMedium=7, kFullDotLarge=8,
  // kFullCircle=20, kFullSquare=21, kFullTriangleUp=22,
  // kFullTriangleDown=23, kOpenCircle=24, kOpenSquare=25,
  // kOpenTriangleUp=26, kOpenDiamond=27, kOpenCross=28,
  // kFullStar=29, kOpenStar=30

  axhist->Draw();

  tg1->SetLineColor(kRed);
  tg1->SetLineWidth(2);
  tg1->SetLineStyle(2);
  tg1->SetMarkerColor(kRed);
  tg1->SetMarkerSize(0.8);
  tg1->SetMarkerStyle(20);
  tg1->Draw("L,same");              // or P for points

  tg2->SetLineColor(kRed);
  tg2->SetLineWidth(2);
  tg2->SetLineStyle(2);
  tg2->SetMarkerColor(kRed);
  tg2->SetMarkerSize(0.8);
  tg2->SetMarkerStyle(20);
  tg2->Draw("L,same");              // or P for points

  tg3->SetLineColor(kRed);
  tg3->SetLineWidth(2);
  tg3->SetLineStyle(2);
  tg3->SetMarkerColor(kRed);
  tg3->SetMarkerSize(0.8);
  tg3->SetMarkerStyle(20);
  // tg3->Draw("L,same");              // or P for points

  tg4->SetLineColor(kBlue);
  tg4->SetLineWidth(2);
  tg4->SetLineStyle(1);
  tg4->SetMarkerColor(kBlue);
  tg4->SetMarkerSize(0.8);
  tg4->SetMarkerStyle(20);
  tg4->Draw("L,same");              // or P for points

  tg5->SetLineColor(kBlue);
  tg5->SetLineWidth(2);
  tg5->SetLineStyle(1);
  tg5->SetMarkerColor(kBlue);
  tg5->SetMarkerSize(0.8);
  tg5->SetMarkerStyle(20);
  tg5->Draw("L,same");              // or P for points

  tg6->SetLineColor(kBlue);
  tg6->SetLineWidth(2);
  tg6->SetLineStyle(1);
  tg6->SetMarkerColor(kBlue);
  tg6->SetMarkerSize(0.8);
  tg6->SetMarkerStyle(20);
  // tg6->Draw("L,same");              // or P for points

  tg7->SetLineColor(kRed);
  tg7->SetLineWidth(2);
  tg7->SetLineStyle(1);
  tg7->SetMarkerColor(kRed);
  tg7->SetMarkerSize(0.8);
  tg7->SetMarkerStyle(21);
  // tg7->Draw("P,same");              // or P for points

  tg8->SetLineColor(kRed);
  tg8->SetLineWidth(2);
  tg8->SetLineStyle(1);
  tg8->SetMarkerColor(kRed);
  tg8->SetMarkerSize(0.8);
  tg8->SetMarkerStyle(21);
  // tg8->Draw("P,same");              // or P for points

  tg9->SetLineColor(kRed);
  tg9->SetLineWidth(2);
  tg9->SetLineStyle(2);
  tg9->SetMarkerColor(kRed);
  tg9->SetMarkerSize(0.8);
  tg9->SetMarkerStyle(21);
  // tg9->Draw("P,same");              // or P for points

  tg10->SetLineColor(kBlack);
  tg10->SetLineWidth(2);
  tg10->SetLineStyle(2);
  tg10->SetMarkerColor(kBlack);
  tg10->SetMarkerSize(0.8);
  tg10->SetMarkerStyle(20);
  tg10->Draw("P,same");              // or P for points

  tg11->SetLineColor(kBlack);
  tg11->SetLineWidth(2);
  tg11->SetLineStyle(2);
  tg11->SetMarkerColor(kBlack);
  tg11->SetMarkerSize(0.8);
  tg11->SetMarkerStyle(20);
  tg11->Draw("P,same");              // or P for points

  tg12->SetLineColor(kBlack);
  tg12->SetLineWidth(2);
  tg12->SetLineStyle(2);
  tg12->SetMarkerColor(kBlack);
  tg12->SetMarkerSize(0.8);
  tg12->SetMarkerStyle(20);
  // tg12->Draw("P,same");              // or P for points

  tg13->SetLineColor(kBlack);
  tg13->SetLineWidth(2);
  tg13->SetLineStyle(3);
  tg13->SetMarkerColor(kBlack);
  tg13->SetMarkerSize(0.8);
  tg13->SetMarkerStyle(20);
  // tg13->Draw("L,same");              // or P for points

  tg14->SetLineColor(kBlack);
  tg14->SetLineWidth(2);
  tg14->SetLineStyle(3);
  tg14->SetMarkerColor(kBlack);
  tg14->SetMarkerSize(0.8);
  tg14->SetMarkerStyle(20);
  // tg14->Draw("L,same");              // or P for points

  tg15->SetLineColor(kBlack);
  tg15->SetLineWidth(2);
  tg15->SetLineStyle(3);
  tg15->SetMarkerColor(kBlack);
  tg15->SetMarkerSize(0.8);
  tg15->SetMarkerStyle(20);
  // tg15->Draw("L,same");              // or P for points


  TLegend* leg = new TLegend(0.53, 0.38, 0.95, .73); // x1, y1, x2, y2
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(tg1, " s / #sqrt{b + #sigma_{b}^{2}}", "l");
  leg->AddEntry(tg4, "  Z_{A}", "l");
  // leg->AddEntry(tg7, " #sqrt{q0}, MC median", "p");
  leg->AddEntry(tg10, " Monte Carlo", "p");
  // leg->AddEntry(tg13, " Z_{bi}" , "l");
  // leg->AddEntry(tg4, "s = 0.03", "l");
  leg->Draw();

  TLatex* tl = new TLatex();
  tl->SetTextAlign(11);
  tl->SetTextSize(0.05);
  tl->SetTextFont(42);
  tl->SetNDC();
  // tl->DrawLatex(.76, 0.33, "s = 2");  
  // tl->DrawLatex(.76, 0.465, "s = 5");  
  // tl->DrawLatex(.76, 0.6, "s = 10");  
  //   tl->DrawLatex(.76, 0.73, "s = 20");  

  tl->DrawLatex(.5, .77, "#sigma_{b}/b = 0.2, 0.5, 1");
  tl->DrawLatex(.5, .85, "s = 2");

  // Fix idiotic problem with frame

  TLine* tli = new TLine();
  tli->SetLineStyle(1);
  tli->SetLineWidth(1);
  tli->DrawLine(xMin, yMin, xMax, yMin);
  tli->DrawLine(xMax, yMin, xMax, yMax);
  tli->DrawLine(xMin, yMax, xMax, yMax);
  tli->DrawLine(xMin, yMin, xMin, yMax);

  TPostScript psfile("medsig_s2_rel_bi.eps", 113);     // 113 makes eps
  canvas->Draw();
  psfile.Close();
  // canvas->Print("plot.gif", "gif");

}

