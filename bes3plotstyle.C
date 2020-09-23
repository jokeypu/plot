/**********************************************************
 *                                                        *
 *         BES III Plotstyle: format functions            *
 *                                                        *
 *         August 2009, Niklaus Berger                    *
 *         nberger@ihep.ac.cn                             *
 *                                                        *
 *********************************************************/

#include "bes3plotstyle.h"

#include <TH1.h>
#include <TH1F.h>
#include <TF1.h>
#include <TLatex.h>
#include <TList.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TGraph.h>

#include <iostream>

// Format for data points
void FormatData(TH1 * datahist){
  datahist->SetMarkerStyle(20);
  datahist->SetMarkerSize(1);
  datahist->SetLineWidth(2);

  FormatAxis(datahist->GetXaxis());
  FormatAxis(datahist->GetYaxis());
}

// Format for graph data points
void FormatData(TGraph * datahist){
  datahist->SetMarkerStyle(20);
  datahist->SetMarkerSize(1);
  datahist->SetLineWidth(2);
}

void FormatAxis(TAxis * axis){
  axis->SetLabelFont(42);
  axis->SetLabelSize(0.04);
  axis->SetLabelOffset(0.01);
  axis->SetNdivisions(510);
  axis->SetTitleFont(42);
  axis->SetTitleColor(1);
  axis->SetTitleSize(0.045);
//  axis->SetTitleOffset(1.2);
  axis->SetTitleOffset(1.1);
  axis->CenterTitle();
}

void NameAxes(TH1 * datahist, char * xname, char * yname){
  if(xname)
    datahist->GetXaxis()->SetTitle(xname);
  if(yname)
    datahist->GetYaxis()->SetTitle(yname);
}

// Format for main MC (red line)
void FormatMC1(TH1 * mc1hist){
  mc1hist->SetLineColor(2);
  mc1hist->SetLineWidth(2);
}

// Graph Format for main MC (red line)
void FormatMC1(TGraph * mc1hist){
  mc1hist->SetLineColor(2);
  mc1hist->SetLineWidth(2);
}


// Format for second MC or background
// (Blue shaded area)
void FormatMC2(TH1 * mc2hist){
  mc2hist->SetLineColor(4);
  mc2hist->SetFillColor(4);
  mc2hist->SetLineWidth(2);
  mc2hist->SetFillStyle(3001);
}

// Graph Format for second MC or background
// (Blue line)
void FormatMC2(TGraph * mc2hist){
  mc2hist->SetLineColor(4);
  mc2hist->SetLineWidth(2);
}

// Graph Format for third MC or background
// (Blue line)
void FormatMC3(TGraph * mc3hist){
  mc3hist->SetLineColor(6);
  mc3hist->SetLineWidth(2);
}

// Write "BESIII" in the upper right corner
void WriteBes3(){
  TLatex * bes3 = new TLatex(0.94,0.94, "BESIII");
  bes3->SetNDC();
  bes3->SetTextFont(72);
  bes3->SetTextSize(0.1);
  bes3->SetTextAlign(33);
//  bes3->Draw("SAME");
}

// Write "Preliminary" below BESIII -
// to be used together with WriteBes3()
void WritePreliminary(){
  TLatex * prelim = new TLatex(0.94,0.86, "Preliminary");
  prelim->SetNDC();
  prelim->SetTextFont(62);
  prelim->SetTextSize(0.055);
  prelim->SetTextAlign(33);
//   prelim->Draw();
}

// Make a legend; 
// position will have to change depending on the data shape
void MakeLegend(TH1 * datahist,   // Histogram with data
		char * dataname,  // Description of data
		TH1 * mc1hist, // Histogram with first MC
		char * mc1name, // Description of first MC
		TH1 * mc2hist, // Histogram with 2nd MC/BG
		char * mc2name, // Description of second MC/BG
		double xlow,      // Left edge of legend 
		                  //(fraction of canavas width)
		double ylow,       // Bottom edge of legend
		                  //(fraction of canavas height)
		double xhi,       // Right edge of legend 
		                  //(fraction of canavas width)
		double yhi){       // Top edge of legend
		                  //(fraction of canavas height)

  TLegend * leg = new TLegend(xlow, ylow, xhi, yhi);
  if(datahist && dataname)
    leg->AddEntry(datahist, dataname, "LEP");
  if(mc1hist && mc1name)
    leg->AddEntry(mc1hist, mc1name, "L");
  if(mc2hist && mc2name)
    leg->AddEntry(mc2hist, mc2name, "LF");
  
  leg->SetFillColor(0);
  leg->SetTextFont(42);
//   leg->SetMargin(0.15)
  leg->SetFillStyle(0);
  leg->Draw();

}


// Make a legend; 
// position will have to change depending on the data shape
void MakeLegend(TGraph * datahist,   // Graph with data
		char * dataname,  // Description of data
		TGraph * mc1hist, // Graph with first MC
		char * mc1name, // Description of first MC
		TGraph * mc2hist, // Graph with 2nd MC/BG
		char * mc2name, // Description of second MC/BG
		TGraph * mc3hist, // Graph with 3rd MC/BG
		char * mc3name, // Description of third MC/BG

		double xlow,      // Left edge of legend 
		                  //(fraction of canavas width)
		double ylow,       // Bottom edge of legend
		                  //(fraction of canavas height)
		double xhi,       // Right edge of legend 
		                  //(fraction of canavas width)
		double yhi){       // Top edge of legend
		                  //(fraction of canavas height)

  TLegend * leg = new TLegend(xlow, ylow, xhi, yhi);
  if(datahist && dataname)
    leg->AddEntry(datahist, dataname, "LEP");
  if(mc1hist && mc1name)
    leg->AddEntry(mc1hist, mc1name, "L");
  if(mc2hist && mc2name)
    leg->AddEntry(mc2hist, mc2name, "L");
  if(mc3hist && mc3name)
    leg->AddEntry(mc3hist, mc3name, "L");
  
  leg->SetFillColor(0);
  leg->SetTextFont(42);
//   leg->SetMargin(0.15);
  leg->SetFillStyle(0);
  leg->Draw();

}



// Make a legend (version for fit functions
// position will have to change depending on the data shape
void MakeLegend(TH1 * datahist,   // Histogram with data
		char * dataname,  // Description of data
		char ** functionnames, // list of function names
		double xlow,      // Left edge of legend 
		                  //(fraction of canavas width)
		double ylow,       // Bottom edge of legend
		                  //(fraction of canavas height)
		double xhi,       // Right edge of legend 
		                  //(fraction of canavas width)
		double yhi){       // Top edge of legend
		                  //(fraction of canavas height)

  TLegend * leg = new TLegend(xlow, ylow, xhi, yhi);
  if(datahist && dataname)
    leg->AddEntry(datahist, dataname, "LEP");

  TList* list = datahist->GetListOfFunctions();
  unsigned int nfun = list->GetEntries();

  for(unsigned int i =0;  i < nfun; i++){
    TF1* f1 = (TF1*)(list->At(i));
    leg->AddEntry(f1, functionnames[i], "L");
  }
  leg->SetFillColor(0);
  leg->SetTextFont(42);
//   leg->SetMargin(0.15);
  leg->SetFillStyle(0);
  leg->Draw();

}





// Set the general style options
void SetStyle(){
  // No Canvas Border
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasBorderSize(0);
  // White BG
  gStyle->SetCanvasColor(10);
  // Format for axes
  gStyle->SetLabelFont(42,"xyz");
  gStyle->SetLabelSize(0.06,"xyz");
  gStyle->SetLabelOffset(0.01,"xyz");
  gStyle->SetNdivisions(510,"xyz");
  gStyle->SetTitleFont(42,"xyz");
  gStyle->SetTitleColor(1,"xyz");
  gStyle->SetTitleSize(0.07,"xyz");
  gStyle->SetTitleOffset(1.2,"xyz");
  // No pad borders
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadBorderSize(0);
  // White BG
  gStyle->SetPadColor(10);
  // Margins for labels etc.
  gStyle->SetPadLeftMargin(0.17);
  gStyle->SetPadBottomMargin(0.17);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadTopMargin(0.065);
  // No error bars in x direction
  gStyle->SetErrorX(0);
  // Format legend
  gStyle->SetLegendBorderSize(0);
}

// Style options for "final" plots
// (no stat/fit box)
void SetPrelimStyle(){
  gStyle->SetOptDate(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(0);
}

// Style options for internal meetings
// (stat/fit box)
void SetMeetingStyle(){
  gStyle->SetOptDate(0);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(1111);
  gStyle->SetStatBorderSize(1);
  gStyle->SetStatColor(10);
  gStyle->SetStatFont(42);
  gStyle->SetStatFontSize(0.03);
  gStyle->SetOptFit(1111);
}


// Plot a data MC plot
void PlotDataMC(char * filename,  // Name for the output files, 
		                  // without extension 
		TH1 * datahist,   // Histogram with data
		char * dataname,  // Description of data
		TH1 * mc1hist, // Histogram with first MC
		char * mc1name, // Description of first MC
		TH1 * mc2hist, // Histogram with 2nd MC/BG
		char * mc2name, // Description of second MC/BG
		int prelim,    // Use 1 for Preliminary plot
		                   // 2 for a publication plot
                                   // and 0 for a meeting plot with 
		                   // stat and fit box
		double xlow, // Left edge of legend 
		                    //(fraction of canavas width)
		double ylow,  // Bottom edge of legend
		                    //(fraction of canavas height)
		double xhi,  // Right edge of legend 
		                    //(fraction of canavas width)
		double yhi){  // Top edge of legend
		                    //(fraction of canavas height)
      
  SetStyle();
  if(prelim)
    SetPrelimStyle();
  else
    SetMeetingStyle();

  TCanvas * c1 = new TCanvas("bes3plots","BESIII Plots", 800,600);

  FormatData(datahist);
  if(mc1hist)
    FormatMC1(mc1hist);
  if(mc2hist)
    FormatMC2(mc2hist);
  
  
  datahist->Draw("axis");
  if(mc2hist)
    mc2hist->Draw("same");
  if(mc1hist)
    mc1hist->Draw("same");
  datahist->Draw("Esame");
  datahist->Draw("axissame");
  if(prelim){
    WriteBes3();
    if(prelim == 1)
      WritePreliminary();
  }
 MakeLegend(datahist, dataname,
	    mc1hist,  mc1name,
	    mc2hist,  mc2name);


 char filenameall[256];
 sprintf(filenameall,"%s.eps", filename);
 c1->SaveAs(filenameall);
 sprintf(filenameall,"%s.png", filename);
 c1->SaveAs(filenameall);
}


// Plot data with one or more (fit) functions
// Functions should be part of the data histograms list of functions
// (i.e. perform fits with the "+" option or add other functions via
// datahist->GetListOfFunctions->Add(TF1 * function))
// functionnames should have at least as many elements as the function
// list
void PlotDataFit(char * filename,  // Name for the output files, 
		                  // without extension 
		TH1F * datahist,   // Histogram with data
		char * dataname,  // Description of data
		char ** functionnames,// Names of associated functions
		int prelim,    // Use 1 for Preliminary plot
		                   // 2 for a publication plot
                                   // and 0 for a meeting plot with 
		                   // stat and fit box
		double xlow, // Left edge of legend 
		                    //(fraction of canavas width)
		double ylow,  // Bottom edge of legend
		                    //(fraction of canavas height)
		double xhi,  // Right edge of legend 
		                    //(fraction of canavas width)
		double yhi){  // Top edge of legend
		                    //(fraction of canavas height)


 SetStyle();
 if(prelim)
   SetPrelimStyle();
 else
   SetMeetingStyle();

  TCanvas * c1 = new TCanvas("bes3plots","BESIII Plots", 800,600);

  FormatData(datahist);

  int linestyles[] = {1,2,3,7,9,10};
  //int linecolors[]   = {2,4,kGreen+2,kOrange+7,kMagenta,2};
  int linecolors[]   = {2,4,6, 9,8,2};

 TList* list = datahist->GetListOfFunctions();
  TH1F * datacopy = new TH1F(*datahist); 
  datacopy->Draw("axis");


  unsigned int nfun = list->GetEntries();
  
  if(nfun > 6){
    std::cout << "ERROR: More than six associated functions not forseen" << std::endl;
    return;
  }


  for(unsigned int i =0;  i < nfun; i++){
    TF1* f1 = (TF1*)(list->At(i));
    f1->SetLineColor(linecolors[i]);
    f1->SetLineStyle(linestyles[i]);
    f1->Draw("same");
  }

  MakeLegend(datahist, dataname, functionnames,xlow, ylow, xhi, yhi);

  datacopy->Draw("Esame");
  datacopy->Draw("axissame");

  if(prelim){
    WriteBes3();
    if(prelim==1)
      WritePreliminary();
  }

  char filenameall[256];
  sprintf(filenameall,"%s.eps", filename);
  c1->SaveAs(filenameall);
  sprintf(filenameall,"%s.png", filename);
  c1->SaveAs(filenameall);


}


// Scatter plot
void PlotScatter(char * filename,  // Name for the output files, 
		                   // without extension 
		 TH1 * datahist,   // Histogram with data
		 int prelim       // preliminary plot
		 ){
		                  
      
  SetStyle();
  if(prelim)
    SetPrelimStyle();
  else
    SetMeetingStyle();

  TCanvas * c1 = new TCanvas("bes3plots","BESIII Plots", 800,600);

  FormatData(datahist);
  
  if(datahist->Integral() > 5000)
    datahist->SetMarkerStyle(1);
  else if(datahist->Integral() > 500)
    datahist->SetMarkerSize(0.5);
  
  
  datahist->Draw("");

  if(prelim){
    WriteBes3();
    if(prelim==1)
      WritePreliminary();
  }


 char filenameall[256];
 sprintf(filenameall,"%s.eps", filename);
 c1->SaveAs(filenameall);
 sprintf(filenameall,"%s.png", filename);
 c1->SaveAs(filenameall);
}
