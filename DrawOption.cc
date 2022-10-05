{
  // Use times new roman, precision 2 
  Int_t snFont        = 132;
  // Line thickness
  Double_t snWidth    = 2.00;
  // Text size
  Double_t snTSize    = 0.08; 
  
  // use plain black on white colors
  // gROOT->SetStyle("Plain"); 
  TStyle *snStyle= new TStyle("snStyle","Sn plots style");

  snStyle->SetFillStyle(1001);   // solid
  snStyle->SetFrameFillColor(0);
  snStyle->SetFrameBorderMode(0);
  snStyle->SetPadBorderMode(0);
  snStyle->SetPadColor(0);
  snStyle->SetCanvasBorderMode(0);
  snStyle->SetCanvasColor(0);
  snStyle->SetStatColor(0);
  snStyle->SetLegendBorderSize(0);
  
  // If you want the usual gradient palette (blue -> red)
  // https://root.cern.ch/doc/master/classTColor.html
  snStyle->SetPalette(55);
  
  // set the paper & margin sizes
  snStyle->SetPaperSize(20,26);
  snStyle->SetPadTopMargin(0.05);
  snStyle->SetPadRightMargin(0.12); // increase for colz plots
  snStyle->SetPadBottomMargin(0.12);
  snStyle->SetPadLeftMargin(0.12);
  
  // use large fonts
  snStyle->SetTextFont(snFont);
  snStyle->SetTextSize(snTSize);
  snStyle->SetLabelFont(snFont,"x");
  snStyle->SetLabelFont(snFont,"y");
  snStyle->SetLabelFont(snFont,"z");
  snStyle->SetLabelSize(snTSize,"x");
  snStyle->SetLabelSize(snTSize,"y");
  snStyle->SetLabelSize(snTSize,"z");
  snStyle->SetTitleFont(snFont);
  snStyle->SetTitleFont(snFont,"x");
  snStyle->SetTitleFont(snFont,"y");
  snStyle->SetTitleFont(snFont,"z");
  snStyle->SetTitleSize(snTSize,"x");
  snStyle->SetTitleSize(snTSize,"y");
  snStyle->SetTitleSize(snTSize,"z");

  // use medium bold lines and thick markers
  snStyle->SetLineWidth(snWidth);
  snStyle->SetFrameLineWidth(snWidth);
  snStyle->SetHistLineWidth(snWidth);
  snStyle->SetFuncWidth(snWidth);
  snStyle->SetGridWidth(snWidth);
  snStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
  snStyle->SetMarkerStyle(20);
  snStyle->SetMarkerSize(1.0);
  
  // label offsets
  snStyle->SetLabelOffset(0.010,"X");
  snStyle->SetLabelOffset(0.010,"Y");

  // by default, do not display histogram decorations:
  snStyle->SetOptStat(0);
  // show only nent -e , mean - m , rms -r, name -n
  //snStyle->SetOptStat("emrn");
  // full opts at http://root.cern.ch/root/html/TStyle.html#TStyle:SetOptStat
  snStyle->SetStatFormat("6.3g"); // specified as c printf options
  //snStyle->SetOptTitle(0);
  snStyle->SetOptFit(1);
  // order is probability, Chi2, errors, parameter
  // snStyle->SetOptFit(1011); 
  // titles
  snStyle->SetTitleOffset(0.95,"X");
  snStyle->SetTitleOffset(0.95,"Y");
  snStyle->SetTitleOffset(1.2,"Z");
  snStyle->SetTitleFillColor(0);
  snStyle->SetTitleStyle(0);
  snStyle->SetTitleBorderSize(0);
  snStyle->SetTitleFont(snFont,"title");
  snStyle->SetTitleX(0.0);
  snStyle->SetTitleY(1.0); 
  snStyle->SetTitleW(1.0);
  snStyle->SetTitleH(0.05);
  
  // look of the statistics box:
  snStyle->SetStatBorderSize(0);
  snStyle->SetStatFont(snFont);
  snStyle->SetStatFontSize(0.08);
  //snStyle->SetStatX(0.9);
  snStyle->SetStatX(0.99);
  snStyle->SetStatY(0.99);
  snStyle->SetStatW(0.22);
  snStyle->SetStatH(0.2);

  // put tick marks on top and RHS of plots
  // snStyle->SetPadTickX(1);
  // snStyle->SetPadTickY(1);

  // histogram divisions: only 5 in x to avoid label overlaps
  snStyle->SetNdivisions(509,"x");
  snStyle->SetNdivisions(509,"y");

  gROOT->SetStyle("snStyle");
  gROOT->ForceStyle();
  
  // add Sn label
  TPaveText *snName = new TPaveText(gStyle->GetPadLeftMargin() + 0.05,
				    0.87 - gStyle->GetPadTopMargin(),
				    gStyle->GetPadLeftMargin() + 0.20,
				    0.95 - gStyle->GetPadTopMargin(),
				    "BRNDC");
  snName->AddText("Sn");
  snName->SetFillColor(0);
  snName->SetTextAlign(12);
  snName->SetBorderSize(0);

  TText *snLabel = new TText();
  snLabel->SetTextFont(snFont);
  snLabel->SetTextColor(1);
  snLabel->SetTextSize(snTSize);
  snLabel->SetTextAlign(12);

  TLatex *snLatex = new TLatex();
  snLatex->SetTextFont(snFont);
  snLatex->SetTextColor(1);
  snLatex->SetTextSize(snTSize);
  snLatex->SetTextAlign(12);

  /// xji
  snStyle->SetEndErrorSize(4);
}

