void rootlogon ()
{
  cout << "~/.rootrcdir/rootlogon.C is executing.\n";
  gStyle->SetCanvasColor(kWhite); // background is no longer mouse-dropping white
  gStyle->SetPalette(1,0); // blue to red false color palette. Use 9 for b/w
  gStyle->SetCanvasBorderMode(0); // turn off canvas borders
  gStyle->SetFrameBorderMode(0); // turn off frame borders
  gStyle->SetPadBorderMode(0);
  gStyle->SetPaintTextFormat("5.2f"); // What precision to put numbers if plotted with "TEXT"

  // For publishing:
  gStyle->SetTitleOffset(0.7,"Y");
  gStyle->SetTitleOffset(0.9,"X");
  gStyle->SetOptStat(kFALSE);
  gStyle->SetLineWidth(3);
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.05,"XY");
  gStyle->SetTitleFont(62,"XY");
  gStyle->SetTextSize(1.1);
  gStyle->SetLabelSize(0.03,"XY");
  gStyle->SetLabelFont(62,"XY");
  gStyle->SetPadTopMargin(0.004);
  gStyle->SetPadRightMargin(0.003);
  gStyle->SetPadBottomMargin(0.1);
  gStyle->SetPadLeftMargin(0.075);
  
  gStyle->SetCanvasDefH(830);
  gStyle->SetCanvasDefW(1015);

  gStyle->SetTextFont(62);
  
}
