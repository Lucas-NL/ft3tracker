#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <TCanvas.h>
#include "TLatex.h"
#include "TROOT.h"
#include <TStyle.h>
#include "TSystem.h"
#include <TPaveText.h>
#include <TLegend.h>


#endif

//_________________________________________________________________________________________________
template <typename H> void exportHisto(const H &histo) {
  // gStyle->SetImageScaling(3.);
  TCanvas *c = new TCanvas;
  c->SetBatch();
  std::string imgpath{"images/"};
  gSystem->MakeDirectory(imgpath.c_str());
  H *h = new H(histo);
  h->Draw();
  gSystem->ProcessEvents();
  c->Update();
  for (std::string type : {".pdf", ".png"})
    c->Print((imgpath + std::string(h->GetName()) + type).c_str());
}

//_________________________________________________________________________________________________
void FitSlicesy(TH2F &histo1, TH2F &histo2, std::string title_1 = "-", std::string title_2 = "Resolution") {
  // ## FitSlicesy
  TH2F *h1 = &histo1;
  TH2F *h2 = &histo2;

  h1->FitSlicesY(0, 0, -1, 1);
  // Create a canvas and divide it
  auto CanvasName = std::string(h1->GetName()) + std::string("FitSlicesY");
  TCanvas *c1 = new TCanvas(CanvasName.c_str(), CanvasName.c_str(), 700, 500);
  c1->Divide(2, 1);
  TPad *leftPad = (TPad *)c1->cd(1);
  ;
  leftPad->Divide(1, 2);

  // Draw 2-d original histogram histo1
  leftPad->cd(1);
  gPad->SetTopMargin(0.05);
  h1->Draw();

  // Draw histo2
  leftPad->cd(2);
  h2->Draw();
  TPad *rightPad = (TPad *)c1->cd(2);
  rightPad->Divide(1, 2);
  rightPad->cd(1);
  TH2F *h1_Fit1 = (TH2F *)gDirectory->Get(
      (std::string(h1->GetName()) + std::string("_1")).c_str());
  h1_Fit1->SetStats(0);
  h1_Fit1->SetTitle(title_1.c_str());
  h1_Fit1->Draw();

  // Show fitted "sigma" for each slice
  rightPad->cd(2);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.05);
  TH2F *h1_Fit2 = (TH2F *)gDirectory->Get(
      (std::string(h1->GetName()) + std::string("_2")).c_str());
  h1_Fit2->SetStats(0);
  h1_Fit2->SetTitle(title_2.c_str());
  h1_Fit2->Draw();

  ////h2_InvPtResolution_0->Write();
  // h1_Fit1->Write();
  // h1_Fit2->Write();
  // c1->Write();
}

//_________________________________________________________________________________________________
template <typename H1, typename H2, typename H3, typename H4>
TCanvas summary_report(H1 &histo1, H2 &histo2, H3 &histo3, H4 &histo4,
                       std::string CanvasName, std::string tlt = "Summary",
                       int h1log = 0, // 0 = linear y scale; 1 = log y scale
                       int h2log = 0, int h3log = 0, int h4log = 0,
                       std::string h1_foot = "-", std::string h2_foot = "-",
                       std::string h3_foot = "-", std::string h4_foot = "-") {
  H1 *h1 = &histo1;
  H2 *h2 = &histo2;
  H3 *h3 = &histo3;
  H4 *h4 = &histo4;

  // h1->FitSlicesY(0,0,-1,1);
  // Create a canvas and divide it
  TCanvas *c1 = new TCanvas(CanvasName.c_str(), CanvasName.c_str(), 1920, 1080);
  c1->UseCurrentStyle();
  // c1->SetCanvasSize(1920,1080);
  // gROOT->SetStyle("Bold");

  TLatex Title = TLatex();
  Title.SetTextSize(0.035);
  Title.SetTextAlign(23);
  Title.DrawLatex(.5, .995, tlt.c_str());
  Title.Draw();

  c1->Divide(2, 1);
  TPad *leftPad = (TPad *)c1->cd(1);
  leftPad->SetPad(0.000, 0.000, 0.5, .96);
  leftPad->Divide(1, 2);

  // Top left
  leftPad->cd(1);

  gPad->SetBottomMargin(0.15);
  // gPad->SetRightMargin(0.15);
  gPad->SetLogy(h1log);

  h1->Draw();
  h1->SetMarkerColor(4);
  h1->SetMarkerStyle(21);
  h1->SetMarkerSize(1);
  TLatex h1_entries = TLatex();
  h1_entries.SetNDC();
  h1_entries.SetTextSize(0.035);
  h1_entries.SetTextAlign(23);
  h1_entries.DrawLatex(.08, .08, h1_foot.c_str());
  h1_entries.Draw();

  // Bottom left
  leftPad->cd(2);

  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.10);
  gPad->SetLogy(h2log);

  h2->Draw();
  h2->SetMarkerColor(4);
  h2->SetMarkerStyle(21);
  h2->SetMarkerSize(1);
  TLatex h2_entries = TLatex();
  h2_entries.SetNDC();
  h2_entries.SetTextSize(0.035);
  h2_entries.SetTextAlign(23);
  h2_entries.DrawLatex(.08, .08, h2_foot.c_str());
  h2_entries.Draw();

  TPad *rightPad = (TPad *)c1->cd(2);
  rightPad->SetPad(0.5, 0, 1, 0.97);
  rightPad->Divide(1, 2);

  // Top right
  rightPad->cd(1);

  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.10);
  gPad->SetLogy(h3log);

  h3->Draw();
  h3->SetMarkerColor(4);
  h3->SetMarkerStyle(21);
  h3->SetMarkerSize(1);
  TLatex h3_entries = TLatex();
  h3_entries.SetNDC();
  h3_entries.SetTextSize(0.035);
  h3_entries.SetTextAlign(23);
  h3_entries.DrawLatex(.08, .08, h3_foot.c_str());
  h3_entries.Draw();

  // Bottom right
  rightPad->cd(2);

  gPad->SetBottomMargin(0.15);
  gPad->SetLogy(h4log);

  h4->Draw();
  h4->SetMarkerColor(4);
  h4->SetMarkerStyle(21);
  h4->SetMarkerSize(1);
  TLatex h4_entries = TLatex();
  h4_entries.SetNDC();
  h4_entries.SetTextSize(0.035);
  h4_entries.SetTextAlign(23);
  h4_entries.DrawLatex(.08, .08, h4_foot.c_str());
  h4_entries.Draw();
  c1->cd();

  c1->Update();
  c1->Write();
  return c1;
}

//_________________________________________________________________________________________________
template <typename H1, typename H2, typename H3, typename H4, typename H5,
          typename H6, typename H7, typename H8, typename H9>
TCanvas
summary_report_3x3(H1 &histo1, H2 &histo2, H3 &histo3, H4 &histo4, H5 &histo5,
                   H6 &histo6, H7 &histo7, H8 &histo8, H9 &histo9,
                   std::string CanvasName, std::string tlt = "Summary",
                   int h1log = 0, // 0 = linear y scale; 1 = log y scale
                   int h2log = 0, int h3log = 0, int h4log = 0, int h5log = 0,
                   int h6log = 0, int h7log = 0, int h8log = 0, int h9log = 0,
                   std::string h1_foot = "-", std::string h2_foot = "-",
                   std::string h3_foot = "-", std::string h4_foot = "-",
                   std::string h5_foot = "-", std::string h6_foot = "-",
                   std::string h7_foot = "-", std::string h8_foot = "-",
                   std::string h9_foot = "-") {

  H1 *h1 = &histo1;
  H2 *h2 = &histo2;
  H3 *h3 = &histo3;
  H4 *h4 = &histo4;
  H5 *h5 = &histo5;
  H6 *h6 = &histo6;
  H7 *h7 = &histo7;
  H8 *h8 = &histo8;
  H9 *h9 = &histo9;

  TCanvas *c1 = new TCanvas(CanvasName.c_str(), CanvasName.c_str(), 1920, 1080);
  c1->UseCurrentStyle();

  TLatex Title = TLatex();
  Title.SetTextSize(0.035);
  Title.SetTextAlign(23);
  Title.DrawLatex(.5, .995, tlt.c_str());
  Title.Draw();

  c1->Divide(1, 2);
  TPad *topPad = (TPad *)c1->cd(1);
  topPad->SetPad(0.000, 0.0, 1, .8);

  TPad *bottomPad = (TPad *)c1->cd(2);
  bottomPad->SetPad(0.000, 0.0, 1, .95);
  bottomPad->Divide(3, 3, 0.005, 0.005);

  // New Pad
  bottomPad->cd(1);
  // gPad->SetTopMargin(0.15);
  gPad->SetBottomMargin(0.15);
  // gPad->SetLeftMargin(0.15);
  // gPad->SetRightMargin(0.15);
  gPad->SetLogy(h1log);

  h1->Draw();
  h1->SetMarkerColor(4);
  h1->SetMarkerStyle(21);
  h1->SetMarkerSize(1);
  TLatex h1_entries = TLatex();
  h1_entries.SetNDC();
  h1_entries.SetTextSize(0.045);
  h1_entries.SetTextAlign(11);
  h1_entries.DrawLatex(.08, .08, h1_foot.c_str());
  h1_entries.Draw();

  // New Pad
  bottomPad->cd(2);
  // pad->SetPad(0.33,0.66,0.66,.96);

  gPad->SetBottomMargin(0.15);
  // gPad->SetRightMargin(0.15);
  gPad->SetLogy(h2log);

  h2->Draw();
  h2->SetMarkerColor(4);
  h2->SetMarkerStyle(21);
  h2->SetMarkerSize(1);
  TLatex h2_entries = TLatex();
  h2_entries.SetNDC();
  h2_entries.SetTextSize(0.045);
  h2_entries.SetTextAlign(11);
  h2_entries.DrawLatex(.08, .08, h2_foot.c_str());
  h2_entries.Draw();

  // New Pad
  bottomPad->cd(3);
  // pad->SetPad(0.66,0.66,1,.96);
  gPad->SetBottomMargin(0.15);
  // gPad->SetRightMargin(0.15);
  gPad->SetLogy(h3log);

  h3->Draw();
  h3->SetMarkerColor(4);
  h3->SetMarkerStyle(21);
  h3->SetMarkerSize(1);
  TLatex h3_entries = TLatex();
  h3_entries.SetNDC();
  h3_entries.SetTextSize(0.045);
  h3_entries.SetTextAlign(11);
  h3_entries.DrawLatex(.08, .08, h3_foot.c_str());
  h3_entries.Draw();

  // New Pad
  bottomPad->cd(4);
  // pad->SetPad(0.000,0.5,0.33,.96);
  gPad->SetBottomMargin(0.15);
  // gPad->SetRightMargin(0.15);
  gPad->SetLogy(h4log);

  h4->Draw();
  h4->SetMarkerColor(4);
  h4->SetMarkerStyle(21);
  h4->SetMarkerSize(1);
  TLatex h4_entries = TLatex();
  h4_entries.SetNDC();
  h4_entries.SetTextSize(0.045);
  h4_entries.SetTextAlign(11);
  h4_entries.DrawLatex(.08, .08, h4_foot.c_str());
  h4_entries.Draw();

  // New Pad
  bottomPad->cd(5);
  gPad->SetBottomMargin(0.15);
  // gPad->SetRightMargin(0.15);
  gPad->SetLogy(h5log);

  h5->Draw();
  h5->SetMarkerColor(4);
  h5->SetMarkerStyle(21);
  h5->SetMarkerSize(1);
  TLatex h5_entries = TLatex();
  h5_entries.SetNDC();
  h5_entries.SetTextSize(0.045);
  h5_entries.SetTextAlign(11);
  h5_entries.DrawLatex(.08, .08, h5_foot.c_str());
  h5_entries.Draw();

  // New Pad
  bottomPad->cd(6);
  // pad->SetPad(0.000,0.5,0.33,.96);

  gPad->SetBottomMargin(0.15);
  // gPad->SetLeftMargin(0.15);
  gPad->SetLogy(h6log);

  h6->Draw();
  h6->SetMarkerColor(4);
  h6->SetMarkerStyle(21);
  h6->SetMarkerSize(1);
  TLatex h6_entries = TLatex();
  h6_entries.SetNDC();
  h6_entries.SetTextSize(0.045);
  h6_entries.SetTextAlign(11);
  h6_entries.DrawLatex(.08, .08, h6_foot.c_str());
  h6_entries.Draw();

  // New Pad
  bottomPad->cd(7);
  // pad->SetPad(0.000,0.5,0.33,.96);

  gPad->SetTopMargin(0.15);
  gPad->SetBottomMargin(0.15);
  // gPad->SetLeftMargin(0.15);
  // gPad->SetRightMargin(0.15);
  gPad->SetLogy(h7log);

  h7->Draw();
  h7->SetMarkerColor(4);
  h7->SetMarkerStyle(21);
  h7->SetMarkerSize(1);
  TLatex h7_entries = TLatex();
  h7_entries.SetNDC();
  h7_entries.SetTextSize(0.045);
  h7_entries.SetTextAlign(11);
  h7_entries.DrawLatex(.08, .08, h7_foot.c_str());
  h7_entries.Draw();

  // New Pad
  bottomPad->cd(8);
  // pad->SetPad(0.000,0.5,0.33,.96);

  gPad->SetTopMargin(0.15);
  gPad->SetBottomMargin(0.15);
  // gPad->SetLeftMargin(0.15);
  // gPad->SetRightMargin(0.15);
  gPad->SetLogy(h8log);

  h8->Draw();
  h8->SetMarkerColor(4);
  h8->SetMarkerStyle(21);
  h8->SetMarkerSize(1);
  TLatex h8_entries = TLatex();
  h8_entries.SetNDC();
  h8_entries.SetTextSize(0.045);
  h8_entries.SetTextAlign(11);
  h8_entries.DrawLatex(.08, .08, h8_foot.c_str());
  h8_entries.Draw();

  // New Pad
  bottomPad->cd(9);
  // pad->SetPad(0.000,0.5,0.33,.96);

  gPad->SetTopMargin(0.15);
  gPad->SetBottomMargin(0.15);
  // gPad->SetLeftMargin(0.15);
  // gPad->SetRightMargin(0.15);
  gPad->SetLogy(h9log);

  h9->Draw();
  h9->SetMarkerColor(4);
  h9->SetMarkerStyle(21);
  h9->SetMarkerSize(1);
  TLatex h9_entries = TLatex();
  h9_entries.SetNDC();
  h9_entries.SetTextSize(0.045);
  h9_entries.SetTextAlign(11);
  h9_entries.DrawLatex(.08, .08, h9_foot.c_str());
  h9_entries.Draw();

  c1->cd();

  c1->Update();
  c1->Write();
  CanvasName += ".png";
  c1->Print(CanvasName.c_str());
  return c1;
}

//_________________________________________________________________________________________________
template <typename H1, typename H2, typename H3, typename H4, typename H5,
          typename H6>
TCanvas summary_report_3x2(H1 &histo1, H2 &histo2, H3 &histo3, H4 &histo4,
                           H5 &histo5, H6 &histo6, std::string CanvasName,
                           std::string tlt = "Summary",
                           int h1log = 0, // 0 = linear y scale; 1 = log y scale
                           int h2log = 0, int h3log = 0, int h4log = 0,
                           int h5log = 0, int h6log = 0,
                           std::string h1_foot = "-", std::string h2_foot = "-",
                           std::string h3_foot = "-", std::string h4_foot = "-",
                           std::string h5_foot = "-",
                           std::string h6_foot = "-") {
  H1 *h1 = &histo1;
  H2 *h2 = &histo2;
  H3 *h3 = &histo3;
  H4 *h4 = &histo4;
  H5 *h5 = &histo5;
  H6 *h6 = &histo6;

  // h1->FitSlicesY(0,0,-1,1);
  // Create a canvas and divide it
  TCanvas *c1 = new TCanvas(CanvasName.c_str(), CanvasName.c_str(), 1920, 1080);
  c1->UseCurrentStyle();
  // c1->SetCanvasSize(1920,1080);
  // gROOT->SetStyle("Bold");

  TLatex Title = TLatex();
  Title.SetTextSize(0.035);
  Title.SetTextAlign(23);
  Title.DrawLatex(.5, .995, tlt.c_str());
  Title.Draw();

  c1->Divide(1, 2);
  TPad *topPad = (TPad *)c1->cd(1);
  topPad->SetPad(0.000, 0.0, 1, .8);

  TPad *bottomPad = (TPad *)c1->cd(2);
  bottomPad->SetPad(0.000, 0.0, 1, .95);
  bottomPad->Divide(3, 2, 0.005, 0.005);

  // New Pad
  bottomPad->cd(1);
  // gPad->SetTopMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gPad->SetRightMargin(0.10);
  gPad->SetLogy(h1log);

  h1->Draw();
  h1->SetMarkerColor(4);
  h1->SetMarkerStyle(21);
  h1->SetMarkerSize(1);
  TLatex h1_entries = TLatex();
  h1_entries.SetNDC();
  h1_entries.SetTextSize(0.035);
  h1_entries.SetTextAlign(23);
  h1_entries.DrawLatex(.08, .08, h1_foot.c_str());
  h1_entries.Draw();

  // New Pad
  bottomPad->cd(2);
  // pad->SetPad(0.33,0.66,0.66,.96);

  gPad->SetBottomMargin(0.15);
  // gPad->SetRightMargin(0.15);
  // gPad->SetLeftMargin(0.15);
  gPad->SetLogy(h2log);

  h2->Draw();
  h2->SetMarkerColor(4);
  h2->SetMarkerStyle(21);
  h2->SetMarkerSize(1);
  TLatex h2_entries = TLatex();
  h2_entries.SetNDC();
  h2_entries.SetTextSize(0.035);
  h2_entries.SetTextAlign(23);
  h2_entries.DrawLatex(.08, .08, h2_foot.c_str());
  h2_entries.Draw();

  // New Pad
  bottomPad->cd(3);
  // pad->SetPad(0.66,0.66,1,.96);
  gPad->SetBottomMargin(0.15);
  // gPad->SetRightMargin(0.15);
  gPad->SetLogy(h3log);

  h3->Draw();
  h3->SetMarkerColor(4);
  h3->SetMarkerStyle(21);
  h3->SetMarkerSize(1);
  TLatex h3_entries = TLatex();
  h3_entries.SetNDC();
  h3_entries.SetTextSize(0.035);
  h3_entries.SetTextAlign(23);
  h3_entries.DrawLatex(.08, .08, h3_foot.c_str());
  h3_entries.Draw();

  // New Pad
  bottomPad->cd(4);
  // pad->SetPad(0.000,0.5,0.33,.96);
  gPad->SetBottomMargin(0.15);
  // gPad->SetRightMargin(0.15);
  gPad->SetLogy(h4log);

  h4->Draw();
  h4->SetMarkerColor(4);
  h4->SetMarkerStyle(21);
  h4->SetMarkerSize(1);
  TLatex h4_entries = TLatex();
  h4_entries.SetNDC();
  h4_entries.SetTextSize(0.035);
  h4_entries.SetTextAlign(23);
  h4_entries.DrawLatex(.08, .08, h4_foot.c_str());
  h4_entries.Draw();

  // New Pad
  bottomPad->cd(5);
  gPad->SetBottomMargin(0.15);
  // gPad->SetRightMargin(0.15);
  gPad->SetLogy(h5log);

  h5->Draw();
  h5->SetMarkerColor(4);
  h5->SetMarkerStyle(21);
  h5->SetMarkerSize(1);
  TLatex h5_entries = TLatex();
  h5_entries.SetNDC();
  h5_entries.SetTextSize(0.035);
  h5_entries.SetTextAlign(23);
  h5_entries.DrawLatex(.08, .08, h5_foot.c_str());
  h5_entries.Draw();

  // New Pad
  bottomPad->cd(6);
  // pad->SetPad(0.000,0.5,0.33,.96);

  gPad->SetBottomMargin(0.15);
  // gPad->SetLeftMargin(0.15);
  gPad->SetLogy(h6log);

  h6->Draw();
  h6->SetMarkerColor(4);
  h6->SetMarkerStyle(21);
  h6->SetMarkerSize(1);
  TLatex h6_entries = TLatex();
  h6_entries.SetNDC();
  h6_entries.SetTextSize(0.035);
  h6_entries.SetTextAlign(23);
  h6_entries.DrawLatex(.08, .08, h6_foot.c_str());
  h6_entries.Draw();

  c1->cd();

  c1->Update();
  c1->Write();
  return c1;
}

//_________________________________________________________________________________________________
void TH3HistoSlicer(string cname, string ctitle, std::unique_ptr<TH3F>& Hist3DDBG, std::vector<double> list, double window, string value, float marker_size=0.7) {

  auto canva = new TCanvas(cname.c_str(), ctitle.c_str(), 1080, 1080);
  int marker = kFullCircle;
  int colormarker = 1;

  Hist3DDBG->GetYaxis()->SetRange(0, 0);
  Hist3DDBG->GetXaxis()->SetRange(0, 0);

  if (cname.find("VsEta")<50) {
    marker = kFullCircle;
    for (auto ptmin : list) {
      auto ptmax = ptmin + window;
      Hist3DDBG->GetXaxis()->SetRangeUser(ptmin, ptmax);

      string ytitle = "\\sigma (";
      ytitle += Hist3DDBG->GetZaxis()->GetTitle();
      ytitle += ")";
      auto title = Form("_%1.2f_%1.2f_yz", ptmin, ptmax);
      auto aDBG = (TH2F*)Hist3DDBG->Project3D(title);
      aDBG->GetXaxis()->SetRangeUser(0, 0);

      aDBG->FitSlicesX(0, 0, -1, 1);
      auto th1DBG = (TH1F*)gDirectory->Get((std::string(aDBG->GetName()) + value).c_str());
      th1DBG->SetTitle(Form("MC_DBG: %1.2f < p_t < %1.2f", ptmin, ptmax));
      th1DBG->SetMarkerStyle(marker);
      if (colormarker == 5) colormarker ++;
      th1DBG->SetMarkerColor(colormarker);
      th1DBG->SetMarkerSize(marker_size);
      th1DBG->SetLineColorAlpha(colormarker++,0.0);
      th1DBG->SetStats(0);
      th1DBG->SetYTitle(ytitle.c_str());
//      th1DBG->GetXaxis()->SetRangeUser(2.39, 3.6);
//      th1DBG->SetEntryLabel(Form("MC_DBG: %1.2f < p_t < %1.2f", ptmin, ptmax));
//      th1DBG->SetTitle("");
      th1DBG->Draw("same");
    }
  }else if (cname.find("VsPt")<50) {
    marker = kFullSquare;
//  gPad->SetLogy(1);
    for (auto etamin:list) {
      auto etamax = etamin + window;
      Hist3DDBG->GetYaxis()->SetRangeUser(etamin, etamax);

      string ytitle = "\\sigma (";
      ytitle += Hist3DDBG->GetZaxis()->GetTitle();
      ytitle += ")";
      auto title = Form("_%1.2f_%1.2f_xz", etamin, etamax);
      auto aDBG = (TH2F*)Hist3DDBG->Project3D(title);
      aDBG->GetXaxis()->SetRangeUser(0, 0);

      aDBG->FitSlicesX(0, 0, -1, 1);
      auto th1DBG = (TH1F*)gDirectory->Get((std::string(aDBG->GetName()) + value).c_str());
      th1DBG->SetTitle(Form("MC_DBG: %1.2f < \\eta < %1.2f", etamin, etamax));
      th1DBG->SetMarkerStyle(marker);
      if (colormarker == 5) colormarker ++;
      th1DBG->SetMarkerColor(colormarker);
      th1DBG->SetMarkerSize(marker_size);
      th1DBG->SetLineColorAlpha(colormarker++,0.0);
      th1DBG->SetStats(0);
      th1DBG->SetYTitle(ytitle.c_str());
      th1DBG->Draw("same");
    }
  } else {
    cout << " Could not made Histograms Vs Pt/Eta. Please check first argument (canvas name) "; exit(1);
  }
TPaveText *t = new TPaveText(0.2223748,0.9069355,0.7776252,0.965, "brNDC"); // left-up
t->SetBorderSize(0);
t->SetFillColor(gStyle->GetTitleFillColor());
t->AddText(ctitle.c_str());
t->Draw();

  cname += ".png";
//  canva->SetTitle("titulo");
  canva->BuildLegend();
  canva->SetTicky();
  canva->SetGridy();
  canva->Write();
  canva->Print(cname.c_str());

}
//_________________________________________________________________________________________________
TCanvas canvas_summary_3x2(TCanvas* canvas1, TCanvas* canvas2, TCanvas* canvas3, TCanvas* canvas4,
                           TCanvas* canvas5, TCanvas* canvas6, std::string CanvasName) {
                           
  // Create a canvas and divide it
  TCanvas *c1 = new TCanvas(CanvasName.c_str(), CanvasName.c_str(), 1920, 1080);
  c1->UseCurrentStyle();

  // gROOT->SetStyle("Bold");

//std::string tlt = CanvasName.c_str();
  TLatex Title = TLatex();
  Title.SetTextSize(0.035);
  Title.SetTextAlign(23);
  Title.DrawLatex(.5, .995, CanvasName.c_str());
  Title.Draw();

  c1->Divide(1, 2);
  TPad *topPad = (TPad *)c1->cd(1);
  topPad->SetPad(0.000, 0.0, 1, .8);

  TPad *bottomPad = (TPad *)c1->cd(2);
  bottomPad->SetPad(0.000, 0.0, 1, .95);
  bottomPad->Divide(3, 2, 0.005, 0.005);


  // New Pad
  bottomPad->cd(1);
  // gPad->SetTopMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gPad->SetRightMargin(0.10);
  canvas1->DrawClonePad();
  gPad->SetLogy(1);
  
  // New Pad
  bottomPad->cd(2);
  // pad->SetPad(0.33,0.66,0.66,.96);
  gPad->SetBottomMargin(0.15);
  // gPad->SetRightMargin(0.15);
  // gPad->SetLeftMargin(0.15);
  //  gPad->SetLogy(h2log);
  canvas2->DrawClonePad();

  // New Pad
  bottomPad->cd(3);
  // pad->SetPad(0.66,0.66,1,.96);
  gPad->SetBottomMargin(0.15);
  // gPad->SetRightMargin(0.15);
  //  gPad->SetLogy(h3log);
  canvas3->DrawClonePad();

  // New Pad
  bottomPad->cd(4);
  // pad->SetPad(0.000,0.5,0.33,.96);
  gPad->SetBottomMargin(0.15);
  // gPad->SetRightMargin(0.15);
  canvas4->DrawClonePad();
  
  // New Pad
  bottomPad->cd(5);
  gPad->SetBottomMargin(0.15);
  // gPad->SetRightMargin(0.15);
  canvas5->DrawClonePad();

  // New Pad
  bottomPad->cd(6);
  // pad->SetPad(0.000,0.5,0.33,.96);
  gPad->SetBottomMargin(0.15);
  // gPad->SetLeftMargin(0.15);
  canvas6->DrawClonePad();

  c1->cd();

  c1->Update();
  c1->Write();
  return c1;
}
//_________________________________________________________________________________________________
TCanvas canvas_summary_2x2(TCanvas* canvas1, TCanvas* canvas2, TCanvas* canvas3, TCanvas* canvas4,
                           std::string CanvasName) {

  // Create a canvas and divide it
  TCanvas *c1 = new TCanvas(CanvasName.c_str(), CanvasName.c_str(), 1920, 1080);
  c1->UseCurrentStyle();

  // gROOT->SetStyle("Bold");

//std::string tlt = CanvasName.c_str();
  TLatex Title = TLatex();
  Title.SetTextSize(0.035);
  Title.SetTextAlign(23);
  Title.DrawLatex(.5, .995, CanvasName.c_str());
  Title.Draw();

  c1->Divide(1, 2);
  TPad *topPad = (TPad *)c1->cd(1);
  topPad->SetPad(0.000, 0.0, 1, .8);

  TPad *bottomPad = (TPad *)c1->cd(2);
  bottomPad->SetPad(0.000, 0.0, 1, .95);
  bottomPad->Divide(2, 2, 0.005, 0.005);


  // New Pad
  bottomPad->cd(1);
  // gPad->SetTopMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gPad->SetRightMargin(0.10);
  canvas1->DrawClonePad();

  // New Pad
  bottomPad->cd(2);
  // pad->SetPad(0.33,0.66,0.66,.96);
  gPad->SetBottomMargin(0.15);
  // gPad->SetRightMargin(0.15);
  // gPad->SetLeftMargin(0.15);
  //  gPad->SetLogy(h2log);
  canvas2->DrawClonePad();

  // New Pad
  bottomPad->cd(3);
  // pad->SetPad(0.66,0.66,1,.96);
  gPad->SetBottomMargin(0.15);
  // gPad->SetRightMargin(0.15);
  //  gPad->SetLogy(h3log);
  canvas3->DrawClonePad();

  // New Pad
  bottomPad->cd(4);
  // pad->SetPad(0.000,0.5,0.33,.96);
  gPad->SetBottomMargin(0.15);
  // gPad->SetRightMargin(0.15);
  canvas4->DrawClonePad();

  c1->cd();

  c1->Update();
  c1->Write();
  return c1;
}
