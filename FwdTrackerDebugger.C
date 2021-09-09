#if !defined(__CLING__) || defined(__ROOTCLING__)

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "MFTTracking/Cluster.h"
#include "MFTTracking/MFTTrackingParam.h"
#include "ITSMFTSimulation/Hit.h"
#include "Math/SMatrix.h"
#endif

#include "FT3Track.h"
#include "TrackFitter.h"

#pragma link C++ class o2::ft3::FT3Track + ;

using SMatrix55Sym = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using o2::ft3::FT3Track;

enum TH3HistosCodes {
  kFT3TrackDeltaXVertexPtEta,
  kFT3TrackDeltaYVertexPtEta,
  kFT3TrackInvQPtResolutionPtEta,
  kFT3TrackInvQPtPullPtEta,
  kFT3TrackXPullPtEta,
  kFT3TrackYPullPtEta,
  kFT3TrackPhiPullPtEta,
  kFT3TrackTanlPullPtEta,
  kFT3TrackPhi2PullPtEta,
  kFT3TrackReducedChi2PtEta
};

bool DEBUG_VERBOSE = !true;
bool DEBUG_qpt = !true;
bool DEBUG_fitter = !true;

Float_t EtaToTheta(Float_t arg);
Float_t EtaToTanl(Float_t arg);
std::vector<float_t> loadx2X0fromFile(std::string configFileName);

void DrawTrees(const char* treefile = "fwdtrackdebugger.root");
void printCanvas(TCanvas* c, const char* filename);
void exportHisto(TH2F* histo, const char* filename);
void exportHisto(TH1F* histo, const char* filename);
void setSeedCovariances(FT3Track& track);
void th2Hists_EtaPt(string name, std::unique_ptr<TH3F>& Hist3DDBG, std::vector<double> list, double window, string value);

void simFwdTracks(size_t nTracks, float ptMinCut, float ptMax, float etaMin, float etaMax, float zField, std::vector<float> zPositionsVec, std::vector<float> layResVec, TTree* trkDBGTree, std::vector<std::unique_ptr<TH3F>>& TH3Histos);
void simFT3Tracks(size_t nTracks, float ptMinCut, float ptMax, float etaMin, float etaMax, float zField, TTree* trkDBGTree, std::vector<std::unique_ptr<TH3F>>& TH3Histos);

o2::ft3::TrackFitter fitter;
o2::track::TrackParFwd MCTrack_;
o2::mft::TrackMFT probe;
TRandom3 rnd(0);
float pixelPitch;
float clSigma;
bool enableClusterSmearing;

void FwdTrackerDebugger(size_t nTracks = 20000,
                        float zField = -5.0,
                        float clSigma_ = 8.44e-4, //~( o2::itsmft::SegmentationAlpide::PitchCol / std::sqrt(12))
                        bool enableClusterSmearing_ = true,
                        float pixelPitch_ = 0.003, // o2::itsmft::SegmentationAlpide::PitchCol = 0.00292400f
                        const char* treeFile = "fwdtrackdebugger.root")
{
  pixelPitch = pixelPitch_;
  clSigma = clSigma_;
  enableClusterSmearing = enableClusterSmearing_;

  Double_t ptMax = 21.0;
  Double_t ptMin = 0.1;
  Double_t etaMin = -3.8;
  Double_t etaMax = -0.85;
  auto ptmincut = 0.0001;
  fitter.setBz(zField);
  fitter.setLayersx2X0(std::vector<double>(10, 0.0)); // Disable MCS effects
  fitter.mVerbose = DEBUG_fitter;
  //fitter.setTrackModel(o2::mft::MFTTrackModel::Helix);
  //fitter.setTrackModel(o2::mft::MFTTrackModel::Quadratic);

  auto fT = TFile::Open(treeFile, "RECREATE");
  TTree* trkDBGTree = new TTree("fwdTrackDBG", "FwDTrackDBGTree");
  trkDBGTree->Branch("MCTrack", &MCTrack_);
  trkDBGTree->Branch("probe", &probe);

  auto etaInner = -3.6;
  auto etaCenter = -3.06;
  auto etaOuter = -2.45;
  auto etaTransition = -0.881;

  std::vector<std::unique_ptr<TH3F>> TH3Histos;

  std::map<int, const char*> TH3Names{
    {kFT3TrackDeltaXVertexPtEta, "FT3DBGTrackDeltaXVertexPtEta"},
    {kFT3TrackDeltaYVertexPtEta, "FT3DBGTrackDeltaYVertexPtEta"},
    {kFT3TrackInvQPtResolutionPtEta, "FT3DBGTrackInvQPtResolutionPtEta"},
    {kFT3TrackInvQPtPullPtEta, "FT3DBGTrackInvQPtPullPtEta"},
    {kFT3TrackXPullPtEta, "FT3DBGTrackXPullPtEta"},
    {kFT3TrackYPullPtEta, "FT3DBGTrackYPullPtEta"},
    {kFT3TrackPhiPullPtEta, "FT3DBGPhiPullPtEta"},
    {kFT3TrackPhi2PullPtEta, "FT3DBGPhi2PullPtEta"},
    {kFT3TrackTanlPullPtEta, "FT3DBGTrackTanlPullPtEta"},
    {kFT3TrackInvQPtPullPtEta, "FT3DBGTrackInvQPtPullPtEta"},
    {kFT3TrackReducedChi2PtEta, "FT3TrackReducedChi2PtEta"}};

  std::map<int, const char*> TH3Titles{
    {kFT3TrackDeltaXVertexPtEta, "FT3DBGTrackDeltaXVertexPtEta"},
    {kFT3TrackDeltaYVertexPtEta, "FT3DBGTrackDeltaYVertexPtEta"},
    {kFT3TrackInvQPtResolutionPtEta, "FT3DBGTrackInvQPtResolutionPtEta;P_t;\\eta;Res Q/pt"},
    {kFT3TrackInvQPtPullPtEta, "FT3DBGTrackInvQPtPullPtEta"},
    {kFT3TrackXPullPtEta, "FT3DBGTrackXPullPtEta;p_t;\\eta;(\\Delta x)/\\sigma_{x}"},
    {kFT3TrackYPullPtEta, "FT3DBGTrackYPullPtEta;p_t;\\eta;(\\Delta y)/\\sigma_{y}"},
    {kFT3TrackPhiPullPtEta, "FT3DBGTrackPhiPullPtEta;P_t;\\eta;(\\Delta \\phi)/\\sigma_{\\phi}"},
      {kFT3TrackPhi2PullPtEta, "FT3DBGPhi2PullPtEta;p_t;\\eta;(\\Delta Phi2)/\\sigma_{Phi}"},
    {kFT3TrackTanlPullPtEta, "FT3DBGTrackTanlPullPtEta;P_t;\\eta;(\\Delta tan\\lambda)/\\sigma_{tan\\lambda}"},
    {kFT3TrackInvQPtPullPtEta, "FT3DBGTrackInvQPtPullPtEta;P_t;\\eta;(\\Delta q/pt)/\\sigma_{q/pt}"},
    {kFT3TrackReducedChi2PtEta, "FT3DBGTrackReducedChi2PtEta;P_t;\\eta;\\Reduced Chi2"}};

  std::map<int, std::array<double, 9>> TH3Binning{
    {kFT3TrackDeltaXVertexPtEta, {50, 0, 10, 84, 1.7, 4.5, 500, -500, 500}},
    {kFT3TrackDeltaYVertexPtEta, {50, 0, 10, 84, 1.7, 4.5, 500, -500, 500}},
    {kFT3TrackInvQPtResolutionPtEta, {200, 0, 20, 28, 1.7, 4.5, 100, -1, 1}},
    {kFT3TrackXPullPtEta, {200, 0, 20, 28, 1.7, 4.5, 100, -10, 10}},
    {kFT3TrackYPullPtEta, {200, 0, 20, 28, 1.7, 4.5, 100, -10, 10}},
    {kFT3TrackPhiPullPtEta, {200, 0, 20, 28, 1.7, 4.5, 100, -10, 10}},
    {kFT3TrackPhi2PullPtEta, {50, 0, 20, 28, 1.7, 4.5, 100, -10, 50}},
    {kFT3TrackTanlPullPtEta, {200, 0, 20, 28, 1.7, 4.5, 100, -10, 10}},
    {kFT3TrackInvQPtPullPtEta, {200, 0, 20, 28, 1.7, 4.5, 100, -10, 10}},
    {kFT3TrackReducedChi2PtEta, {200, 0, 20, 28, 1.7, 4.5, 50, -0.5, 30}}};

  std::map<int, const char*> TH3XaxisTitles{
    {kFT3TrackDeltaXVertexPtEta, "p_t"},
    {kFT3TrackDeltaYVertexPtEta, "p_t"},
    {kFT3TrackInvQPtPullPtEta, "p_t"},
    {kFT3TrackInvQPtResolutionPtEta, "p_t"}};

  std::map<int, const char*> TH3YaxisTitles{
    {kFT3TrackDeltaXVertexPtEta, "\\eta"},
    {kFT3TrackDeltaYVertexPtEta, "\\eta"},
    {kFT3TrackInvQPtPullPtEta, "\\eta"},
    {kFT3TrackInvQPtResolutionPtEta, "\\eta"}};

  std::map<int, const char*> TH3ZaxisTitles{
    {kFT3TrackDeltaXVertexPtEta, "X residual at vertex (um)"},
    {kFT3TrackDeltaYVertexPtEta, "Y residual at vertex (um)"},
    {kFT3TrackInvQPtPullPtEta, "(\\Delta q/p_t)/\\sigma_{q/pt}"},
    {kFT3TrackInvQPtResolutionPtEta, "(q/p_t residual)/(q/pt)"}};

  const int nTH3Histos = TH3Names.size();

  std::cout << " nTH3Histos = " << nTH3Histos << std::endl;
  //for (auto& h : TH3Histos) {
  for (auto n3Histo = 0; n3Histo < nTH3Histos; n3Histo++) {

    TH3Histos.emplace_back(std::make_unique<TH3F>(TH3Names[n3Histo], TH3Titles[n3Histo],
                                                  (int)TH3Binning[n3Histo][0],
                                                  TH3Binning[n3Histo][1],
                                                  TH3Binning[n3Histo][2],
                                                  (int)TH3Binning[n3Histo][3],
                                                  TH3Binning[n3Histo][4],
                                                  TH3Binning[n3Histo][5],
                                                  (int)TH3Binning[n3Histo][6],
                                                  TH3Binning[n3Histo][7],
                                                  TH3Binning[n3Histo][8]));
  }

  simFT3Tracks(nTracks, ptmincut, ptMax, etaMin, etaMax, zField, trkDBGTree, TH3Histos);

  fT->mkdir("MoreHistos");
  fT->cd("MoreHistos");

  for (auto& h : TH3Histos) {
    h->Write();
  }

  fT->Write();
}

Float_t EtaToTheta(Float_t arg)
{
  return (180. / TMath::Pi()) * 2. * atan(exp(-arg));
}

Float_t EtaToTanl(Float_t eta)
{
  return tan(TMath::Pi() / 2 - 2 * atan(exp(-eta)));
}

void simFT3Tracks(size_t nTracks, float ptMinCut, float ptMax, float etaMin, float etaMax, float zField, TTree* trkDBGTree, std::vector<std::unique_ptr<TH3F>>& TH3Histos)
{
  std::vector<float> zPositionsVec{-16., -20., -24., -77., -100., -122., -150., -180., -220., -279.};
  std::vector<float> layResVec(zPositionsVec.size(), clSigma);
  simFwdTracks(nTracks, ptMinCut, ptMax, etaMin, etaMax, zField, zPositionsVec, layResVec, trkDBGTree, TH3Histos);
}

void simFwdTracks(size_t nTracks, float ptMinCut, float ptMax, float etaMin, float etaMax, float zField, std::vector<float> zPositionsVec, std::vector<float> layResVec, TTree* trkDBGTree, std::vector<std::unique_ptr<TH3F>>& TH3Histos)
{

//just a quick chi2 histogram
  TH1D * chi_hist = new TH1D("h1_chi2","chi2",100,0,50);
    chi_hist->SetFillColor(kRed);
    chi_hist->SetFillStyle(3005);
  TH1D * chiRed_hist = new TH1D("h1_redu_chi2","reduced_chi2",100,0,6);
    chiRed_hist->SetFillColor(kRed);
    chiRed_hist->SetFillStyle(3005);

  for (size_t i = 0; i < nTracks; i++) {
    if (i > 10)
      fitter.mVerbose = !true;
    o2::track::TrackParFwd MCTrack;
    MCTrack.setX(0);
    MCTrack.setY(0);
    MCTrack.setPhi(rnd.Uniform(-TMath::Pi(), TMath::Pi()));
    //MCTrack.setPhi(TMath::Pi() / 4);

    auto eta = rnd.Uniform(etaMin, etaMax);

    MCTrack.setTanl(EtaToTanl(eta));
    auto qpt = 0.0;
    while (1) {
      qpt = rnd.Uniform(-ptMax, ptMax);
      if (std::abs(qpt) > ptMinCut)
        break;
    }
    MCTrack.setInvQPt(1.0 / qpt);
    MCTrack_ = MCTrack;
    if (DEBUG_VERBOSE) {
      o2::track::TrackParFwd tempTrk = MCTrack;
      auto z = zPositionsVec.back();
      tempTrk.propagateParamToZhelix(z, zField);
      std::cout << "\n\n *** New gen track ***\n  MCTrack: X = " << MCTrack.getX() << " Y = " << MCTrack.getY() << " Z = " << MCTrack.getZ() << " Tgl = " << MCTrack.getTanl() << "  Phi = " << MCTrack.getPhi() << "  q/pt = " << MCTrack.getInvQPt() << std::endl;
      std::cout << "    MCTrack at last disk: X = " << tempTrk.getX() << " Y = " << tempTrk.getY() << " Z = " << tempTrk.getZ() << " Tgl = " << tempTrk.getTanl() << "  Phi = " << tempTrk.getPhi() << "  (" << o2::math_utils::toPMPiGen(tempTrk.getPhi()) << ") q/pt = " << tempTrk.getInvQPt() << std::endl;
    }
    FT3Track ft3ProbeTr;
    int nLayer = 0;
    for (auto z : zPositionsVec) {
      MCTrack.propagateParamToZhelix(z, zField);
      //std::cout << " AddHit: MCTrack.getX() = " << MCTrack.getX() << " ; MCTrack.getY() =  " << MCTrack.getY() << "  MCTrack.getZ() = " << MCTrack.getZ() << std::endl;

      o2::itsmft::Hit hit(0, 0, {MCTrack.getX(), MCTrack.getY(), MCTrack.getZ()}, {0, 0, 0}, {0, 0, 0}, 0, 0, 0, 0, 0);
      ft3ProbeTr.addHit(hit, nLayer, clSigma);
      nLayer++;
    }
    //std::cout << std::endl;
    //if (DEBUG_VERBOSE)
    // std::cout << std::endl;
    ft3ProbeTr.sort();

    fitter.initTrack(ft3ProbeTr);
    setSeedCovariances(ft3ProbeTr);
//    fitter.MinuitFit(ft3ProbeTr); 
    fitter.fit(ft3ProbeTr);
    if (DEBUG_qpt)
      std::cout << "Track " << i << ": \n    q/pt_MC = " << MCTrack.getInvQPt() << "\n    q/pt_Seed = " << ft3ProbeTr.getInvQPtSeed() << "\n    q/pt_fit = " << ft3ProbeTr.getInvQPt() << std::endl;
    ft3ProbeTr.propagateToZhelix(0.0, zField);
    if (DEBUG_VERBOSE)
      std::cout << "  Track at vertex: X = " << ft3ProbeTr.getX() << " Y = " << ft3ProbeTr.getY() << " Z = " << ft3ProbeTr.getZ() << " Tgl = " << ft3ProbeTr.getTanl() << "  Phi = " << ft3ProbeTr.getPhi() << "  (" << o2::math_utils::toPMPiGen(ft3ProbeTr.getPhi()) << ")  q/pt = " << ft3ProbeTr.getInvQPt() << std::endl;

    probe = (o2::mft::TrackMFT)ft3ProbeTr;
    trkDBGTree->Fill();

    auto vx_MC = MCTrack_.getX();
    auto vy_MC = MCTrack_.getY();
    auto invQPt_MC = MCTrack_.getInvQPt();
    auto Pt_MC = MCTrack_.getPt();
    auto eta_MC = MCTrack_.getEta();
    auto phi_MC = MCTrack_.getPhi();
    auto tanl_MC = MCTrack_.getTanl();

    auto dx = ft3ProbeTr.getX() - vx_MC;
    auto dy = ft3ProbeTr.getY() - vy_MC;
    auto invQPt_Fit = ft3ProbeTr.getInvQPt();
    auto d_invQPt = invQPt_Fit - invQPt_MC;
    auto d_Phi = ft3ProbeTr.getPhi() - phi_MC;
    auto d_tanl = ft3ProbeTr.getTanl() - tanl_MC;
    auto trackChi2 = ft3ProbeTr.getTrackChi2();
    auto nClusters = ft3ProbeTr.getNumberOfPoints();

    TH3Histos[kFT3TrackDeltaXVertexPtEta]->Fill(Pt_MC, std::abs(eta_MC), 1e4 * dx);
    TH3Histos[kFT3TrackDeltaYVertexPtEta]->Fill(Pt_MC, std::abs(eta_MC), 1e4 * dy);
    TH3Histos[kFT3TrackInvQPtResolutionPtEta]->Fill(Pt_MC, std::abs(eta_MC), (invQPt_Fit - invQPt_MC) / invQPt_MC);
    TH3Histos[kFT3TrackXPullPtEta]->Fill(Pt_MC, std::abs(eta_MC), dx / sqrt(ft3ProbeTr.getCovariances()(0, 0)));
    TH3Histos[kFT3TrackYPullPtEta]->Fill(Pt_MC, std::abs(eta_MC), dy / sqrt(ft3ProbeTr.getCovariances()(1, 1)));
    TH3Histos[kFT3TrackPhiPullPtEta]->Fill(Pt_MC, std::abs(eta_MC), d_Phi / sqrt(ft3ProbeTr.getCovariances()(2, 2)));
    TH3Histos[kFT3TrackTanlPullPtEta]->Fill(Pt_MC, std::abs(eta_MC), d_tanl / sqrt(ft3ProbeTr.getCovariances()(3, 3)));
    TH3Histos[kFT3TrackInvQPtPullPtEta]->Fill(Pt_MC, std::abs(eta_MC), d_invQPt / sqrt(ft3ProbeTr.getCovariances()(4, 4)));
    TH3Histos[kFT3TrackReducedChi2PtEta]->Fill(Pt_MC, std::abs(eta_MC), trackChi2 / (2*nClusters - 5)); // 5: number of fitting paramet
    chi_hist->Fill(trackChi2);
    chiRed_hist->Fill(trackChi2/(2*nClusters - 5));
  }
    chi_hist->Draw();
    chiRed_hist->Draw();

  // InvPtResolution vs eta MC Debuger
  auto& FT3TrackInvQPtResolutionPtEtaDBG = TH3Histos[kFT3TrackInvQPtResolutionPtEta];
  FT3TrackInvQPtResolutionPtEtaDBG->GetYaxis()->SetRange(0, 0);

  int marker = kFullCircle;
  std::vector<double> ptList({1, 5, 9.});
  double ptWindow = 1.;
  std::vector<double> etaList({2., 3., 3.8});
  double etaWindow = 0.2;

/* /last argument is a string: _1 for mean, _2 for std dev
  th2Hists_EtaPt("InvQPtResVsEta", TH3Histos[kFT3TrackInvQPtResolutionPtEta], ptList, ptWindow, "_2");
  th2Hists_EtaPt("InvQPtResVsPt", TH3Histos[kFT3TrackInvQPtResolutionPtEta], etaList, etaWindow, "_2");
  TH3Histos[kFT3TrackInvQPtResolutionPtEta]->GetYaxis()->SetRangeUser(0, 0);
  TH3Histos[kFT3TrackInvQPtResolutionPtEta]->GetXaxis()->SetRangeUser(0, 0);

  th2Hists_EtaPt("Redu_Chi2VsEta", TH3Histos[kFT3TrackReducedChi2PtEta], ptList, ptWindow, "_2");
  th2Hists_EtaPt("Redu_Chi2VsPt", TH3Histos[kFT3TrackReducedChi2PtEta], etaList, etaWindow, "_2");
  TH3Histos[kFT3TrackReducedChi2PtEta]->GetYaxis()->SetRangeUser(0, 0);
  TH3Histos[kFT3TrackReducedChi2PtEta]->GetXaxis()->SetRangeUser(0, 0);

  th2Hists_EtaPt("PhiPullVsEta", TH3Histos[kFT3TrackPhiPullPtEta], ptList, ptWindow, "_2");
  th2Hists_EtaPt("PhiPullVsPt", TH3Histos[kFT3TrackPhiPullPtEta], etaList, etaWindow, "_2");
  TH3Histos[kFT3TrackPhiPullPtEta]->GetYaxis()->SetRangeUser(0, 0);
  TH3Histos[kFT3TrackPhiPullPtEta]->GetXaxis()->SetRangeUser(0, 0);

  th2Hists_EtaPt("TanlPullVsEta", TH3Histos[kFT3TrackTanlPullPtEta], ptList, ptWindow, "_2");
  th2Hists_EtaPt("TanlPullVsPt", TH3Histos[kFT3TrackTanlPullPtEta], etaList, etaWindow, "_2");
  TH3Histos[kFT3TrackTanlPullPtEta]->GetYaxis()->SetRangeUser(0, 0);
  TH3Histos[kFT3TrackTanlPullPtEta]->GetXaxis()->SetRangeUser(0, 0);

  th2Hists_EtaPt("DeltaXVertexVsEta",TH3Histos[kFT3TrackDeltaXVertexPtEta],ptList, ptWindow, "_2");
  th2Hists_EtaPt("DeltaXVertexVsPt",TH3Histos[kFT3TrackDeltaXVertexPtEta],etaList, etaWindow, "_2");
  TH3Histos[kFT3TrackDeltaXVertexPtEta]->GetYaxis()->SetRangeUser(0, 0);
  TH3Histos[kFT3TrackDeltaXVertexPtEta]->GetXaxis()->SetRangeUser(0, 0);

  th2Hists_EtaPt("DeltaYVertexVsEta",TH3Histos[kFT3TrackDeltaYVertexPtEta],ptList, ptWindow, "_2");
  th2Hists_EtaPt("DeltaYVertexVsPt",TH3Histos[kFT3TrackDeltaYVertexPtEta],etaList, etaWindow, "_2");
  TH3Histos[kFT3TrackDeltaYVertexPtEta]->GetYaxis()->SetRangeUser(0, 0);
  TH3Histos[kFT3TrackDeltaYVertexPtEta]->GetXaxis()->SetRangeUser(0, 0);

  th2Hists_EtaPt("InvQPtPullVsEta", TH3Histos[kFT3TrackInvQPtPullPtEta], ptList, ptWindow, "_2");
  th2Hists_EtaPt("InvQPtPullVsPt", TH3Histos[kFT3TrackInvQPtPullPtEta], etaList, etaWindow, "_2");
  TH3Histos[kFT3TrackInvQPtPullPtEta]->GetYaxis()->SetRangeUser(0, 0);
  TH3Histos[kFT3TrackInvQPtPullPtEta]->GetXaxis()->SetRangeUser(0, 0);

  th2Hists_EtaPt("XPullVsEta", TH3Histos[kFT3TrackXPullPtEta], ptList, ptWindow, "_2");
  th2Hists_EtaPt("XPullVsPt", TH3Histos[kFT3TrackXPullPtEta], etaList, etaWindow, "_2");
  TH3Histos[kFT3TrackXPullPtEta]->GetYaxis()->SetRangeUser(0, 0);
  TH3Histos[kFT3TrackXPullPtEta]->GetXaxis()->SetRangeUser(0, 0);

  th2Hists_EtaPt("YPullVsEta", TH3Histos[kFT3TrackYPullPtEta], ptList, ptWindow, "_2");
  th2Hists_EtaPt("YPullVsPt", TH3Histos[kFT3TrackYPullPtEta], etaList, etaWindow, "_2");
*/
}

void th2Hists_EtaPt(string name, std::unique_ptr<TH3F>& Hist3DDBG, std::vector<double> list, double window, string value)
{
  auto canva = new TCanvas(name.c_str(), name.c_str(), 1080, 1080);
  int marker = kFullCircle;
  int colormarker = 1;

  Hist3DDBG->GetYaxis()->SetRange(0, 0);
  Hist3DDBG->GetXaxis()->SetRange(0, 0);
  if (name.find("VsEta")<50) {
    for (auto ptmin : list) {
      auto ptmax = ptmin + window;
      Hist3DDBG->GetXaxis()->SetRangeUser(ptmin, ptmax);

      auto title = Form("_%1.2f_%1.2f_yz", ptmin, ptmax);
      auto aDBG = (TH2F*)Hist3DDBG->Project3D(title);
//      aDBG->GetXaxis()->SetRangeUser(0, 0);

      aDBG->FitSlicesX(0, 0, -1, 1);
      auto th2DBG = (TH1F*)gDirectory->Get((std::string(aDBG->GetName()) + value).c_str());
      th2DBG->SetTitle(Form("MC_DBG: %1.2f < p_t < %1.2f (w.o.MCS)", ptmin, ptmax));
      th2DBG->SetMarkerStyle(marker++);
      th2DBG->SetMarkerColor(colormarker);
      th2DBG->SetLineColorAlpha(colormarker++,0.0);
      th2DBG->SetStats(0);
      th2DBG->Draw("same");
      Hist3DDBG->GetXaxis()->SetRangeUser(0, 0);
    }
  }else if (name.find("VsPt")<50) {
    for (auto etamin:list) {
      auto etamax = etamin + window;
      Hist3DDBG->GetYaxis()->SetRangeUser(etamin, etamax);

      auto title = Form("_%1.2f_%1.2f_xz", etamin, etamax);
      auto aDBG = (TH2F*)Hist3DDBG->Project3D(title);
//      aDBG->GetXaxis()->SetRangeUser(0, 0);

      aDBG->FitSlicesX(0, 0, -1, 1);
      auto th2DBG = (TH1F*)gDirectory->Get((std::string(aDBG->GetName()) + value).c_str());
      th2DBG->SetTitle(Form("MC_DBG: %1.2f < \\eta < %1.2f (w.o.MCS)", etamin, etamax));
      th2DBG->SetMarkerStyle(marker++);
      th2DBG->SetMarkerColor(colormarker);
      th2DBG->SetLineColorAlpha(colormarker++,0.0);
      th2DBG->SetStats(0);
      th2DBG->Draw("same");
    }
  } else {
    cout << " SOMETHINGWENTWRONG! GOCHECKITOUT! "; exit(1);
  }

  
  name += ".png";
  canva->BuildLegend();
  canva->SetTicky();
  canva->SetGridy();
  canva->Write();
  canva->Print(name.c_str());
}

void setSeedCovariances(FT3Track& track)
{

  auto tan_k = 10.0;
  auto q2pt_c = 10.;
  auto phi_c = 1 / 16.;

  SMatrix55Sym Covariances{}; ///< \brief Covariance matrix of track parameters
  float q2ptcov = TMath::Max(std::abs(track.getInvQPt()), .5);
  float tanlerr = TMath::Max(std::abs(track.getTanl()), .5);

  Covariances(0, 0) = 1;
  Covariances(1, 1) = 1;
  Covariances(2, 2) = phi_c * TMath::Pi() * TMath::Pi();
  Covariances(3, 3) = tan_k * tanlerr * tanlerr;
  Covariances(4, 4) = q2pt_c * q2ptcov * q2ptcov;
  track.setCovariances(Covariances);
}

void printCanvas(TCanvas* c, const char* filename)
{
  //gStyle->SetImageScaling(3.);
  c->SetBatch();
  gSystem->ProcessEvents();
  c->Update();
  c->Print(filename);
}

void exportHisto(TH2F* histo, const char* filename)
{
  TCanvas* C1_ = new TCanvas(filename, filename, 1080, 1080);
  histo->Draw("colz");
  gSystem->ProcessEvents();
  C1_->Update();
  C1_->Print(filename);
}

void exportHisto(TH1F* histo, const char* filename)
{
  TCanvas* C1_ = new TCanvas(filename, filename, 1080, 1080);
  histo->Draw("");
  gSystem->ProcessEvents();
  C1_->Update();
  C1_->Print(filename);
}

std::vector<float_t> loadx2X0fromFile(std::string configFileName = "FT3_layout.cfg")
{
  std::vector<float_t> Layersx2X0;
  std::ifstream ifs(configFileName.c_str());
  if (!ifs.good()) {
    LOG(FATAL) << " Invalid FT3Base.configFile!";
  }
  std::string tempstr;
  float z_layer, r_in, r_out, Layerx2X0;
  char delimiter;
  int layerNumber = 0;
  while (std::getline(ifs, tempstr)) {
    if (tempstr[0] == '#') {
      LOG(INFO) << " Comment: " << tempstr;
      continue;
    }
    LOG(INFO) << " Line: " << tempstr;
    std::istringstream iss(tempstr);
    iss >> z_layer;
    iss >> r_in;
    iss >> r_out;
    iss >> Layerx2X0;

    Layersx2X0.push_back(Layerx2X0);
    LOG(INFO) << " loadx2X0fromFile z =  " << z_layer << " ; x/X0 = " << Layerx2X0 << std::endl;
  }
  return Layersx2X0;
}
