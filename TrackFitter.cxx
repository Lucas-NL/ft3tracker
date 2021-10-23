// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file TrackFitter.cxx

#include "TrackFitter.h"
#include "CommonConstants/MathConstants.h"
#include "MathUtils/Utils.h"
#include "MathUtils/fit.h"
#include <TF1.h>
#include <TF2.h>
#include <TMath.h>
#include <TMatrixD.h>

using o2::math_utils::fitGaus;

ROOT::Math::SMatrix<double,10,10,ROOT::Math::MatRepSym<double,10> > invCy_10, invCx_10;
ROOT::Math::SVector<double, 10> xPositionsTrack_10, yPositionsTrack_10;
ROOT::Math::SMatrix<double, 10, 3> G_10;

namespace o2
{
namespace ft3
{

  double k;// = TMath::Abs(o2::constants::math::B2C * TrackFitter::mBZField);
  int Hz;// = std::copysign(1, TrackFitter::mBZField);
//_________________________________________________________________________________________________
void TrackFitter::setBz(Double_t bZ)
{

  /// Set the constant magnetic field
  mBZField = bZ;

  if (mVerbose) {
    LOG(INFO) << "Setting Fitter field = " << bZ;
  }
}

//_________________________________________________________________________________________________
bool TrackFitter::fit(FT3Track& track, bool outward)
{

  /// Fit a track using its attached clusters
  /// Returns false in case of failure

  auto nClusters = track.getNumberOfPoints();

  if (mVerbose) {
    std::cout << "Seed Parameters: X = " << track.getX() << " Y = " << track.getY() << " Z = " << track.getZ() << " Tgl = " << track.getTanl() << "  Phi = " << track.getPhi() << " q/pt = " << track.getInvQPt() << std::endl;
    std::cout << "Seed covariances: \n"
              << track.getCovariances() << std::endl
              << std::endl;
  }

  // recursively compute clusters, updating the track parameters
  if (!outward) { // Inward for vertexing
    //nClusters--;
    while (nClusters-- > 0) {
      if (!computeCluster(track, nClusters)) {
        return false;
      }
    }
  } else {       // Outward for MCH matching
    int ncl = 0; //1;
    while (ncl < nClusters) {
      if (!computeCluster(track, ncl)) {
        return false;
      }
      ncl++;
    }
  }
  if (mVerbose) {
    //  Print final covariances? std::cout << "Track covariances:";
    //  track->getCovariances().Print();
    std::cout << "Track Chi2 = " << track.getTrackChi2() << std::endl;
    cout << "\n[nClusters = " << track.getNumberOfPoints() << "] LayerIDs, zCoord => ";
    for (auto i = 0; i < track.getNumberOfPoints(); i++) {
      std::cout << " " << track.getLayers()[i] << ", "
                << track.getZCoordinates()[i] << " ";
    }
    std::cout << std::endl
              << std::endl
              << std::endl;
    std::cout << " ***************************** Done fitting *****************************\n";
  }

  return true;
}

//_________________________________________________________________________________________________
bool TrackFitter::initTrack(FT3Track& track, bool outward)
{

  // initialize the starting track parameters
  Double_t x0;
  Double_t y0;
  Double_t z0;
  Double_t sigmainvQPtsq;
  Double_t chi2invqptquad = 0;
  auto invQPt0 = invQPtFromFCF(track, mBZField, sigmainvQPtsq);

  auto nPoints = track.getNumberOfPoints();
  k = TMath::Abs(o2::constants::math::B2C * mBZField);
  Hz = std::copysign(1, mBZField);

  track.setInvQPtSeed(invQPt0);
  track.setChi2QPtSeed(chi2invqptquad);
  track.setInvQPt(invQPt0);

  /// Compute the initial track parameters to seed the Kalman filter
  int first_cls, last_cls;
  if (outward) { // MCH matching
    first_cls = 0;
    last_cls = 1;
    x0 = track.getXCoordinates()[0];
    y0 = track.getYCoordinates()[0];
    z0 = track.getZCoordinates()[0];
  } else { // Vertexing
    x0 = track.getXCoordinates()[nPoints - 1];
    y0 = track.getYCoordinates()[nPoints - 1];
    z0 = track.getZCoordinates()[nPoints - 1];
    first_cls = nPoints - 2;
    last_cls = nPoints - 1;
  }

  //Compute tanl using first two clusters
  auto deltaX = track.getXCoordinates()[1] - track.getXCoordinates()[0];
  auto deltaY = track.getYCoordinates()[1] - track.getYCoordinates()[0];
  auto deltaZ = track.getZCoordinates()[1] - track.getZCoordinates()[0];
  auto deltaR = TMath::Sqrt(deltaX * deltaX + deltaY * deltaY);
  auto tanl0 = 0.5 * TMath::Sqrt2() * (deltaZ / deltaR) *
               TMath::Sqrt(TMath::Sqrt((invQPt0 * deltaR * k) * (invQPt0 * deltaR * k) + 1) + 1);

  // Compute phi at the last cluster using two last clusters
  deltaX = track.getXCoordinates()[last_cls] - track.getXCoordinates()[first_cls];
  deltaY = track.getYCoordinates()[last_cls] - track.getYCoordinates()[first_cls];
  deltaZ = track.getZCoordinates()[last_cls] - track.getZCoordinates()[first_cls];
  deltaR = TMath::Sqrt(deltaX * deltaX + deltaY * deltaY);
  auto phi0 = TMath::ATan2(deltaY, deltaX) - 0.5 * Hz * invQPt0 * deltaZ * k / tanl0;

  track.setX(x0);
  track.setY(y0);
  track.setZ(z0);
  track.setPhi(phi0);
  track.setTanl(tanl0);
  if (mVerbose) {
    std::cout << "  Init track: X = " << track.getX() << " Y = " << track.getY() << " Z = " << track.getZ() << " Tgl = " << track.getTanl() << "  Phi = " << track.getPhi() << "  (" << o2::math_utils::toPMPiGen(track.getPhi()) << ") q/pt = " << track.getInvQPt() << std::endl;
  }
  SMatrix55Sym lastParamCov;
  Double_t qptsigma = TMath::Max(std::abs(track.getInvQPt()), .5);
  Double_t tanlsigma = TMath::Max(std::abs(track.getTanl()), .5);

  lastParamCov(0, 0) = 1;                              // <X,X>
  lastParamCov(1, 1) = 1;                              // <Y,X>
  lastParamCov(2, 2) = TMath::Pi() * TMath::Pi() / 16; // <PHI,X>
  lastParamCov(3, 3) = 10 * tanlsigma * tanlsigma;     // <TANL,X>
  lastParamCov(4, 4) = 10 * qptsigma * qptsigma;       // <INVQPT,X>

  track.setCovariances(lastParamCov);
  track.setTrackChi2(0.);

  return true;
}

//_________________________________________________________________________________________________
bool TrackFitter::computeCluster(FT3Track& track, int cluster)
{
  /// Propagate track to the z position of the new cluster
  /// accounting for MCS dispersion in the current layer and the other(s)
  /// crossed Recompute the parameters adding the cluster constraint with the
  /// Kalman filter Returns false in case of failure

  const auto& clx = track.getXCoordinates()[cluster];
  const auto& cly = track.getYCoordinates()[cluster];
  const auto& clz = track.getZCoordinates()[cluster];
  const auto& sigmaX2 = track.getSigmasX2()[cluster];
  const auto& sigmaY2 = track.getSigmasY2()[cluster];

  /*if (track.getZ() == clz) {
    LOG(INFO) << "AddCluster ERROR: The new cluster must be upstream!"
              << (track.isCA() ? " CATrack" : "LTFTrack");
    LOG(INFO) << "track.getZ() = " << track.getZ() << " ; newClusterZ = " << clz
              << " ==> Skipping point.";
    return true;
  }*/
  if (mVerbose) {
    std::cout << "computeCluster:     X = " << clx << " Y = " << cly
              << " Z = " << clz << " nCluster = " << cluster << std::endl;
  }

  if (mVerbose) {
    std::cout << "  BeforeExtrap: X = " << track.getX() << " Y = " << track.getY() << " Z = " << track.getZ() << " Tgl = " << track.getTanl() << "  Phi = " << track.getPhi() << " q/pt = " << track.getInvQPt() << std::endl;
  }

  // Propagate track to the z position of the new cluster
  track.propagateToZhelix(clz, mBZField);
  //track.propagateToZ(clz, mBZField);

  if (mVerbose) {
    std::cout << "   AfterExtrap: X = " << track.getX() << " Y = " << track.getY() << " Z = " << track.getZ() << " Tgl = " << track.getTanl() << "  Phi = " << track.getPhi() << " q/pt = " << track.getInvQPt() << std::endl;
    std::cout << " Track covariances after extrap: \n"
              << track.getCovariances() << std::endl
              << ((track.getCovariances()(4, 4) < 0.) ? " NEGATIVE q/Pt VARIANCE!" : "")
              << std::endl;
  }

  // add MCS effects for the new cluster
  auto Layerx2X0 = mLayersx2X0[cluster];
  //track.addMCSEffect(0.5 * Layerx2X0);

  //if (mVerbose) {
  //  std::cout << "   After MCS Effects I:" << std::endl
  //            << track.getCovariances() << std::endl
  //            << std::endl;
  //}

  // recompute parameters
  const std::array<Float_t, 2>& pos = {clx, cly};
  const std::array<Float_t, 2>& cov = {sigmaX2, sigmaY2};

  if (track.update(pos, cov)) {
    if (mVerbose) {
      std::cout << "   New Cluster: X = " << clx << " Y = " << cly << " Z = " << clz << std::endl;
      std::cout << "   AfterKalman: X = " << track.getX() << " Y = " << track.getY() << " Z = " << track.getZ() << " Tgl = " << track.getTanl() << "  Phi = " << track.getPhi() << " q/pt = " << track.getInvQPt() << std::endl;
      std::cout << "   Track covariances after Kalman update: \n"
                << track.getCovariances() << std::endl
                << std::endl;
    }
    track.addMCSEffect(1.0 * Layerx2X0);

    if (mVerbose) {
      std::cout << "  After MCS Effects II:  mLayersx2X0[cluster] = " << Layerx2X0 << std::endl;
      std::cout << " " << track.getCovariances() << std::endl
                << std::endl;
    }
    return true;
  }
  return false;
}

//_________________________________________________________________________________________________
void TrackFitter::MinuitFit(FT3Track& track)
{
  initTrack(track, 1);

  TrackFitter::PosX = track.getXCoordinates();
  TrackFitter::PosY = track.getYCoordinates();
  TrackFitter::PosZ = track.getZCoordinates();
  TrackFitter::ErrorsX = track.getSigmasX2();
  TrackFitter::ErrorsY = track.getSigmasY2();

  int aux = 0;
  bool debug_minuit = !true;
  bool MCSeffects = true;

  if (debug_minuit) {
    std::ofstream log("MinuitFitter_Log",std::ofstream::app);
    log << " Track Coordinates: "<<std::endl;
    for (auto a:TrackFitter::PosZ) {
      log << " Z Pos: " << a << " X Pos: " << TrackFitter::PosX[aux] << " Y Pos: " << TrackFitter::PosY[aux] << endl;
      aux ++;
    }
    log << " Track done. Number of clusters: "<< track.getNumberOfPoints() << " One more check (aux count): " << aux <<" \n "<<std::endl;
  }

  auto csclambda = TMath::Abs(TMath::Sqrt(1 + track.getTanl() * track.getTanl()) / track.getTanl());
  auto pathLengthOverX0 = mLayersx2X0[0] * csclambda; //
  //  Angular dispersion square of the track (variance) in a plane perpendicular to the trajectory
  sigmathetasq = 0.0136 * track.getInverseMomentum();
  sigmathetasq *= sigmathetasq * pathLengthOverX0;

//for (int i=0;i<10;i++) {
//      log << " Z Pos: " << PosZ[i] << " X Pos: " << TrackFitter::PosX[i] << endl;
//    }
 auto nhits = track.getNumberOfPoints();
 if (nhits == 10) {
   for (int n = 0; n < nhits; n++) {
      auto z = PosZ[n] - PosZ[0];
      G_10(n,0) = 1;
      G_10(n,1) = z;
      G_10(n,2) = z * z * 0.5;
   } 
   ROOT::Math::SMatrix<double,10,10,ROOT::Math::MatRepSym<double,10> > idM = ROOT::Math::SMatrixIdentity();
   if (MCSeffects) {
    invCy_10(0,0) = 0;
    for (int m = 1; m < nhits; m++) {
      for (int n = m; n < nhits; n++) {
        invCy_10(m,n) = 0;
        for (int j = 0; j < m; j++) { // j from 0 to Min[n.m] - 1
          invCy_10(m,n) += TrackFitter::sigmathetasq * (PosZ[m] - PosZ[j]) * (PosZ[n] - PosZ[j]);
        }
      } // end columns loop
    } //end of lines loop
    invCx_10 = invCy_10; //assuming sigmatheta is the same for Y and X
    invCx_10 += idM * TrackFitter::ErrorsX[0];
    invCy_10 += idM * TrackFitter::ErrorsY[0];
   } else {
      invCx_10 = idM * TrackFitter::ErrorsX[0];
      invCy_10 = idM * TrackFitter::ErrorsY[0];
   }
    invCx_10.Invert();
    invCy_10.Invert();

  xPositionsTrack_10.SetElements(TrackFitter::PosX.begin(), TrackFitter::PosX.end());
  yPositionsTrack_10.SetElements(TrackFitter::PosY.begin(), TrackFitter::PosY.end());
 }

  TVirtualFitter::SetDefaultFitter("Minuit");
  TVirtualFitter* minuit = TVirtualFitter::Fitter(0, 5);
  minuit->SetParameter(0, "X", track.getX(), 1., 0, 0);
  minuit->SetParameter(1, "Y", track.getY(), 1., 0, 0);
  minuit->SetParameter(2, "Phi", track.getPhi(), 1., 0, 0);
  minuit->SetParameter(3, "Tanl", track.getTanl(), 1., 0, 0);
  minuit->SetParameter(4, "invQPt", track.getInvQPt(), 1., 0, 0); //track.getInvQPt()*0.5

  minuit->SetFCN(myFitFcn);

  //log << "\n\n Beginning Minuit Fit with initial parameters: \n";
  //for (int i = 0; i <= 4; ++i) {
  //log << minuit->GetParName(i) << " = " << minuit->GetParameter(i)<<std::endl;
  //}

  Double_t arglist[100];
  arglist[0] = 0;
  // set print level
  minuit->ExecuteCommand("SET PRINT", arglist, 2);

  // minimize
//  arglist[0] = 10000; // number of function calls
  arglist[1] = 0.001; // tolerance
  minuit->ExecuteCommand("MIGRAD", arglist, 2);

  if (debug_minuit) {
    std::ofstream log("MinuitFitter_Log",std::ofstream::app);
    log << " Inv Q/Pt comparasion: \n";
    log << " Before fit \t\t After Fit(Minuit) \n";
    log << " " << track.getInvQPt() << "\t\t"  << minuit->GetParameter(4) << endl;
    log << " End of track comparasion\n " <<endl;
  }
  //log << "\n\n/************** Minuit Fit *************/  " << std::endl;

  //for (int i = 0; i <= 4; ++i) {
  //log << minuit->GetParName(i) << " = " << minuit->GetParameter(i) << " +/- " << minuit->GetParError(i)<<std::endl;
  //}
  double chi2, edm, errdef;
  int nvpar, nparx;
  minuit->GetStats(chi2,edm,errdef,nvpar,nparx);

  track.setTrackChi2(chi2);

  track.setX(minuit->GetParameter(0));
  track.setY(minuit->GetParameter(1));
  track.setPhi(minuit->GetParameter(2));
  track.setTanl(minuit->GetParameter(3));
  track.setInvQPt(minuit->GetParameter(4));

  SMatrix55Sym mCovariances{};
  mCovariances(0, 0) = minuit->GetParError(0)*minuit->GetParError(0);
  mCovariances(1, 1) = minuit->GetParError(1)*minuit->GetParError(1);
  mCovariances(2, 2) = minuit->GetParError(2)*minuit->GetParError(2);
  mCovariances(3, 3) = minuit->GetParError(3)*minuit->GetParError(3);
  mCovariances(4, 4) = minuit->GetParError(4)*minuit->GetParError(4);
  track.setCovariances(mCovariances);
/*
    log << "   Track covariances after MCS update: \n"
                << track.getCovariances() << std::endl
                << std::endl;
// exit(1);
*/  
//  initTrack(track,0);
  //  fit(track);

  //  log << " Kalman: " << std::endl;
  //  log << " X coord = " << track.getX() << " Y coord= " << track.getY() << " Z coord= " << track.getZ() << "Phi = " << track.getPhi() << "  Tanl= " << track.getTanl() << " InQPt =  " << track.getInvQPt() << "\n"
  //<< std::endl;
}

//_________________________________________________________________________________________________
void myFitFcn(Int_t&, Double_t*, Double_t& fval, Double_t* p, Int_t)
{
  Double_t chi2 = 0;
  Double_t tmp;
  auto zPositionsTrack = TrackFitter::PosZ;
  const int nhits = zPositionsTrack.size();

  // a will be the vector of paramters, a0, a1, a2, defined by the trajectory model
  auto& X = p[0];
  auto& Y = p[1];
  auto& phi = p[2];
  auto& tanl = p[3];
  auto& invqpt = p[4];

  ROOT::Math::SVector<double,3>  ax;
  ax[0] = X;
  ax[1] = TMath::Cos(phi) / tanl;
  ax[2] = Hz * k * invqpt * TMath::Sin(phi) / ( tanl * tanl);

  ROOT::Math::SVector<double,3>  ay;
  ay[0] = Y;
  ay[1] = TMath::Sin(phi) / tanl;
  ay[2] = -Hz * k * invqpt * TMath::Cos(phi) / ( tanl * tanl);

  if (nhits == 10) {
    chi2 += ROOT::Math::Similarity ( (xPositionsTrack_10 - G_10 * ax), invCx_10);
    chi2 += ROOT::Math::Similarity ( (yPositionsTrack_10 - G_10 * ay), invCy_10);	
  }/*
  else if (nhits == 9) {
    ROOT::Math::SMatrix<double, 9, 3> G_n3;
    for (int n = 0; n < nhits; n ++) {
      auto z = zPositionsTrack[n] - zPositionsTrack[0];
      G_n3(n,0) = 1;
      G_n3(n,1) = z;
      G_n3(n,2) = z * z * 0.5;
    } 
//    cout << " Printing 9x3 G matrix: \n" << G_n3 <<endl; 
    ROOT::Math::SMatrix<double,9,9,ROOT::Math::MatRepSym<double,9> > idM = ROOT::Math::SMatrixIdentity();
    ROOT::Math::SMatrix<double,9,9,ROOT::Math::MatRepSym<double,9> > invCy, invCx;
//    invCy = invCy * sigmaResYsq; 
   if (MCSeffects) {
    for (int m = 1; m < nhits; m++) {
      for (int n = m; n < nhits; n++) {
        for (int j = 0; j < m; j++) { // j from 0 to Min[n.m] - 1
          invCy(m,n) += TrackFitter::sigmathetasq * (zPositionsTrack[m] - zPositionsTrack[j]) * (zPositionsTrack[n] - zPositionsTrack[j]);
        }
      } // end columns loop
    } //end of lines loop

    invCx = invCy; //assuming sigmatheta is the same for Y and X
    invCx += idM * sigmaResXsq;
    invCy += idM * sigmaResYsq;
   } else {
      invCx = idM * sigmaResXsq;
      invCy = idM * sigmaResYsq;
   }
    invCx.Invert();
    invCy.Invert();

//    cout << " Printing 10x10 inv Cy matrix: \n" << invCy <<endl; 
    ROOT::Math::SVector<double, 9> xPositionsTrack(posX_array,9);
    ROOT::Math::SVector<double, 9> yPositionsTrack(posY_array,9);

    chi2 += ROOT::Math::Similarity ( (xPositionsTrack - G_n3 * ax), invCx);
    chi2 += ROOT::Math::Similarity ( (yPositionsTrack - G_n3 * ay), invCy);	
  }
  else if (nhits == 8) {
    ROOT::Math::SMatrix<double, 8, 3> G_n3;
    for (int n = 0; n < nhits; n ++) {
      auto z = zPositionsTrack[n] - zPositionsTrack[0];
      G_n3(n,0) = 1;
      G_n3(n,1) = z;
      G_n3(n,2) = z * z * 0.5;
    } 
//    cout << " Printing 8x3 G matrix: \n" << G_n3 <<endl; 
    ROOT::Math::SMatrix<double,8,8,ROOT::Math::MatRepSym<double,8> > idM = ROOT::Math::SMatrixIdentity();
    ROOT::Math::SMatrix<double,8,8,ROOT::Math::MatRepSym<double,8> > invCy, invCx;
//    invCy = invCy * sigmaResYsq; 
   if (MCSeffects) {
    for (int m = 1; m < nhits; m++) {
      for (int n = m; n < nhits; n++) {
        for (int j = 0; j < m; j++) { // j from 0 to Min[n.m] - 1
          invCy(m,n) += TrackFitter::sigmathetasq * (zPositionsTrack[m] - zPositionsTrack[j]) * (zPositionsTrack[n] - zPositionsTrack[j]);
        }
      } // end columns loop
    } //end of lines loop

    invCx = invCy; //assuming sigmatheta is the same for Y and X
    invCx += idM * sigmaResXsq;
    invCy += idM * sigmaResYsq;
   } else {
      invCx = idM * sigmaResXsq;
      invCy = idM * sigmaResYsq;
   }
    invCx.Invert();
    invCy.Invert();

//    cout << " Printing 10x10 inv Cy matrix: \n" << invCy <<endl; 
    ROOT::Math::SVector<double, 8> xPositionsTrack(posX_array,8);
    ROOT::Math::SVector<double, 8> yPositionsTrack(posY_array,8);

    chi2 += ROOT::Math::Similarity ( (xPositionsTrack - G_n3 * ax), invCx);
    chi2 += ROOT::Math::Similarity ( (yPositionsTrack - G_n3 * ay), invCy);
  }
  else if (nhits == 7) {
    ROOT::Math::SMatrix<double, 7, 3> G_n3;
    for (int n = 0; n < nhits; n ++) {
      auto z = zPositionsTrack[n] - zPositionsTrack[0];
      G_n3(n,0) = 1;
      G_n3(n,1) = z;
      G_n3(n,2) = z * z * 0.5;
    }
//cout << " Printing 7x3 G matrix: \n" << G_n3 <<endl; 
    ROOT::Math::SMatrix<double,7,7,ROOT::Math::MatRepSym<double,7> > idM = ROOT::Math::SMatrixIdentity();
    ROOT::Math::SMatrix<double,7,7,ROOT::Math::MatRepSym<double,7> > invCy, invCx;
//    invCy = invCy * sigmaResYsq; 
   if (MCSeffects) {
    for (int m = 1; m < nhits; m++) {
      for (int n = m; n < nhits; n++) {
        for (int j = 0; j < m; j++) { // j from 0 to Min[n.m] - 1
          invCy(m,n) += TrackFitter::sigmathetasq * (zPositionsTrack[m] - zPositionsTrack[j]) * (zPositionsTrack[n] - zPositionsTrack[j]);
        }
      } // end columns loop
    } //end of lines loop

    invCx = invCy;
    invCx += idM * sigmaResXsq;
    invCy += idM * sigmaResYsq;
   } else {
      invCx = idM * sigmaResXsq;
      invCy = idM * sigmaResYsq;
   }
    invCx.Invert();
    invCy.Invert();

//    cout << " Printing 10x10 inv Cy matrix: \n" << invCy <<endl; 
    ROOT::Math::SVector<double, 7> xPositionsTrack(posX_array,7);
    ROOT::Math::SVector<double, 7> yPositionsTrack(posY_array,7);

    chi2 += ROOT::Math::Similarity ( (xPositionsTrack - G_n3 * ax), invCx);
    chi2 += ROOT::Math::Similarity ( (yPositionsTrack - G_n3 * ay), invCy);
  }
*/
  fval = chi2;
}

//_________________________________________________________________________________________________
Double_t invQPtFromFCF(const FT3Track& track, Double_t bFieldZ,
                       Double_t& sigmainvqptsq)
{

  auto& xPositions = track.getXCoordinates();
  auto& yPositions = track.getYCoordinates();
  auto& zPositions = track.getZCoordinates();
  auto& SigmasX2 = track.getSigmasX2();
  auto& SigmasY2 = track.getSigmasY2();

  // Fast Circle Fit (Hansroul, Jeremie, Savard, 1987)
  auto nPoints = track.getNumberOfPoints();
  Double_t* xVal = new Double_t[nPoints];
  Double_t* yVal = new Double_t[nPoints];
  Double_t* zVal = new Double_t[nPoints];
  Double_t* xErr = new Double_t[nPoints];
  Double_t* yErr = new Double_t[nPoints];
  Double_t* uVal = new Double_t[nPoints - 1];
  Double_t* vVal = new Double_t[nPoints - 1];
  Double_t* vErr = new Double_t[nPoints - 1];
  Double_t* fweight = new Double_t[nPoints - 1];
  Double_t* Rn = new Double_t[nPoints - 1];
  Double_t* Pn = new Double_t[nPoints - 1];
  Double_t A, Aerr, B, Berr, x2, y2, invx2y2, a, b, r, sigmaRsq, u2, sigma;
  Double_t F0, F1, F2, F3, F4, SumSRn, SumSPn, SumRn, SumUPn, SumRP;

  SumSRn = SumSPn = SumRn = SumUPn = SumRP = 0.0;
  F0 = F1 = F2 = F3 = F4 = 0.0;

  for (auto np = 0; np < nPoints; np++) {
    xErr[np] = SigmasX2[np];
    yErr[np] = SigmasY2[np];
    if (np > 0) {
      xVal[np] = xPositions[np] - xVal[0];
      yVal[np] = yPositions[np] - yVal[0];
      xErr[np] *= std::sqrt(2.);
      yErr[np] *= std::sqrt(2.);
    } else {
      xVal[np] = 0.;
      yVal[np] = 0.;
    }
    zVal[np] = zPositions[np];
  }
  for (int i = 0; i < (nPoints - 1); i++) {
    x2 = xVal[i + 1] * xVal[i + 1];
    y2 = yVal[i + 1] * yVal[i + 1];
    invx2y2 = 1. / (x2 + y2);
    uVal[i] = xVal[i + 1] * invx2y2;
    vVal[i] = yVal[i + 1] * invx2y2;
    vErr[i] =
      std::sqrt(8. * xErr[i + 1] * xErr[i + 1] * x2 * y2 +
                2. * yErr[i + 1] * yErr[i + 1] * (x2 - y2) * (x2 - y2)) *
      invx2y2 * invx2y2;
    u2 = uVal[i] * uVal[i];
    fweight[i] = 1. / vErr[i];
    F0 += fweight[i];
    F1 += fweight[i] * uVal[i];
    F2 += fweight[i] * u2;
    F3 += fweight[i] * uVal[i] * u2;
    F4 += fweight[i] * u2 * u2;
  }

  Double_t Rn_det1 = F2 * F4 - F3 * F3;
  Double_t Rn_det2 = F1 * F4 - F2 * F3;
  Double_t Rn_det3 = F1 * F3 - F2 * F2;
  Double_t Pn_det1 = Rn_det2;
  Double_t Pn_det2 = F0 * F4 - F2 * F2;
  Double_t Pn_det3 = F0 * F3 - F1 * F2;

  for (int j = 0; j < (nPoints - 1); j++) {
    Rn[j] = fweight[j] *
            (Rn_det1 - uVal[j] * Rn_det2 + uVal[j] * uVal[j] * Rn_det3);
    SumSRn += Rn[j] * Rn[j] * vErr[j] * vErr[j];
    SumRn += Rn[j];

    Pn[j] = fweight[j] *
            (-Pn_det1 + uVal[j] * Pn_det2 - uVal[j] * uVal[j] * Pn_det3);
    SumSPn += Pn[j] * Pn[j] * vErr[j] * vErr[j];
    SumUPn += uVal[j] * Pn[j];

    SumRP += Rn[j] * Pn[j] * vErr[j] * vErr[j] * vErr[j];
  }

  Double_t invqpt_fcf;
  Int_t qfcf;
  //  chi2 = 0.;
  if (LinearRegression((nPoints - 1), uVal, vVal, vErr, B, Berr, A, Aerr)) {
    // v = a * u + b
    // circle passing through (0,0):
    // (x - rx)^2 + (y - ry)^2 = r^2
    // ---> a = - rx / ry;
    // ---> b = 1 / (2 * ry)
    b = 1. / (2. * A);
    a = -B * b;
    r = std::sqrt(a * a + b * b);
    double_t invR = 1. / r;

    // pt --->
    Double_t invpt = 1. / (o2::constants::math::B2C * bFieldZ * r);

    // sign(q) --->
    // rotate around the first point (0,0) to bring the last point
    // on the x axis (y = 0) and check the y sign of the rotated
    // center of the circle
    Double_t x = xVal[1], y = yVal[1],
             z = zVal[1];
    Double_t slope = TMath::ATan2(y, x);
    Double_t cosSlope = TMath::Cos(slope);
    Double_t sinSlope = TMath::Sin(slope);
    Double_t rxRot = a * cosSlope + b * sinSlope;
    Double_t ryRot = a * sinSlope - b * cosSlope;
    qfcf = (ryRot > 0.) ? -1 : +1;

    Double_t alpha = 2.0 * std::abs(TMath::ATan2(rxRot, ryRot));
    Double_t x0 = xVal[0], y0 = yVal[0], z0 = zVal[0];
    Double_t dxyz2 =
      (x - x0) * (x - x0) + (y - y0) * (y - y0) + (z - z0) * (z - z0);
    Double_t cst = 1000.;
    Double_t c_alpha = cst * alpha;
    Double_t p, pt, pz;
    pt = 1. / invpt;
    p = std::sqrt(dxyz2) / c_alpha;
    pz = std::sqrt(p * p - pt * pt);

    invqpt_fcf = qfcf * invpt;

    // error calculations:
    Double_t invA2 = 1. / (A * A);

    Double_t sigmaAsq = SumSRn / (SumRn * SumRn);
    Double_t sigmaBsq = SumSPn / (SumUPn * SumUPn);
    Double_t sigmaAB = SumRP / (SumRn * SumUPn);

    Double_t sigmaasq_FCF =
      TMath::Abs(0.25 * invA2 * invA2 *
                 (B * B * sigmaAsq + A * A * sigmaBsq - A * B * sigmaAB));
    Double_t sigmabsq_FCF = TMath::Abs(0.25 * invA2 * invA2 * sigmaAsq);
    Double_t sigma2R =
      invR * invR *
      (b * b * sigmaasq_FCF + a * a * sigmabsq_FCF +
       2 * a * b * TMath::Sqrt(sigmaasq_FCF) * TMath::Sqrt(sigmabsq_FCF));

    sigmainvqptsq = sigma2R * invpt * invpt * invR * invR;

  } else { // the linear regression failed...
    LOG(WARN) << "LinearRegression failed!";
    invqpt_fcf = 1. / 100.;
  }

  return invqpt_fcf;
}

////_________________________________________________________________________________________________
Bool_t LinearRegression(Int_t nVal, Double_t* xVal, Double_t* yVal,
                        Double_t* yErr, Double_t& B, Double_t& Berr,
                        Double_t& A, Double_t& Aerr)
{
  // linear regression y = B * x + A

  Double_t S1, SXY, SX, SY, SXX, SsXY, SsXX, SsYY, Xm, Ym, s, delta, difx;
  Double_t invYErr2;

  S1 = SXY = SX = SY = SXX = 0.0;
  SsXX = SsYY = SsXY = Xm = Ym = 0.;
  difx = 0.;
  for (Int_t i = 0; i < nVal; i++) {
    invYErr2 = 1. / (yErr[i] * yErr[i]);
    S1 += invYErr2;
    SXY += xVal[i] * yVal[i] * invYErr2;
    SX += xVal[i] * invYErr2;
    SY += yVal[i] * invYErr2;
    SXX += xVal[i] * xVal[i] * invYErr2;
    if (i > 0) {
      difx += TMath::Abs(xVal[i] - xVal[i - 1]);
    }
    Xm += xVal[i];
    Ym += yVal[i];
    SsXX += xVal[i] * xVal[i];
    SsYY += yVal[i] * yVal[i];
    SsXY += xVal[i] * yVal[i];
  }
  delta = SXX * S1 - SX * SX;
  if (delta == 0.) {
    return kFALSE;
  }
  B = (SXY * S1 - SX * SY) / delta;
  A = (SY * SXX - SX * SXY) / delta;

  Ym /= (Double_t)nVal;
  Xm /= (Double_t)nVal;
  SsYY -= (Double_t)nVal * (Ym * Ym);
  SsXX -= (Double_t)nVal * (Xm * Xm);
  SsXY -= (Double_t)nVal * (Ym * Xm);
  Double_t eps = 1.E-24;
  if ((nVal > 2) && (TMath::Abs(difx) > eps) &&
      ((SsYY - (SsXY * SsXY) / SsXX) > 0.)) {
    s = TMath::Sqrt((SsYY - (SsXY * SsXY) / SsXX) / (nVal - 2));
    Aerr = s * TMath::Sqrt(1. / (Double_t)nVal + (Xm * Xm) / SsXX);
    Berr = s / TMath::Sqrt(SsXX);
  } else {
    Aerr = 0.;
    Berr = 0.;
  }
  return kTRUE;
}

} // namespace ft3
} // namespace o2
