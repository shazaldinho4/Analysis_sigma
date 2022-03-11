#include "./DataEvent.cpp"
#include <iostream>
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TF2.h"
#include "TEllipse.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TError.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "PolFunc.h"
#include <sstream>
using namespace std;

vector<double> Phi;
vector<double> GammaP;
vector<double> ProtMom;

//BINNED FUNCTION
double MaxLikeBIN(const double *xx )
{
  double sum=0;
  const Double_t x = xx[0];
  typedef vector<int>::size_type vec_sz;
  vec_sz m=Phi.size();
  for (vec_sz i=0; i<m; i++) {
    sum+=-2.0*TMath::Log(1-GammaP[i]*x*TMath::Cos(2*Phi[i]));
  }
  return sum;
}

//UNBINNED FUNCTION
double MaxLike(const double *xx)
{
  double sum=0;
  const Double_t x1 = xx[0];
  const Double_t x2 = xx[1];
  const Double_t x3 = xx[2];
  typedef vector<int>::size_type vec_sz;

  vec_sz m=Phi.size();
  for (vec_sz i=0; i<m; i++) {
    sum+=-2.0*TMath::Log(1-GammaP[i]*(x1+x2*(ProtMom[i]*ProtMom[i])+x3*(ProtMom[i]*ProtMom[i]*ProtMom[i]))*TMath::Cos(2*Phi[i]));
  }
  return sum;
}


TF2 *fEllipse = new TF2("fEllipse","TMath::Power((x-[0])*cos([4])+(y-[2])*sin([4]),2)/TMath::Power([1],2)+TMath::Power((x-[0])*sin([4])-(y-[2])*cos([4]),2)/TMath::Power([3],2)",0,5,0,5);
TF2 *fEllipseN = new TF2("fEllipseN","TMath::Power((x-[0])*cos([4])+(y-[2])*sin([4]),2)/TMath::Power([1],2)+TMath::Power((x-[0])*sin([4])-(y-[2])*cos([4]),2)/TMath::Power([3],2)",0,5,0,5);


int myNominalCuts(DataEvent *myData);
TLorentzVector getNeutron_track(DataEvent *myData);


void Analysis(){
  gROOT->ProcessLine(".L ./Loader.C+");

  vector< vector<double> > INGammaPBIN;
  vector< vector<double> > INPhiBIN;
  INGammaPBIN.resize(10);
  INPhiBIN.resize(10);

  vector<double>  vMEASGammaP;
  vector<double>  vMEASPhi;
  vector<double>  vMEASPhiSigma;
  vector<double>  vMEASPhiProt;
  vector<double>  vPx;

  int    totalbins=10;
  double Pxbin[]={0.00, 0.03, 0.06, 0.09, 0.12, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50};
  vector<vector<double> > vvMEASPhi;
  vector<vector<double> > vvMEASPhiSigma;
  vector<vector<double> > vvMEASPhiProt;
  vector<vector<double> > vvMEASGammaP;
  vvMEASPhi.resize(totalbins);
  vvMEASPhiSigma.resize(totalbins);
  vvMEASPhiProt.resize(totalbins);

  vvMEASGammaP.resize(totalbins);

  //ELLIPSE
  double radx=0.08, rady=0.015, offsetx=0.945, offsety=1.2, angle=45*TMath::DegToRad(); //MM(k x) vs MM (k pi x)
  fEllipse->SetParameters(offsetx,radx, offsety, rady, angle);

  radx=0.008, rady=0.045, offsetx=1.195, offsety=0.944, angle=5*TMath::DegToRad(); //MM(proton) vs IM(sigma)
  fEllipseN->SetParameters(offsetx,radx, offsety, rady, angle);

  TEllipse *myEllipse=new TEllipse(offsetx, offsety, radx, rady, 0, 360, 5);
  myEllipse->SetFillStyle(0);
  myEllipse->SetLineColor(kRed);



  //FILE INPUT

  string treeName="g13b";

  //For 1 root file use constructor below
  //string fileNamePERP="/Users/shazeabayub/Documents/msc/projek/skimmedFSI/skim1_ce13_33_AUTO_v3.root";
  //DataEvent *myDataPERP=new DataEvent(fileNamePERP,treeName);

  vector<string> filename;
  vector<string>PolTableName;

  // No corresponding polarisation tables
  // filename.push_back("/Users/shazeabayub/Documents/msc/projek/skimmedFSI/skim1_ce13_33_AUTO_v3.root");//1
  // filename.push_back("/Users/shazeabayub/Documents/msc/projek/skimmedFSI/skim1_ce13_39_AUTO_v3.root");//2
  // filename.push_back("/Users/shazeabayub/Documents/msc/projek/skimmedFSI/skim1_ce13_39_PARA_v3.root");//3

  filename.push_back("/Users/shazeabayub/Documents/msc/projek/skimmedFSI/skim1_ce13_39_PERP_v3.root");//4
  filename.push_back("/Users/shazeabayub/Documents/msc/projek/skimmedFSI/skim1_ce13_42_PARA_v3.root");//5
  filename.push_back("/Users/shazeabayub/Documents/msc/projek/skimmedFSI/skim1_ce13_42_PERP_v3.root");//6
  filename.push_back("/Users/shazeabayub/Documents/msc/projek/skimmedFSI/skim1_ce15_41_PARA_v3.root");//7
  filename.push_back("/Users/shazeabayub/Documents/msc/projek/skimmedFSI/skim1_ce15_41_PERP_v3.root");//8
  filename.push_back("/Users/shazeabayub/Documents/msc/projek/skimmedFSI/skim1_ce15_45_PARA_v3.root");//9
  filename.push_back("/Users/shazeabayub/Documents/msc/projek/skimmedFSI/skim1_ce15_45_PERP_v3.root");//10
  filename.push_back("/Users/shazeabayub/Documents/msc/projek/skimmedFSI/skim1_ce17_41_PARA_v3.root");//11
  filename.push_back("/Users/shazeabayub/Documents/msc/projek/skimmedFSI/skim1_ce17_41_PERP_v3.root");//12
  filename.push_back("/Users/shazeabayub/Documents/msc/projek/skimmedFSI/skim1_ce17_47_PARA_v3.root");//13
  filename.push_back("/Users/shazeabayub/Documents/msc/projek/skimmedFSI/skim1_ce17_47_PERP_v3.root");//14
  filename.push_back("/Users/shazeabayub/Documents/msc/projek/skimmedFSI/skim1_ce19_51_PARA_v3.root");//15
  filename.push_back("/Users/shazeabayub/Documents/msc/projek/skimmedFSI/skim1_ce19_51_PERP_v3.root");//16
  filename.push_back("/Users/shazeabayub/Documents/msc/projek/skimmedFSI/skim1_ce21_51_PARA_v3.root");//17
  filename.push_back("/Users/shazeabayub/Documents/msc/projek/skimmedFSI/skim1_ce21_51_PERP_v3.root");//18
  filename.push_back("/Users/shazeabayub/Documents/msc/projek/skimmedFSI/skim1_ce21_52_PARA_v3.root");//19
  filename.push_back("/Users/shazeabayub/Documents/msc/projek/skimmedFSI/skim1_ce21_52_PERP_v3.root");//20
  filename.push_back("/Users/shazeabayub/Documents/msc/projek/skimmedFSI/skim1_ce23_52_PARA_v3.root");//21
  filename.push_back("/Users/shazeabayub/Documents/msc/projek/skimmedFSI/skim1_ce23_52_PERP_v3.root");//22

  PolTableName.push_back("/Users/shazeabayub/Documents/msc/projek/FSIAnalysis/Neutron/PolTables_ce13_39_PERP.dat");//4
  PolTableName.push_back("/Users/shazeabayub/Documents/msc/projek/FSIAnalysis/Neutron/PolTables_ce13_42_PARA.dat");//5
  PolTableName.push_back("/Users/shazeabayub/Documents/msc/projek/FSIAnalysis/Neutron/PolTables_ce13_42_PERP.dat");//6
  PolTableName.push_back("/Users/shazeabayub/Documents/msc/projek/FSIAnalysis/Neutron/PolTables_ce15_41_PARA.dat");//7
  PolTableName.push_back("/Users/shazeabayub/Documents/msc/projek/FSIAnalysis/Neutron/PolTables_ce15_41_PERP.dat");//8
  PolTableName.push_back("/Users/shazeabayub/Documents/msc/projek/FSIAnalysis/Neutron/PolTables_ce15_45_PARA.dat");//9
  PolTableName.push_back("/Users/shazeabayub/Documents/msc/projek/FSIAnalysis/Neutron/PolTables_ce15_45_PERP.dat");//10
  PolTableName.push_back("/Users/shazeabayub/Documents/msc/projek/FSIAnalysis/Neutron/PolTables_ce17_41_PARA.dat");//11
  PolTableName.push_back("/Users/shazeabayub/Documents/msc/projek/FSIAnalysis/Neutron/PolTables_ce17_41_PERP.dat");//12
  PolTableName.push_back("/Users/shazeabayub/Documents/msc/projek/FSIAnalysis/Neutron/PolTables_ce17_47_PARA.dat");//13
  PolTableName.push_back("/Users/shazeabayub/Documents/msc/projek/FSIAnalysis/Neutron/PolTables_ce17_47_PERP.dat");//14
  PolTableName.push_back("/Users/shazeabayub/Documents/msc/projek/FSIAnalysis/Neutron/PolTables_ce19_51_PARAnewTable.dat");//15
  PolTableName.push_back("/Users/shazeabayub/Documents/msc/projek/FSIAnalysis/Neutron/PolTables_ce19_51_PERPnewTable.dat");//16
  PolTableName.push_back("/Users/shazeabayub/Documents/msc/projek/FSIAnalysis/Neutron/PolTables_ce21_51_PARA.dat");//17
  PolTableName.push_back("/Users/shazeabayub/Documents/msc/projek/FSIAnalysis/Neutron/PolTables_ce21_51_PERP.dat");//18
  PolTableName.push_back("/Users/shazeabayub/Documents/msc/projek/FSIAnalysis/Neutron/PolTables_ce21_52_PARA.dat");//19
  PolTableName.push_back("/Users/shazeabayub/Documents/msc/projek/FSIAnalysis/Neutron/PolTables_ce21_52_PERP.dat");//20
  PolTableName.push_back("/Users/shazeabayub/Documents/msc/projek/FSIAnalysis/Neutron/PolTables_ce23_52_PARA.dat");//21
  PolTableName.push_back("/Users/shazeabayub/Documents/msc/projek/FSIAnalysis/Neutron/PolTables_ce23_52_PERP.dat");//22

  const int NumbOfFiles=filename.size();
  DataEvent *myDataPERP[NumbOfFiles];
  for (int i=0;i<NumbOfFiles;i++){
    myDataPERP[i]=new DataEvent(filename.at(i),treeName);
  }

  const int NumbOfPolFiles=PolTableName.size();
  for (int i=0;i<NumbOfPolFiles;i++){
    char cstr[PolTableName.at(i).size()+1];
    std::copy(PolTableName.at(i).begin(), PolTableName.at(i).end(), cstr);
    cstr[PolTableName.at(i).size()] = '\0';
    LoadPolTable(i, cstr);
  }

  TFile *Analysis_Output=new TFile("Analysis_Output.root","recreate","");

  //HISTOGRAM DEFINITIONS
  gStyle->SetTitleSize(.05, "xyz");
  gStyle->SetTitleOffset(.85, "xyz");
  gStyle->SetLabelSize(.05, "xyz");

  TH2F *h_DeltaBe[3];
  h_DeltaBe[0]=new TH2F("h_DeltaBe_0","Kaon ;p [GeV/c];#Delta #beta;",200, 0, 3, 200, -0.1, 0.1);
  h_DeltaBe[1]=new TH2F("h_DeltaBe_1","Pion ;p [GeV/c];#Delta #beta;",200, 0, 3, 200, -0.1, 0.1);
  h_DeltaBe[2]=new TH2F("h_DeltaBe_2","Neutron ;p [GeV/c];#Delta #beta;",200, 0, 3, 200, -0.1, 0.1);
  TH2F *h_BeVSp[3];
  h_BeVSp[0]=new TH2F("h_BeVSp_0","Kaon ;p [GeV/c]; #beta;",200, 0, 3, 200, 0, 1);
  h_BeVSp[1]=new TH2F("h_BeVSp_1","Pion ;p [GeV/c];#beta;",200, 0, 3, 200, 0, 1);
  h_BeVSp[2]=new TH2F("h_BeVSp_2","Neutron ;p [GeV/c]; #beta;",200, 0, 3, 200, 0, 1);
  TH1F *h_DeltaPath;
  h_DeltaPath=new TH1F("h_DeltaPath","Neutron ;#Delta path [cm]; counts;",200, -50, 50);
  TH1F *h_NeutronBeta;
  h_NeutronBeta=new TH1F("h_NeutronBeta","Neutron ;#beta; counts;",200, 0, 1.2);
  TH1F *h_NeutronBeta_cut;
  h_NeutronBeta_cut=new TH1F("h_NeutronBeta_cut","Neutron ;#beta; counts;",200, 0, 1.2);
  TH1F *h_SigmamIM;
  h_SigmamIM=new TH1F("h_SigmamIM","Invariant Mass ;IM(#pi^{-} n) [GeV]; counts;",200, 1.05, 1.35);
  TH1F *h_SigmamIM_cut;
  h_SigmamIM_cut=new TH1F("h_SigmamIM_cut","Invariant Mass ;IM(#pi^{-} n) [GeV]; counts;",200, 1.05, 1.35);
  TH2F *h_SigmamIM_NeutronBeta;
  h_SigmamIM_NeutronBeta=new TH2F("h_SigmamIM_NeutronBeta","Neutron ;#beta;IM(#pi^{-} n) [GeV];",200, 0, 1.2, 200, 1.05, 1.35);
  TH2F *h_SigmamIM_protMissMass;
  h_SigmamIM_protMissMass=new TH2F("h_SigmamIM_protMissMass","Neutron ;IM(#pi^{-} n) [GeV]; MM(#gamma d #rightarrow K^{+} #pi^{-} n) [GeV];",200, 1.05, 1.35, 200, 0.7, 1.2);
  TH2F *h_SigmamIM_protMissMass_cut;
  h_SigmamIM_protMissMass_cut=new TH2F("h_SigmamIM_protMissMass_cut","Neutron ;IM(#pi^{-} n) [GeV]; MM(#gamma d #rightarrow K^{+} #pi^{-} n) [GeV];",200, 1.05, 1.35, 200, 0.7, 1.2);
  TH1F *h_MissingMomentum;
  h_MissingMomentum=new TH1F("h_MissingMomentum","Spectator Proton ;p_{X} [GeV]; counts;",200, 0, 1.0);
  TH1F *h_MissingMomentum_cut;
  h_MissingMomentum_cut=new TH1F("h_MissingMomentum_cut","Spectator Proton ;p_{X} [GeV]; counts;",200, 0, 1.0);
  TH1F *h_photonselection_pi;
  h_photonselection_pi=new TH1F("h_photonselection_pi","Pion ;#Delta t [ns]; counts;",1000, -5.0, 5.0);
  TH1F *h_photonselection_k;
  h_photonselection_k=new TH1F("h_photonselection_k","Kaon ;#Delta t [ns]; counts;",1000, -5.0, 5.0);
  TH2F *h_misID;
  h_misID=new TH2F("h_misID","Misidentification ; MM(#gamma n #rightarrow K^{+} #pi^{-} X) [GeV] ; MM(#gamma n #rightarrow #pi^{+} #pi^{-} X) [GeV];",600,0.5,1.4,600,0.5,1.4);
  TH2F*h_ellipse;
  h_ellipse=new TH2F("h_ellipse","Background Removal ; MM(#gamma n #rightarrow K^{+} #pi^{-} X) [GeV] ; MM(#gamma n #rightarrow K^{+} X) [GeV];",600,0.6,1.4,600,0.6,1.4);

  int testing=0;

  for (int i=0;i<NumbOfFiles;i++){
    cout<<"File has "<<myDataPERP[i]->getEntries()<<" entries."<<endl;
    while (myDataPERP[i]->getEntry()<myDataPERP[i]->getEntries()){
      if (testing && myDataPERP[i]->getEntry()>10000)
      break;
      myDataPERP[i]->getNextEntry();
      if (myDataPERP[i]->getEntry() % 10000 == 0){
        fprintf (stderr, "Loop over [%s] %.2f percent\r", filename.at(i).c_str(), myDataPERP[i]->getEntry()*100.0/myDataPERP[i]->getEntries());
        fflush (stderr);
      }

      double deltbeta[3]; //define deltbeta as array with 3 elements??
      if (myNominalCuts(myDataPERP[i])!=1)
      continue;

      for (int ij=0;ij<myDataPERP[i]->getNum_tracks();ij++){
        deltbeta[ij]=myDataPERP[i]->getEVNT_track(ij).Beta()-myDataPERP[i]->getEVNT_bem(ij); //Beta()= p/E   bem= Beta from timeof flight
        h_DeltaBe[ij]->Fill(myDataPERP[i]->getEVNT_track(ij).Rho(),deltbeta[ij]);
        h_BeVSp[ij]->Fill(myDataPERP[i]->getEVNT_track(ij).Rho(),myDataPERP[i]->getEVNT_bem(ij));
        //this block is simply calculating the delta beta by calling each beta and subtracting. it then fills the beta vs momentum plot
      }

      TVector3 mvrt=myDataPERP[i]->getEVNT_vertex(2); //gets event vertex for a neutron
      double calcPath=TMath::Sqrt(TMath::Power(myDataPERP[i]->getECPB_X(2)-mvrt.X(),2)+TMath::Power(myDataPERP[i]->getECPB_Y(2)-mvrt.Y(),2)+TMath::Power(myDataPERP[i]->getECPB_Z(2)-mvrt.Z(),2));//calculated path of the neutron
      h_DeltaPath->Fill(calcPath-myDataPERP[i]->getECPB_path(0));//path difference. also fills the histogram
      h_NeutronBeta->Fill(myDataPERP[i]->getECPB_betacorr(2));

      TLorentzVector neutron=getNeutron_track(myDataPERP[i]);
      TLorentzVector p4PiN=myDataPERP[i]->getEVNT_track(1)+neutron;//1 is pion
      h_SigmamIM->Fill(p4PiN.M());
      h_SigmamIM_NeutronBeta->Fill(myDataPERP[i]->getECPB_betacorr(2),p4PiN.M()); //what is this betacorr?//correction factor within the calorimeter

      TLorentzVector photonbeam, deuteron;
      photonbeam.SetXYZM(0,0,myDataPERP[i]->getTAGR_epho(myDataPERP[i]->getIndex_pi(0)),0); //sets lorentz vector with only a z component
      deuteron.SetXYZM(0,0,0,1.8756);//deutron started at rest so only has M component
      TLorentzVector MissingMass=photonbeam+deuteron-myDataPERP[i]->getEVNT_track(0)-myDataPERP[i]->getEVNT_track(1)-neutron;// missing proton mass
      h_SigmamIM_protMissMass->Fill(p4PiN.M(),MissingMass.M());
      h_MissingMomentum->Fill(MissingMass.Rho());
      if (MissingMass.M()>0.89 && MissingMass.M()<1.0)
      h_SigmamIM_cut->Fill(p4PiN.M());

      h_misID->Fill(myDataPERP[i]->getKpi_mm(myDataPERP[i]->getIndex_pi(0)).M(),myDataPERP[i]->getPipi_mm(myDataPERP[i]->getIndex_pi(0)).M());
      h_ellipse->Fill(myDataPERP[i]->getKpi_mm(myDataPERP[i]->getIndex_pi(0)).M(),myDataPERP[i]->getK_mm(myDataPERP[i]->getIndex_pi(0)).M());
      h_photonselection_pi->Fill(myDataPERP[i]->getDelt_t_pi(myDataPERP[i]->getIndex_pi(0)));
      h_photonselection_k->Fill(myDataPERP[i]->getDelt_t_k(myDataPERP[i]->getIndex_k(0)));

      //Neutron Virtuality Definition
      TLorentzVector mandt=photonbeam-myDataPERP[i]->getEVNT_track(0)-p4PiN;
      double NeutronVirtuality=0.9395654*0.9395654-mandt.M2();

      TLorentzVector neutron_target;
      neutron_target.SetXYZM(0,0,0,0.9395654);
      TLorentzVector PhotonNeutron = photonbeam+neutron_target;
      TVector3 PhotonNeutron_boost = PhotonNeutron.BoostVector();
      TLorentzVector Kaon_copy = myDataPERP[i]->getEVNT_track(0);
      Kaon_copy.Boost(-PhotonNeutron_boost);

      //TCUT
      // if (Kaon_copy.CosTheta()>0)
      // continue;

      if (fEllipseN->Eval(p4PiN.M(),MissingMass.M())>1)
      continue;
      h_SigmamIM_protMissMass_cut->Fill(p4PiN.M(),MissingMass.M());
      h_NeutronBeta_cut->Fill(myDataPERP[i]->getECPB_betacorr(2));
      h_MissingMomentum_cut->Fill(MissingMass.Rho());

      int ik=i;
      double PhotoPol=GetPol(ik, myDataPERP[i]->getCoh_edge(), myDataPERP[i]->getTAGR_epho(myDataPERP[i]->getIndex_pi(0))*1000.0, 8, 0.2,0.3);
      if (PhotoPol<0.5) continue;

      for(Int_t loc_it=0;loc_it<totalbins;loc_it++){
        if (MissingMass.Rho()>Pxbin[loc_it]//missingmass->NeutronVirtuality
        &&MissingMass.Rho()<Pxbin[loc_it+1]){//if its within the bin width //missingmass->NeutronVirtuality
          if (myDataPERP[i]->getCoh_plan()==0){
            vvMEASGammaP[loc_it].push_back(PhotoPol); //polarisation
            vvMEASPhi[loc_it].push_back(myDataPERP[i]->getEVNT_track(0).Phi());
            vMEASGammaP.push_back(PhotoPol);
            vMEASPhi.push_back(myDataPERP[i]->getEVNT_track(0).Phi());
            vvMEASPhiSigma[loc_it].push_back(p4PiN.Phi());
            vMEASPhiSigma.push_back(p4PiN.Phi());
            vvMEASPhiProt[loc_it].push_back(MissingMass.Phi());
            vMEASPhiProt.push_back(MissingMass.Phi());
            vPx.push_back(MissingMass.Rho()); //missingmass->NeutronVirtuality
            break;
          }
          else if (myDataPERP[i]->getCoh_plan()==1){
            vvMEASGammaP[loc_it].push_back(-PhotoPol);
            vvMEASPhi[loc_it].push_back(myDataPERP[i]->getEVNT_track(0).Phi());
            vMEASGammaP.push_back(-PhotoPol);
            vMEASPhi.push_back(myDataPERP[i]->getEVNT_track(0).Phi());
            vvMEASPhiSigma[loc_it].push_back(p4PiN.Phi());
            vMEASPhiSigma.push_back(p4PiN.Phi());
            vvMEASPhiProt[loc_it].push_back(MissingMass.Phi());
            vMEASPhiProt.push_back(MissingMass.Phi());
            vPx.push_back(MissingMass.Rho()); //missingmass->NeutronVirtuality
            //same block as the last but for other polarisation
            break;
          }

        }
      }
    }
  }
  cout<<endl;

  TCanvas *c0=new TCanvas("c0","My plots", 900, 500);
  c0->Divide(2,3);
  c0->cd(1);
  h_DeltaBe[0]->Draw("colz");
  c0->cd(3);
  h_DeltaBe[1]->Draw("colz");
  c0->cd(5);
  h_DeltaBe[2]->Draw("colz");
  c0->cd(2);
  h_BeVSp[0]->Draw("colz");
  c0->cd(4);
  h_BeVSp[1]->Draw("colz");
  c0->cd(6);
  h_BeVSp[2]->Draw("colz");

  TCanvas *c1=new TCanvas("c1","My plots", 900, 500);
  c1->Divide(2,3);
  c1->cd(1);
  h_DeltaPath->Draw("");
  c1->cd(2);
  h_NeutronBeta->Draw("");
  h_NeutronBeta_cut->SetLineColor(kRed);
  h_NeutronBeta_cut->Draw("same");
  c1->cd(3);
  h_SigmamIM->Draw("");
  h_SigmamIM_cut->SetLineColor(kRed);
  h_SigmamIM_cut->Draw("same");
  c1->cd(4);
  h_SigmamIM_NeutronBeta->Draw("colz");
  c1->cd(5);
  h_SigmamIM_protMissMass->Draw("colz");
  myEllipse->Draw("same");
  c1->cd(6);
  h_MissingMomentum->Draw("");
  h_MissingMomentum_cut->SetLineColor(kRed);
  h_MissingMomentum_cut->Draw("same");

  double SigmaArr[10]={0};
  double eSigma[10]={0};
  double SigmaArrSigma[10]={0};
  double eSigmaSigma[10]={0};
  double SigmaArrProt[10]={0};
  double eSigmaProt[10]={0};
  double eX[10]={0};
  double X[10]={0};
  const char *minName = "Minuit2";
  const char *algoName = "";

  //BINNED TECHNIQUE
  for (int ik=0; ik<totalbins; ik++){
    ROOT::Math::Minimizer* minimum = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
    minimum->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
    minimum->SetMaxIterations(10000);  // for GSL
    minimum->SetTolerance(0.001);
    minimum->SetPrintLevel(1);

    //binned kaon
    Phi=vvMEASPhi.at(ik);
    GammaP=vvMEASGammaP.at(ik);
    ROOT::Math::Functor f(&MaxLikeBIN,1);
    double step[1] = {0.01};
    double variable[1] = {0.2};
    minimum->SetFunction(f);
    minimum->SetVariable(0,"Sigma",variable[0], step[0]);
    minimum->Minimize();
    const double *xs = minimum->X();
    const double *xes= minimum->Errors();
    SigmaArr[ik]=xs[0];
    eSigma[ik]=xes[0];

    //binned sigma
    ROOT::Math::Functor fSigma(&MaxLikeBIN,1);
    minimum->SetFunction(fSigma);
    Phi=vvMEASPhiSigma.at(ik);
    GammaP=vvMEASGammaP.at(ik);
    minimum->SetVariable(0,"Sigma",variable[0], step[0]);
    minimum->Minimize();
    xs = minimum->X();
    xes= minimum->Errors();
    SigmaArrSigma[ik]=xs[0];
    eSigmaSigma[ik]=xes[0];

    //binned proton
    ROOT::Math::Functor fProton(&MaxLikeBIN,1);
    minimum->SetFunction(fProton);
    Phi=vvMEASPhiProt.at(ik);
    GammaP=vvMEASGammaP.at(ik);
    minimum->SetVariable(0,"Sigma",variable[0], step[0]);
    minimum->Minimize();
    xs = minimum->X();
    xes= minimum->Errors();
    SigmaArrProt[ik]=xs[0];
    eSigmaProt[ik]=xes[0];

    X[ik]=(Pxbin[ik+1]+Pxbin[ik])/2.0;

  }

  //UNBINNED TECHNIQUE

  //unbinned kaon
  ROOT::Math::Minimizer* minimum = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
  minimum->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
  minimum->SetMaxIterations(10000);  // for GSL
  minimum->SetTolerance(0.001);
  minimum->SetPrintLevel(1);
  Phi=vMEASPhi;//
  GammaP=vMEASGammaP;
  ProtMom=vPx;
  ROOT::Math::Functor f(&MaxLike,3);
  double step[3] = {0.01,0.01,0.01};//steps to find minimisation
  double variable[3] = {0.1,0.1,0.1};//where to start minimisation
  minimum->SetFunction(f);
  minimum->SetVariable(0,"x0",variable[0], step[0]);
  minimum->SetVariable(1,"x1",variable[1], step[1]);
  minimum->SetVariable(2,"x2",variable[2], step[2]);
  minimum->Minimize();
  const double *xs = minimum->X();
  const double *xes = minimum->Errors();
  TF1 *SigmaFunc3=new TF1("SigmaFunc3","[0]+[1]*(x*x)+[2]*(x*x*x)",0, 0.5);
  SigmaFunc3->SetParameters(xs[0],xs[1],xs[2]);
  double mydetValuesKaon[3]={xs[0],xs[1],xs[2]};

  //unbinned sigma
  Phi=vMEASPhiSigma;
  GammaP=vMEASGammaP;
  ProtMom=vPx;
  ROOT::Math::Functor fSigma(&MaxLike,3);
  minimum->SetFunction(fSigma);
  minimum->SetVariable(0,"x0",variable[0], step[0]);
  minimum->SetVariable(1,"x1",variable[1], step[1]);
  minimum->SetVariable(2,"x2",variable[2], step[2]);
  minimum->Minimize();
  xs = minimum->X();
  xes = minimum->Errors();
  TF1 *SigmaFunc3Sigma=new TF1("SigmaFunc3Sigma","[0]+[1]*(x*x)+[2]*(x*x*x)",0, 0.5);
  SigmaFunc3Sigma->SetParameters(xs[0],xs[1],xs[2]);
  double mydetValuesSigma[3]={xs[0],xs[1],xs[2]};

  //unbinned proton
  Phi=vMEASPhiProt;
  GammaP=vMEASGammaP;
  ProtMom=vPx;
  ROOT::Math::Functor fProt(&MaxLike,3);
  minimum->SetFunction(fProt);
  minimum->SetVariable(0,"x0",variable[0], step[0]);
  minimum->SetVariable(1,"x1",variable[1], step[1]);
  minimum->SetVariable(2,"x2",variable[2], step[2]);
  minimum->Minimize();
  xs = minimum->X();
  xes = minimum->Errors();
  TF1 *SigmaFunc3Prot=new TF1("SigmaFunc3Prot","[0]+[1]*(x*x)+[2]*(x*x*x)",0, 0.5);
  SigmaFunc3Prot->SetParameters(xs[0],xs[1],xs[2]);
  double mydetValuesProton[3]={xs[0],xs[1],xs[2]};

  //BOOTSTRAPPING
  TH2D *Extracted=new TH2D("Extracted", "Kaon",100, 0,0.5, 200, -1,1);
  TH2D *ExtractedSigma=new TH2D("ExtractedSigma", "Sigma",100, 0,0.5, 200, -1,1);
  TH2D *ExtractedProt=new TH2D("ExtractedProton", "Proton",100, 0,0.5, 200, -1,1);
  Extracted->SetTitle("Kaon ; Missing Momentum (p) [GeV] ; #Sigma");
  ExtractedSigma->SetTitle("Sigma ; Missing Momentum (p) [GeV] ; #Sigma");
  ExtractedProt->SetTitle("Proton ; Missing Momentum (p) [GeV] ; #Sigma");
  Extracted->GetYaxis()->SetRangeUser(-1,1);
  ExtractedSigma->GetYaxis()->SetRangeUser(-1,1);
  ExtractedProt->GetYaxis()->SetRangeUser(-1,1);

  TRandom3 *rndmn=new TRandom3(0);

  for (int ik=0; ik<200; ik++){
    vector<double> INGammaP2,INPhi2,INVpx2;
    vector<double> INPhiSigma2;
    vector<double> INPhiProton2;
    for (int i=0;i<vMEASPhi.size();i++){
      int k=rndmn->Integer(vMEASPhi.size());
      INGammaP2.push_back(vMEASGammaP[k]);
      INPhi2.push_back(vMEASPhi[k]);
      INVpx2.push_back(vPx[k]);

      INPhiSigma2.push_back(vMEASPhiSigma[k]);
      INPhiProton2.push_back(vMEASPhiProt[k]);
    }

    //bootstrapping kaon
    ROOT::Math::Minimizer* minimum = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
    minimum->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
    minimum->SetMaxIterations(10000);  // for GSL
    minimum->SetTolerance(0.001);
    minimum->SetPrintLevel(1);
    Phi=INPhi2;
    GammaP=INGammaP2;
    ProtMom=INVpx2;
    ROOT::Math::Functor f(&MaxLike,3);
    double step[3] = {0.01,0.01,0.01};
    double variable[3] = {0.1,0.1,0.1};
    minimum->SetFunction(f);
    minimum->SetVariable(0,"x0",variable[0], step[0]);
    minimum->SetVariable(1,"x1",variable[1], step[1]);
    minimum->SetVariable(2,"x2",variable[2], step[2]);
    minimum->Minimize();
    const double *xs = minimum->X();
    const double *xes = minimum->Errors();
    TF1 *SigmaFunc3_BS=new TF1("SigmaFunc3_BS","[0]+[1]*(x*x)+[2]*(x*x*x)",0, 0.5);
    SigmaFunc3_BS->SetParameters(xs[0],xs[1],xs[2]);
    SigmaFunc3_BS->SetParLimits(1,0,5);

    //bootstrapping sigma
    Phi=INPhiSigma2;
    GammaP=INGammaP2;
    ProtMom=INVpx2;
    ROOT::Math::Functor fSigma(&MaxLike,3);
    minimum->SetFunction(fSigma);
    minimum->SetVariable(0,"x0",variable[0], step[0]);
    minimum->SetVariable(1,"x1",variable[1], step[1]);
    minimum->SetVariable(2,"x2",variable[2], step[2]);
    minimum->Minimize();
    xs = minimum->X();
    xes = minimum->Errors();
    TF1 *SigmaFunc3_BSSigma=new TF1("SigmaFunc3_BSSigma","[0]+[1]*(x*x)+[2]*(x*x*x)",0, 0.5);
    SigmaFunc3_BSSigma->SetParameters(xs[0],xs[1],xs[2]);
    SigmaFunc3_BSSigma->SetParLimits(1,0,5);

    //bootstrapping proton
    Phi=INPhiProton2;
    GammaP=INGammaP2;
    ProtMom=INVpx2;
    ROOT::Math::Functor fProt(&MaxLike,3);
    minimum->SetFunction(fProt);
    minimum->SetVariable(0,"x0",variable[0], step[0]);
    minimum->SetVariable(1,"x1",variable[1], step[1]);
    minimum->SetVariable(2,"x2",variable[2], step[2]);
    minimum->Minimize();
    xs = minimum->X();
    xes = minimum->Errors();
    TF1 *SigmaFunc3_BSProt=new TF1("SigmaFunc3_BSProt","[0]+[1]*(x*x)+[2]*(x*x*x)",0, 0.5);
    SigmaFunc3_BSProt->SetParameters(xs[0],xs[1],xs[2]);


    //bootstrapping histograms
    for (int kk=0; kk<Extracted->GetNbinsX(); kk++){
      TAxis *xaxis = Extracted->GetXaxis();
      Double_t binCenter = xaxis->GetBinCenter(kk);
      Extracted->Fill(binCenter, SigmaFunc3_BS->Eval(binCenter));
    }
    for (int kk=0; kk<ExtractedSigma->GetNbinsX(); kk++){
      TAxis *xaxissigma = ExtractedSigma->GetXaxis();
      Double_t binCentersigma = xaxissigma->GetBinCenter(kk);
      ExtractedSigma->Fill(binCentersigma, SigmaFunc3_BSSigma->Eval(binCentersigma));
    }
    for (int kk=0; kk<ExtractedProt->GetNbinsX(); kk++){
      TAxis *xaxisprot = ExtractedProt->GetXaxis();
      Double_t binCenterprot = xaxisprot->GetBinCenter(kk);
      ExtractedProt->Fill(binCenterprot, SigmaFunc3_BSProt->Eval(binCenterprot));
    }




  }


  //DRAWING
  //COMPARISON TO NICK'S ANALYSIS
  double_t xvalue[]= {0};
  double_t exvalue[]= {0};

  double_t oneone_pos[] = {0.807194702};
  double_t oneone_pos_uncr[] = {0.006589945};

  double_t onethree_pos[] = {0.896037586};
  double_t onethree_pos_uncr[] = {0.003370546};

  double_t onefive_pos[] = {0.822728768};
  double_t onefive_pos_uncr[] = {0.00349643};

  double_t oneseven_pos[] = {0.924893477};
  double_t oneseven_pos_uncr[] = {0.003507759};

  double_t onenine_pos[] = {0.899835906};
  double_t onenine_pos_uncr[] = {0.004911022};

  double_t twoone_pos[] = {0.88526966};
  double_t twoone_pos_uncr[] = {0.005447865};

  double_t oneone_neg[] = {0.762523916};
  double_t oneone_neg_uncr[] = {0.009325124};

  double_t onethree_neg[] = {0.875788791};
  double_t onethree_neg_uncr[] = {0.004313865};

  double_t onefive_neg[] = {0.818716011};
  double_t onefive_neg_uncr[] = {0.004936524};

  double_t oneseven_neg[] = {0.841185348};
  double_t oneseven_neg_uncr[] = {0.006330255};

  double_t onenine_neg[] = {0.785482969};
  double_t onenine_neg_uncr[] = {0.015944021};

  double_t twoone_neg[] = {0.646815886};
  double_t twoone_neg_uncr[] = {0.030527737};

  // TGraphErrors *oneone_posg=new TGraphErrors(1.0, xvalue, onefive_pos,exvalue, onefive_pos_uncr);//x,y,ex,ey
  // oneone_posg->SetLineWidth(5);
  // oneone_posg->SetLineColor(kMagenta);

  TGraphErrors *gsigma=new TGraphErrors(totalbins, X, SigmaArr, eX, eSigma);
  gsigma->SetTitle("Kˆ{+} ;Missing Momentum(p) [GeV] ;#Sigma");
  TGraphErrors *gsigmaSigma=new TGraphErrors(totalbins, X, SigmaArrSigma, eX, eSigmaSigma);
  gsigmaSigma->SetTitle("#Sigmaˆ{-} ;Missing Momentum(p) [GeV] ;#Sigma");
  TGraphErrors *gsigmaProt=new TGraphErrors(totalbins, X, SigmaArrProt, eX, eSigmaProt);
  gsigmaProt->SetTitle("Proton ;Missing Momentum(p) [GeV] ;#Sigma");
  TCanvas *ci=new TCanvas("ci","Malakia", 1400,700);
  gsigma->GetYaxis()->SetRangeUser(-1,1);
  gsigmaSigma->GetYaxis()->SetRangeUser(-1,1);
  gsigmaProt->GetYaxis()->SetRangeUser(-1,1);
  gsigma->GetXaxis()->SetRangeUser(0,.5);
  gsigmaSigma->GetXaxis()->SetRangeUser(0,.5);
  gsigmaProt->GetXaxis()->SetRangeUser(0,.5);

  // ci->Divide(2,2);
  // ci->cd(1);
  // gsigma->Draw("PA");
  SigmaFunc3->SetLineColor(kBlue);
  gsigma->Fit("pol2");
  // Extracted->Draw("colz");
  // SigmaFunc3->Draw("A");
  // // oneone_posg->Draw("same");
  // ci->cd(2);
  // gsigmaSigma->Draw("PA");
  SigmaFunc3Sigma->SetLineColor(kBlue);
  gsigmaSigma->Fit("pol2");
  // ExtractedSigma->Draw("colz");
  // SigmaFunc3Sigma->Draw("A");
  // // oneone_posg->Draw("same");
  // ci->cd(3);
  // gsigmaProt->Draw("PA");
  SigmaFunc3Prot->SetLineColor(kBlue);
  gsigmaProt->Fit("pol2");
  // ExtractedProt->Draw("colz");
  // SigmaFunc3Prot->Draw("A");
  // oneone_posg->Draw("same");

  // SigmaFunc3->SetLineColor(kBlue);
  // SigmaFunc3Sigma->SetLineColor(kRed);
  // SigmaFunc3Prot->SetLineColor(kGreen);
  // SigmaFunc3->Draw("");
  // SigmaFunc3Sigma->Draw("same");
  // SigmaFunc3Prot->Draw("same");


  Analysis_Output->Write();
  gsigma->SetName("KaonBSA");
  gsigmaSigma->SetName("SigmaBSA");
  gsigmaProt->SetName("ProtonBSA");

  SigmaFunc3->SetName("KaonUnbinned");
  SigmaFunc3Sigma->SetName("SigmaUnbinned");
  SigmaFunc3Prot->SetName("ProtonUnbinned");

  SigmaFunc3->Write();
  SigmaFunc3Sigma->Write();
  SigmaFunc3Prot->Write();
  gsigma->Write();
  gsigmaSigma->Write();
  gsigmaProt->Write();
  //oneone_posg->Write();

  //   cout<<"Kaon: "<<mydetValuesKaon[0]<<" "<<mydetValuesKaon[1]<<" "<<mydetValuesKaon[2]<<endl;
  //   cout<<"Sigma: "<<mydetValuesSigma[0]<<" "<<mydetValuesSigma[1]<<" "<<mydetValuesSigma[2]<<endl;
  //   cout<<"Proton: "<<mydetValuesProton[0]<<" "<<mydetValuesProton[1]<<" "<<mydetValuesProton[2]<<endl;
  //
  //   cout<<"Kaon: "<<SigmaArr[0]<<" "<<SigmaArr[1]<<" "<<SigmaArr[2]<<" "<<SigmaArr[3]<<" "<<SigmaArr[4]<<" "<<SigmaArr[5]<<" "<<SigmaArr[6]<<" "<<SigmaArr[7]<<" "<<SigmaArr[8]<<" "<<SigmaArr[9]<<endl;
  //   cout<<"Sigma: "<<SigmaArrSigma[0]<<" "<<SigmaArrSigma[1]<<" "<<SigmaArrSigma[2]<<" "<<SigmaArrSigma[3]<<" "<<SigmaArrSigma[4]<<" "<<SigmaArrSigma[5]<<" "<<SigmaArrSigma[6]<<" "<<SigmaArrSigma[7]<<" "<<SigmaArrSigma[8]<<" "<<SigmaArrSigma[9]<<endl;
  //   cout<<"Proton: "<<SigmaArrProt[0]<<" "<<SigmaArrProt[1]<<" "<<SigmaArrProt[2]<<" "<<SigmaArrProt[3]<<" "<<SigmaArrProt[4]<<" "<<SigmaArrProt[5]<<" "<<SigmaArrProt[6]<<" "<<SigmaArrProt[7]<<" "<<SigmaArrProt[8]<<" "<<SigmaArrProt[9]<<endl;
}

//CUTS
int myNominalCuts(DataEvent *myDataPERP){
  if(myDataPERP->getNum_tracks()!=3)
  return 0;
  if(fabs(myDataPERP->getDelt_t_pi(myDataPERP->getIndex_pi(0)))>1.0)//coincidence time pion
  return 0;
  if(fabs(myDataPERP->getDelt_t_k(myDataPERP->getIndex_k(0)))>1.0)//coincidence time kaon
  return 0;
  if (fabs(myDataPERP->getEVNT_track(0).Beta()-myDataPERP->getEVNT_bem(0))>0.030)//KAON PID - sys
  return 0;
  if (fabs(myDataPERP->getEVNT_track(1).Beta()-myDataPERP->getEVNT_bem(1))>0.050)//PION PID - sys
  return 0;
  if (myDataPERP->getEVNT_vertex(0).Z()<-40 &&myDataPERP->getEVNT_vertex(0).Z()>0)//kaon event vertex
  return 0;
  if (myDataPERP->getPipi_mm(myDataPERP->getIndex_pi(0)).M()<1.015)//remove misidentified kaons - sys
  return 0;
  if (myDataPERP->getTAGR_epho(myDataPERP->getIndex_pi(0))>myDataPERP->getCoh_edge()/1000.0)//remove photons that have energies above coherent edge - sys
  return 0;
  if (myDataPERP->getTAGR_epho(myDataPERP->getIndex_pi(0))<myDataPERP->getCoh_edge()/1000.0-0.200)//remove photons that have energies above coherent edge - sys
  return 0;

  // if (myDataPERP->getTAGR_epho(myDataPERP->getIndex_pi(0))<2.1 || myDataPERP->getTAGR_epho(myDataPERP->getIndex_pi(0))>2.3)
  // return 0;
  //ECUT
  return 1;
}

//NEUTRON 4-VECTOR
//calculated from beta
TLorentzVector getNeutron_track(DataEvent *myData){
  TLorentzVector neutron;
  double beta=myData->getECPB_betacorr(2);
  double denominator=0;
  if (TMath::Power(beta,2)<1){
    denominator=TMath::Power(beta,2);
  }
  else denominator=0.99999;

  float momentum=beta*0.9395654/TMath::Sqrt(1-denominator);
  float cx=myData->getEVNT_track(2).X()/myData->getEVNT_track(2).Rho();
  float cy=myData->getEVNT_track(2).Y()/myData->getEVNT_track(2).Rho();
  float cz=myData->getEVNT_track(2).Z()/myData->getEVNT_track(2).Rho();
  neutron.SetXYZM(momentum*cx,momentum*cy,momentum*cz,0.9395654);
  return neutron;
}


/* DESCRIPTION OF METHODS
There are three particles in each event. The first particle is a proton, the second is a positive kaon and the third is a negative pion
You can access the following info
TLorentzVector getEVNT_track(int i){return loc_EVNT_track->at(i);} //TLorentzVector for kaon (i=0) pion (i=1) and neutron (i=2). The TlorentzVector has the nominal masses of the tracks.
int getEVNT_q(int i){return loc_EVNT_q->at(i);} //charge for track i
int getEVNT_scsec(int i){return loc_EVNT_scsec->at(i);} //SC sector for track i
int getEVNT_scpad(int i){return loc_EVNT_scpad->at(i);} //SC paddle for track i
int getEVNT_schit(int i){return loc_EVNT_schit->at(i);} //SC hit for track i
int getEVNT_stsec(int i){return loc_EVNT_stsec->at(i);} //ST sector for track i
int getEVNT_sthit(int i){return loc_EVNT_sthit->at(i);} //ST hit for track i
int getTAGR_eid(int i){return loc_tagr_eid->at(i);} //TAGR eid for photon i
int getTAGR_tid(int i){return loc_tagr_tid->at(i);} //TAGR tid for photon i
int getTAGR_stat(int i){return loc_tagr_stat->at(i);} //TAGR status for photon i
int getSTPB_sthid(int i){return loc_STPB_sthid->at(i);} //STBP hitd for track i
int getSCPB_ScPdHt(int i){return loc_SCPB_ScPdHt->at(i);} //SCBP ScPdHt for track i
float getEVNT_bem(int i){return loc_EVNT_bem->at(i);} //beta measured for track i
float getEVNT_sc_t(int i){return loc_EVNT_sc_t->at(i);} //sc time for track i
float getEVNT_sc_d(int i){return loc_EVNT_sc_d->at(i);} //sc d for track i
float getEVNT_st_t(int i){return loc_EVNT_st_t->at(i);} //st time for track i
float getEVNT_st_d(int i){return loc_EVNT_st_d->at(i);} //st d for track i
float getEVNT_sc_e(int i){return loc_EVNT_sc_e->at(i);} //sc energy for track i
float getTAGR_epho(int i){return loc_tagr_epho->at(i);} //TAGR epho for photon i
float getTAGR_tpho(int i){return loc_tagr_tpho->at(i);} //TAGR tpho for photon i
int getNum_photons(){return loc_tagr_tpho->size();} //number of photons
int getTrip_flag(){return loc_trip_flag;} //trip flag
int getNumofpart(){return loc_numofpart;} //number of particles (including neutrals)
int getNum_pos(){return loc_num_pos;} //number of positive
int getNum_chargedtracks(){return loc_EVNT_track->size();} //number of charged
int getNum_neg(){return loc_num_neg;} //number of negative
int getNum_neu(){return loc_num_neu;} //number of neutrals
int getHEAD_eventnum(){return loc_head_eventnum;} //HEAD event number
int getHEAD_runnum(){return loc_head_runnum;} //HEAD run number
int getNum_deuterons(){return loc_num_deuterons;} //number of deuterons
int getNum_protons(){return loc_num_protons;} //number of protons
int getNum_poskaons(){return loc_num_poskaons;} //number of postive kaons
int getNum_pospions(){return loc_num_pospions;} //number of positive pions
int getNum_negkaons(){return loc_num_negkaons;} //number of negative kaons
int getNum_negpions(){return loc_num_negpions;} ////number of negative pions
float getCoh_edge(){return loc_coh_edge;} //Coherent Edge from EPICS
float getBeam_en(){return loc_beam_en;} //Beam energy from EPICS
float getCoh_edge_nom(){return loc_coh_edge_nom;} //Nominal Coherent Edge
int getCoh_plan_db(){return loc_coh_plan_db;} //Coherent plane from database
int getCoh_radi(){return loc_coh_radi;} //Coherent radiator from EPICS
int getCoh_plan(){return loc_coh_plan;} //Coherent radiator from EPICS
float getDelt_t_k(int i){return loc_delt_t_k->at(i);} //Photon i coincidence time with kaon
float getDelt_t_pi(int i){return loc_delt_t_pi->at(i);} //Photon i coincidence time with pion
int getNumph_k(){return loc_numph_k;} //Number of photons within 1ns when looking at the photon-kaon coincidence time
int getNumph_pi(){return loc_numph_pi;} //Number of photons within 1ns when looking at the photon-pion coincidence time
int getIndex_k(int i){return loc_index_k->at(i);} //Photon index when sorted using the photon-kaon coincidence time
int getIndex_pi(int ip){return loc_index_pi->at(i);} //Photon index when sorted using the photon-pion coincidence time. getIndex_pi(0) returns the photon position that produces the smallest coincidence time
TLorentzVector getKpi_mm(int i){return loc_kpi_mm->at(i);} // 4-vector g n ->K+ pi- X
TLorentzVector getK_mm(int i){return loc_k_mm->at(i);} // 4-vector g n ->K+ X
TLorentzVector getPpi_mm(int i){return loc_ppi_mm->at(i);} // 4-vector g n ->p pi- X when kaon is given proton mass
TLorentzVector getPipi_mm(int i){return loc_pipi_mm->at(i);} // 4-vector g n ->pi+ pi- X when kaon is given pion mass
TLorentzVector getDKppi_mm(int i){return loc_d_kppi_mm->at(i);} // 4-vector g d ->p K+ pi- X
TLorentzVector getDKp_mm(int i){return loc_d_kp_mm->at(i);} // 4-vector g d ->p K+ X
TVector3 getEVNT_vertex(int i){return loc_EVNT_vertex->at(i);} //EVNT vertex of i track
TVector3 getMVRT_vertex(){return *loc_MVRT_vertex;} //MVRT vertex
int getNextEntry();
int getEntry(){return eventno;}
int getEntries();

*/
