////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jun 13 14:33:54 2016 by ROOT version 6.06/04
// from TTree LHEF/Analysis tree
// found on file: WbT_M1.5_hadronic_pythia_lhe_events.root
//////////////////////////////////////////////////////////

#ifndef tprimeAnalisis_h
#define tprimeAnalisis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TClonesArray.h"
#include "TObject.h"

// private data handlers
#include "TLorentzVector.h"
#include <vector>

const string myrootname0("Wbt_tH_M800_hadronic");//"WbT_M1_pythia_lhe_events");//"WbT_M800G_hadronic_pythia_lhe_events");
const string myrootname1("Wbt_tH_M1_hadronic");//"WbT_M1_pythia_lhe_events");//"WbT_M800G_hadronic_pythia_lhe_events");
const string myrootname2("Wbt_tH_M1.2_hadronic");//"WbT_M1_pythia_lhe_events");//"WbT_M800G_hadronic_pythia_lhe_events");
const string myrootname3("Wbt_tH_M1.5_hadronic");//"WbT_M1_pythia_lhe_events");//"WbT_M800G_hadronic_pythia_lhe_events");
const string myrootname4("Wbt_tH_M1.8_hadronic");//"WbT_M1_pythia_lhe_events");//"WbT_M800G_hadronic_pythia_lhe_events");
const string myrootname5("Wbt_tH_M2_hadronic");//"WbT_M1_pythia_lhe_events");//"WbT_M800G_hadronic_pythia_lhe_events");
const string myrootname6("Wbt_tH_M2.2_hadronic");//"WbT_M1_pythia_lhe_events");//"WbT_M800G_hadronic_pythia_lhe_events");
const string myrootname7("Wbt_tH_M2.5_hadronic");//"WbT_M1_pythia_lhe_events");//"WbT_M800G_hadronic_pythia_lhe_events");
const string namebody("_pythia_lhe_events");

class tprimeAnalisis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   string	   rootfilename;
// Fixed size dimensions of array or collections stored in the TTree if any.
   static const Int_t kMaxEvent = 1;
   static const Int_t kMaxRwgt = 1;
   static const Int_t kMaxParticle = 29;

   // Declaration of leaf types
   Int_t           Event_;
   UInt_t          Event_fUniqueID[kMaxEvent];   //[Event_]
   UInt_t          Event_fBits[kMaxEvent];   //[Event_]
   Long64_t        Event_Number[kMaxEvent];   //[Event_]
   Int_t           Event_Nparticles[kMaxEvent];   //[Event_]
   Int_t           Event_ProcessID[kMaxEvent];   //[Event_]
   Double_t        Event_Weight[kMaxEvent];   //[Event_]
   Double_t        Event_ScalePDF[kMaxEvent];   //[Event_]
   Double_t        Event_CouplingQED[kMaxEvent];   //[Event_]
   Double_t        Event_CouplingQCD[kMaxEvent];   //[Event_]
   Int_t           Event_size;
   Int_t           Rwgt_;
   UInt_t          Rwgt_fUniqueID[kMaxRwgt];   //[Rwgt_]
   UInt_t          Rwgt_fBits[kMaxRwgt];   //[Rwgt_]
   Double_t        Rwgt_Weight[kMaxRwgt];   //[Rwgt_]
   Int_t           Rwgt_size;
   Int_t           Particle_;
   UInt_t          Particle_fUniqueID[kMaxParticle];   //[Particle_]
   UInt_t          Particle_fBits[kMaxParticle];   //[Particle_]
   Int_t           Particle_PID[kMaxParticle];   //[Particle_]
   Int_t           Particle_Status[kMaxParticle];   //[Particle_]
   Int_t           Particle_Mother1[kMaxParticle];   //[Particle_]
   Int_t           Particle_Mother2[kMaxParticle];   //[Particle_]
   Int_t           Particle_ColorLine1[kMaxParticle];   //[Particle_]
   Int_t           Particle_ColorLine2[kMaxParticle];   //[Particle_]
   Double_t        Particle_Px[kMaxParticle];   //[Particle_]
   Double_t        Particle_Py[kMaxParticle];   //[Particle_]
   Double_t        Particle_Pz[kMaxParticle];   //[Particle_]
   Double_t        Particle_E[kMaxParticle];   //[Particle_]
   Double_t        Particle_M[kMaxParticle];   //[Particle_]
   Double_t        Particle_PT[kMaxParticle];   //[Particle_]
   Double_t        Particle_Eta[kMaxParticle];   //[Particle_]
   Double_t        Particle_Phi[kMaxParticle];   //[Particle_]
   Double_t        Particle_Rapidity[kMaxParticle];   //[Particle_]
   Double_t        Particle_LifeTime[kMaxParticle];   //[Particle_]
   Double_t        Particle_Spin[kMaxParticle];   //[Particle_]
   Int_t           Particle_size;

   // List of branches
   TBranch        *b_Event_;   //!
   TBranch        *b_Event_fUniqueID;   //!
   TBranch        *b_Event_fBits;   //!
   TBranch        *b_Event_Number;   //!
   TBranch        *b_Event_Nparticles;   //!
   TBranch        *b_Event_ProcessID;   //!
   TBranch        *b_Event_Weight;   //!
   TBranch        *b_Event_ScalePDF;   //!
   TBranch        *b_Event_CouplingQED;   //!
   TBranch        *b_Event_CouplingQCD;   //!
   TBranch        *b_Event_size;   //!
   TBranch        *b_Rwgt_;   //!
   TBranch        *b_Rwgt_fUniqueID;   //!
   TBranch        *b_Rwgt_fBits;   //!
   TBranch        *b_Rwgt_Weight;   //!
   TBranch        *b_Rwgt_size;   //!
   TBranch        *b_Particle_;   //!
   TBranch        *b_Particle_fUniqueID;   //!
   TBranch        *b_Particle_fBits;   //!
   TBranch        *b_Particle_PID;   //!
   TBranch        *b_Particle_Status;   //!
   TBranch        *b_Particle_Mother1;   //!
   TBranch        *b_Particle_Mother2;   //!
   TBranch        *b_Particle_ColorLine1;   //!
   TBranch        *b_Particle_ColorLine2;   //!
   TBranch        *b_Particle_Px;   //!
   TBranch        *b_Particle_Py;   //!
   TBranch        *b_Particle_Pz;   //!
   TBranch        *b_Particle_E;   //!
   TBranch        *b_Particle_M;   //!
   TBranch        *b_Particle_PT;   //!
   TBranch        *b_Particle_Eta;   //!
   TBranch        *b_Particle_Phi;   //!
   TBranch        *b_Particle_Rapidity;   //!
   TBranch        *b_Particle_LifeTime;   //!
   TBranch        *b_Particle_Spin;   //!
   TBranch        *b_Particle_size;   //!

   bool 	  bkgrd;
   int		  cutlevel;
   bool 	  dodrhig;
   bool		  dodrtop;
   int 		  ht_cut;
   int 		  top_cut;
   int 		  hig_cut;

   string higdr;
   string topdr;
   string cuttype; // = "_ht" + itostr( ht_cut ) + "_hig" + itostr( hig_cut ) + "dr" + higdr + "_top" + itostr( top_cut ) + "dr" + topdr;

    int htcut;
    int higgscut;
    int topcut;
    int xjetcut;

    string cutname[5];// = { "_none", "_htcut", "_higgscut", "_topcut", "_xjetcut" };

    Long64_t nentries; // = fChain->GetEntriesFast();


/////////////////////////////////////////

    int dnc;
    int none;
    int down;
    int up;
    int strange;
    int charm;
    int bottom;
    int top;
    int bprime;
    int tprime;

    int electron;
    int enutrino;
    int muon;
    int mnutrino;
    int tau;
    int tnutrino;

    int gluon;
    int photon;
    int zee;
    int w;
    int higgs;
    int hplus;



//////////////////////////////////////

    int qwcnt;  //  quark from W bosun
    int Tcnt; //  Tprime
    int tTcnt; // top from a Tprime
    int tcnt;  //  top ( Parantage not checked )
    int bhTcnt; // bottom from higgs
    int btTcnt; //  bottom from top
    int bcnt; // bottom ( paratage not checked ) 
    int hTcnt;  //  higgs from Tprime
    int wtTcnt; // W from a top
    int qqcnt; // quark from a quark
    int bgcnt; // bottom from a gluon
    int evcnt; //  event count
    int ltjetcnt; //  ljet with W parent
    int otherjetcnt; // ljet ( paratage not checked )

    double Ht;
    double Ht_pe;
    int htjetsize;
    int part;

    vector<TLorentzVector> qwjet;  //  q from w
    vector<TLorentzVector> bhTjet; //  b from h
    vector<TLorentzVector> qqgjet;  //  q from q 
    vector<TLorentzVector> tjet;
    vector<TLorentzVector> bjet;
    vector<TLorentzVector> otherjet;
    vector<TLorentzVector> wjet;
    vector<TLorentzVector> transjet;
    vector<TLorentzVector> htjet;

    vector<int> bglist;
    vector<int> qjet;
    vector<int> qglist;
    vector<int> ljet;

    int ljetpartcnt[17]; // = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    int otherpartcnt[17]; // = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    int allcnt[25]; // = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    vector<TH1D*> thlist;
    vector<TH2D*> th2list;

        //TLorentzVector Declarations
        TLorentzVector tTvec;
        TLorentzVector wtTvec;
        TLorentzVector hTvec;
        TLorentzVector Tvec;
        TLorentzVector btTvec;
        TLorentzVector bgvec;

        TLorentzVector vReturn;
        int charg;

     Long64_t nbytes, nb;

       // delta_R_top
        double dr1;
        double dr2;
        double dr3;
        double dr_tmax;
        double dr_tmin;
        double dr_wb;
        double dr_h;
        double dr_w;

        Long64_t ientry;

        string fileName;
        char* temp;



   tprimeAnalisis( TTree *tree = 0 );
   tprimeAnalisis( string rootname  );
   virtual ~tprimeAnalisis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual bool	    get2Vector( int x, int particle, int mother, int granny, TLorentzVector &vectorReturn, int &charg );
   virtual bool	    get2VectorFromJet( int x, vector<int> &partlist, int mother, int granny, TLorentzVector &vectorReturn, int &charg, int &part );
   virtual bool     get2VectorFromLists( int numpart, vector<int> particle, vector<int> mother, int granny, TLorentzVector &vectorReturn, int &charg, int &part );
   virtual bool     passcut( int cutrun, TLorentzVector &hTvec, double dr_h, TLorentzVector &tTvec, double dr_wb, double Ht, int &htcut, int &higgscut, int &topcut  );
   virtual bool     cktopcut( TLorentzVector &tTvec, double dr_wb );
   virtual bool     ckhiggscut( TLorentzVector &hTvec, double dr_h );
   virtual double   pt30_eta5_HtBounds( TLorentzVector &vec );
   virtual bool     drcheckhig( double dr );
   virtual bool     drchecktop( double dr );

};

#endif

#ifdef tprimeAnalisis_cxx
tprimeAnalisis::tprimeAnalisis( TTree *tree ) : fChain(0) 
{
// WbT_M2.5_hadronic_pythia_events.lhe
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   string myrootfile = myrootname0 + ".root";
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject( myrootfile.c_str() );
      if (!f || !f->IsOpen()) {
         f = new TFile( myrootfile.c_str());
      }
      f->GetObject("LHEF",tree);

   }
   Init(tree);
   rootfilename = myrootname0;
}

tprimeAnalisis::tprimeAnalisis( string rootname ) : fChain(0)
{
// WbT_M2.5aaaa_hadronic_pythia_events.lhe
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   TTree *tree = 0;

   string myrootfile = rootname + "/" + rootname + namebody + ".root";
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject( myrootfile.c_str() );
      if (!f || !f->IsOpen()) {
         f = new TFile( myrootfile.c_str());
      }
      f->GetObject("LHEF",tree);
   }
   Init(tree);
   rootfilename = rootname;

}


tprimeAnalisis::~tprimeAnalisis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t tprimeAnalisis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t tprimeAnalisis::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void tprimeAnalisis::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   bkgrd = false;
   cutlevel = 3;
   dodrhig = true;
   dodrtop = true;
   ht_cut = 1100;
   top_cut = 400;
   hig_cut = 300;	

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Event", &Event_, &b_Event_);
   fChain->SetBranchAddress("Event.fUniqueID", Event_fUniqueID, &b_Event_fUniqueID);
   fChain->SetBranchAddress("Event.fBits", Event_fBits, &b_Event_fBits);
   fChain->SetBranchAddress("Event.Number", Event_Number, &b_Event_Number);
   fChain->SetBranchAddress("Event.Nparticles", Event_Nparticles, &b_Event_Nparticles);
   fChain->SetBranchAddress("Event.ProcessID", Event_ProcessID, &b_Event_ProcessID);
   fChain->SetBranchAddress("Event.Weight", Event_Weight, &b_Event_Weight);
   fChain->SetBranchAddress("Event.ScalePDF", Event_ScalePDF, &b_Event_ScalePDF);
   fChain->SetBranchAddress("Event.CouplingQED", Event_CouplingQED, &b_Event_CouplingQED);
   fChain->SetBranchAddress("Event.CouplingQCD", Event_CouplingQCD, &b_Event_CouplingQCD);
   fChain->SetBranchAddress("Event_size", &Event_size, &b_Event_size);
   fChain->SetBranchAddress("Rwgt", &Rwgt_, &b_Rwgt_);
   fChain->SetBranchAddress("Rwgt.fUniqueID", &Rwgt_fUniqueID, &b_Rwgt_fUniqueID);
   fChain->SetBranchAddress("Rwgt.fBits", &Rwgt_fBits, &b_Rwgt_fBits);
   fChain->SetBranchAddress("Rwgt.Weight", &Rwgt_Weight, &b_Rwgt_Weight);
   fChain->SetBranchAddress("Rwgt_size", &Rwgt_size, &b_Rwgt_size);
   fChain->SetBranchAddress("Particle", &Particle_, &b_Particle_);
   fChain->SetBranchAddress("Particle.fUniqueID", Particle_fUniqueID, &b_Particle_fUniqueID);
   fChain->SetBranchAddress("Particle.fBits", Particle_fBits, &b_Particle_fBits);
   fChain->SetBranchAddress("Particle.PID", Particle_PID, &b_Particle_PID);
   fChain->SetBranchAddress("Particle.Status", Particle_Status, &b_Particle_Status);
   fChain->SetBranchAddress("Particle.Mother1", Particle_Mother1, &b_Particle_Mother1);
   fChain->SetBranchAddress("Particle.Mother2", Particle_Mother2, &b_Particle_Mother2);
   fChain->SetBranchAddress("Particle.ColorLine1", Particle_ColorLine1, &b_Particle_ColorLine1);
   fChain->SetBranchAddress("Particle.ColorLine2", Particle_ColorLine2, &b_Particle_ColorLine2);
   fChain->SetBranchAddress("Particle.Px", Particle_Px, &b_Particle_Px);
   fChain->SetBranchAddress("Particle.Py", Particle_Py, &b_Particle_Py);
   fChain->SetBranchAddress("Particle.Pz", Particle_Pz, &b_Particle_Pz);
   fChain->SetBranchAddress("Particle.E", Particle_E, &b_Particle_E);
   fChain->SetBranchAddress("Particle.M", Particle_M, &b_Particle_M);
   fChain->SetBranchAddress("Particle.PT", Particle_PT, &b_Particle_PT);
   fChain->SetBranchAddress("Particle.Eta", Particle_Eta, &b_Particle_Eta);
   fChain->SetBranchAddress("Particle.Phi", Particle_Phi, &b_Particle_Phi);
   fChain->SetBranchAddress("Particle.Rapidity", Particle_Rapidity, &b_Particle_Rapidity);
   fChain->SetBranchAddress("Particle.LifeTime", Particle_LifeTime, &b_Particle_LifeTime);
   fChain->SetBranchAddress("Particle.Spin", Particle_Spin, &b_Particle_Spin);
   fChain->SetBranchAddress("Particle_size", &Particle_size, &b_Particle_size);
   Notify();
}

Bool_t tprimeAnalisis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void tprimeAnalisis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t tprimeAnalisis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}






#endif // #ifdef tprimeAnalisis_cxx
