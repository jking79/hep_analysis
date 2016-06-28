#define tprimeAnalisis_cxx
#include "tprimeAnalyis.h"
#include <TCanvas.h>
#include "jwk_ku_tdr_style.h"
#include "CombineHistograms.h"


bool tprimeAnalisis::get2Vector( int numpart, int particle, int mother, int granny, TLorentzVector &vectorReturn, int &charg )
{

    	bool results = false;
//    	cout << "Looking for 2 " << particle << " "<< mother << " " << granny << endl;

	int part = Particle_PID[numpart];
	int mom = Particle_PID[Particle_Mother1[numpart]];
	int gran = Particle_PID[Particle_Mother1[Particle_Mother1[numpart]]];
	int level = Particle_Status[numpart];	

//    	cout << "Found       " << level << " " << part << " " << mom << " " << gran << endl;
    	if( level == 2  ){
        	if( abs(part) == particle ){  
        		if( ( mother == 0 ) || ( abs( mom ) == mother ) ){ 
        			if( ( granny == 0 ) || ( abs( gran ) == granny )){
        				vectorReturn.SetPtEtaPhiM(Particle_PT[numpart], Particle_Eta[numpart], Particle_Phi[numpart], Particle_M[numpart]);
        				results = true;
        				if( part < 0 ) charg = -1; else charg = 1;
    				}	
			}	
		}
    	}		
//  	cout << "Match? :  " << results << endl;
    	return results;
}

bool tprimeAnalisis::get2VectorFromJet( int x, vector<int> &partlist, int mother, int granny, TLorentzVector &vectorReturn, int &charg )
{

	int numparts = partlist.size();
	bool results = false;
	int jetloop = 0;
	while( (!results) && (jetloop < numparts) ){
//		cout << "In jet Looking for : " << partlist[jetloop]<<endl;
		results = get2Vector( x, partlist[jetloop], mother, granny, vectorReturn, charg );
//		if( results ){results = false; cout << "jet loop result for " << x << endl;} 
		jetloop++;
	}
	
//	cout << "jet Match? : " << results << endl;
	return results;

}

bool tprimeAnalisis::get2VectorFromLists( int numpart, vector<int> particle, vector<int> mother, int granny, TLorentzVector &vectorReturn, int &charg )
{

	int numparts = mother.size();
        bool results = false;
        int jetloop = 0;
        while( (!results) && (jetloop < numparts) ){
//              cout << "In jet Looking for : " << partlist[jetloop]<<endl;
                results = get2VectorFromJet( numpart, particle, mother[jetloop], granny, vectorReturn, charg );
//              if( results ){results = false; cout << "jet loop result for " << x << endl;} 
                jetloop++;
        }

//      cout << "jet Match? : " << results << endl;
        return results;

}

TCanvas* drawPrint( TH1D* hist, string type, string title, string xtitle, string ytitle )
{
	TCanvas *c = getTDRCanvas( title );
	tdrHistDraw( hist, c, xtitle, ytitle );
	string savetitle( title + type );
	c->SaveAs( savetitle.c_str() );
	return c;
}

TCanvas* drawPrint( TH1D* hist, string title, string xtitle, string ytitle )
{
        return drawPrint( hist, ".png", title, xtitle, ytitle );
}

TCanvas* drawPrint2D( TH2D* hist, string type, string title, string xtitle, string ytitle )
{
        TCanvas *c = getTDRCanvas( title );
        tdrHist2DDraw( hist, c, xtitle, ytitle );
        string savetitle( title + type );
        c->SaveAs( savetitle.c_str() );
        return c;
}

TCanvas* drawPrint2D( TH2D* hist, string title, string xtitle, string ytitle )
{
        return drawPrint2D( hist, ".png", title, xtitle, ytitle );
}

TCanvas* drawPrintComb( TH1D* hist1, TH1D* hist2, TH2D* hist3, string type, string title )
{
        TCanvas *c = getTDRCanvas( title );
	CombineHistograms( hist1, hist2, hist3, c );
        string savetitle( title + type );
        c->SaveAs( savetitle.c_str() );
        return c;
}

TCanvas* drawPrintComb( TH1D* hist1, TH1D* hist2, TH2D* hist3, string title )
{
        return drawPrintComb( hist1, hist2, hist3, ".png", title );
}

void tprimeAnalisis::Loop()
{
   if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();

    int none = 0;

/////////////////////////////////////////

    int down = 1;
    int up = 2;
    int strange = 3;
    int charm = 4;
    int bottom = 5;
    int top = 6;
    int bprime = 7;
    int tprime = 8000001;

    int electron = 11;
    int enutrino = 12; //5555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555
    int muon = 13;
    int mnutrino = 14;
    int tau = 15;
    int tnutrino = 16;

    int gluon = 21;
    int photon = 22;
    int zee = 23;
    int w = 24;
    int higgs = 25;
    int hplus = 37;

//////////////////////////////////////

    int qwcnt = 0;
    int Tcnt = 0;
    int tTcnt = 0;
    int bhTcnt = 0;
    int btTcnt = 0;
    int hTcnt = 0;
    int wtTcnt = 0;
    int qqcnt = 0;
    int bgcnt = 0;
    int evcnt = 0;

//    int qqstrgcnt = 0;
//    int qqcharmcnt = 0;
//    int qwstrgcnt = 0;
//    int qwcharmcnt = 0;

    vector<TLorentzVector> qwjet;
    vector<TLorentzVector> bhTjet;
    vector<TLorentzVector> qwtjet;
    
    vector<int> bglist;
    bglist.push_back(bottom);

    vector<int> qjet;
    qjet.push_back(down);
    qjet.push_back(up);
    qjet.push_back(strange);
    qjet.push_back(charm);
//    qjet.push_back(top);
//    qjet.push_back(bottom);

    vector<int> qglist;
    qglist.push_back(down);
    qglist.push_back(up);
    qglist.push_back(strange);
    qglist.push_back(charm);
    qglist.push_back(gluon);
//    qglist.push_back(top);
//    qglist.push_back(bottom);

//bg        
    TH1D *pt_bg_hist = new TH1D("Pt of b from Gluon", "Pt of b from Gluon", 50, 0, 1000);
    TH1D *im_bg_hist = new TH1D( "Mass of b from Gluon","Mass of b from Gluon", 100, 0, 10 );
    TH1D *eta_bg_hist = new TH1D("Eta of b from Gluon", "Eta of b from Gluon", 50, -5, 5);
//qq        
    TH1D *pt_qq_hist = new TH1D("Pt of q from qW", "Pt of q from qW", 50, 0, 1000);
    TH1D *im_qq_hist = new TH1D( "Mass of q from qW","Mass of q from qW", 100, 0, 2 );
    TH1D *eta_qq_hist = new TH1D("Eta of q from qW", "Eta of q from qW", 50, -5, 5);
//tT        
    TH1D *pt_top_hist = new TH1D("Pt Top", "Pt Top",100, 0, 2000);
    TH1D *im_top_hist = new TH1D( "Mass Top","Mass of Top", 50, 0, 1000 );
    TH1D *eta_top_hist = new TH1D("Eta Top", "Eta Top", 50, -5, 5);
//wt
    TH1D *pt_w_hist = new TH1D("Pt W", "Pt W", 50, 0, 1000);
    TH1D *eta_w_hist = new TH1D("Eta W", "Eta W", 50, -5, 5);    
    TH1D *im_w_hist = new TH1D( "Mass W","Mass of W", 100, 0, 200 );
//hT
    TH1D *pt_higgs_hist = new TH1D("Pt Higgs","Pt Higgs", 100, 0, 2000 );
    TH1D *eta_higgs_hist = new TH1D("Eta of Higgs", "Eta of HIggs", 50, -5,5 );
    TH1D *im_higgs_hist = new TH1D("Mass Higgs","Mass of the Higgs", 100, 100, 150 );
//T
    TH1D *im_tprime_hist = new TH1D( "Mass TPrime","Mass tPrime", 100, 2000, 4000 );
    TH1D *pt_tprime_hist = new TH1D( "Pt TPrime","Pt of tPrime", 50, 0, 1000 );    
    TH1D *eta_tprime_hist = new TH1D( "Eta TPrime","Eta of tPrime", 50, -5, 5 );
//qwt   <<<<<<<<<<<<
    TH1D *pt_qwt_hist = new TH1D("Pt qw Jet", "Pt qw Jet", 80, 0, 800);
    TH1D *im_qwt_hist = new TH1D( "Mass qw Jet","Mass of qw Jet", 100, 0, 0.2 );
    TH1D *eta_qwt_hist = new TH1D("Eta qw Jet", "Eta qw Jet", 50, -5, 5);
//btT
    TH1D *pt_btT_hist = new TH1D("Pt btT", "Pt btT", 50, 0, 1000);
    TH1D *im_btT_hist = new TH1D( "Mass btT","Mass of btT", 100, 0, 10 );
    TH1D *eta_btT_hist = new TH1D("Eta btT", "Eta btT", 50, -5, 5);
//bhT
    TH1D *pt_bhT_hist = new TH1D("Pt bh Jet", "Pt bh Jet", 50, 0, 1000);
    TH1D *im_bhT_hist = new TH1D( "Mass bh Jet","Mass of bh Jet", 100, 0, 10 );
    TH1D *eta_bhT_hist = new TH1D("Eta bh Jet", "Eta bh Jet", 50, -5, 5);
//dr
    TH1D *drmin_top_hist = new TH1D( "Min DeltaR Top", "Min DeltaR of Top", 40, 0, 4 );
    TH1D *drmax_top_hist = new TH1D( "Max DeltaR Top", "Max DeltaR of Top", 40, 0, 4 );
    TH1D *dr_higgs_hist = new TH1D( "DelaR Higgs", "DeltaR Higgs", 40, 0, 4 );
    TH1D *dr_w_hist = new TH1D( "DeltaR W+", "DeltaR of W+", 40, 0, 4 );
    TH1D *dr_wb_hist = new TH1D( "DeltaR of wb Jet", "DeltaR of wb Jet", 40, 0, 4 );
//dr vs Pt
    TH2D *drvpt_top_hist = new TH2D( "DeltaR vs Pt of Top", "DeltaR vs Pt of Top", 50, 0, 1000, 40, 0, 4 );
    TH2D *drvpt_higgs_hist = new TH2D( "DeltaR vs Pt of Higgs", "DeltaR vs Pt of Higgs", 50, 0, 1000, 40, 0, 4 );
    TH2D *drvpt_w_hist = new TH2D( "DeltaR vs Pt of W", "DeltaR vs Pt of W", 50, 0, 1000, 40, 0, 4 );
//Ht
    TH1D *Ht_hist = new TH1D( "Ht", "Ht", 50, 0, 2000 );


    Long64_t nbytes = 0, nb = 0;
    //loop through events
    cout << "Number of events : " << nentries << endl;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;

        //TLorentzVector Declarations
        TLorentzVector tTvec;
        TLorentzVector wtTvec;
        TLorentzVector hTvec;
        TLorentzVector Tvec;
        TLorentzVector btTvec;
	TLorentzVector bgvec;
	TLorentzVector qqvec;

        TLorentzVector vReturn;
        int charg = 0;

        //Loop Through Particles
//	cout << "Number of particles in an event : " << *Event_Nparticles << endl;
        for(int partnum=0; partnum < *Event_Nparticles; partnum++){ // partnum = particle number

	    vReturn.SetPtEtaPhiM(0,0,0,0);
	    charg = 0;
            if( get2VectorFromLists( partnum, bglist, qglist, none, vReturn, charg )){ bgcnt++; bgvec = vReturn;}
            if( get2VectorFromLists( partnum, qglist, qglist, none, vReturn, charg )){ qqcnt++; qqvec = vReturn;}
            if( get2Vector( partnum, top, tprime, none, vReturn, charg )){ tTcnt++; tTvec = vReturn;}
            if( get2Vector( partnum, w, top, tprime, vReturn, charg )) { wtTcnt++; wtTvec = vReturn; }            
	    if( get2Vector( partnum, higgs, tprime, none, vReturn, charg )) { hTcnt++; hTvec = vReturn;}
	    if( get2Vector( partnum, tprime, none, none, vReturn, charg )) { Tcnt++; Tvec = vReturn;}
            if( get2Vector( partnum, bottom, top, tprime, vReturn, charg )) {btTcnt++; btTvec = vReturn;}
            if( get2Vector( partnum, bottom, higgs, tprime, vReturn, charg )){ bhTcnt++; bhTjet.push_back( vReturn ); } 
	    if( get2VectorFromJet( partnum, qjet, w, none, vReturn, charg )){ qwcnt++; qwjet.push_back( vReturn );}

        }// find particles in above loop

// calc varibles
	evcnt++;
	// delta_R_top
	double dr1 = (qwjet[0]).DeltaR( qwjet[1] );
	double dr2 = (qwjet[1]).DeltaR( btTvec );
	double dr3 = (qwjet[0]).DeltaR( btTvec );
	double dr_tmax = 0;
        double dr_tmin = 0;
        double dr_wb = wtTvec.DeltaR( btTvec );
	if( dr1 > dr2 ){ if( dr1 > dr3 ){ dr_tmax = dr1;} else { if( dr2 > dr3 ){ dr_tmax = dr2;} else{ dr_tmax = dr3; }}}	
	if( dr1 < dr2 ){ if( dr1 < dr3 ){ dr_tmin = dr1;} else { if( dr2 < dr3 ){ dr_tmin = dr2;} else{ dr_tmin = dr3; }}} 
	double dr_h = (bhTjet[0]).DeltaR( bhTjet[1] );
	double dr_w = dr1;

	double Ht = bgvec.Pt() + qqvec.Pt() + tTvec.Pt() + wtTvec() + hTvec() + Tvec.Pt() + btTvec().Pt() + bhTjet[0].Pt() + bhTjet[1].Pt()+ qwjet[0].Pt() + qwjet[1].Pt();
	
	qwjet.clear();
	bhTjet.clear();
//	cout << "momentum " << tVec.Pt() << endl;

        pt_bg_hist->Fill( bgvec.Pt() );
        eta_bg_hist->Fill( bgvec.Eta() );
        im_bg_hist->Fill( bgvec.M() );

        pt_qq_hist->Fill( qqvec.Pt() );
        eta_qq_hist->Fill( qqvec.Eta() );
        im_qq_hist->Fill( qqvec.M() );

        pt_top_hist->Fill( tTvec.Pt() );
        eta_top_hist->Fill( tTvec.Eta() );
        im_top_hist->Fill( tTvec.M() );

        pt_w_hist->Fill( wtTvec.Pt() );
        eta_w_hist->Fill( wtTvec.Eta() );
        im_w_hist->Fill( wtTvec.M() );

	pt_higgs_hist->Fill( hTvec.Pt() );
        eta_higgs_hist->Fill( hTvec.Eta() );
  	im_higgs_hist->Fill( hTvec.M() );

	pt_tprime_hist->Fill( Tvec.Pt() );
        eta_tprime_hist->Fill( Tvec.Eta() );
        im_tprime_hist->Fill( Tvec.M() );

	for( int c = 0; c < 2; c++ ){
        	pt_qwt_hist->Fill( qwjet[c].Pt() );
        	eta_qwt_hist->Fill( qwjet[c].Eta() );
        	im_qwt_hist->Fill( qwjet[c].M() );

        	pt_bhT_hist->Fill( bhTjet[c].Pt() );
        	eta_bhT_hist->Fill( bhTjet[c].Eta() );
        	im_bhT_hist->Fill( bhTjet[c].M() );
	}

        pt_btT_hist->Fill( btTvec.Pt() );
        eta_btT_hist->Fill( btTvec.Eta() );
        im_btT_hist->Fill( btTvec.M() );

	dr_top_hist->Fill( dr_tmin );
	dr_top_hist->Fill( dr_tmax );
	dr_wb_hist->Fill( dr_wb );
        dr_w_hist->Fill( dr_w );
	dr_higgs_hist->Fill( dr_h );

        drvpt_top_hist->Fill( tTvec.Pt(), dr_t );
	drvpt_higgs_hist->Fill( hTvec.Pt(), dr_h );
        drvpt_w_hist->Fill( wtTvec.Pt(), dr_w );
	
	Ht_hist->Fill( Ht );

    }// calc values and fill histograms

    setTDRStyle();

    TCanvas *c1 = drawPrint( pt_top_hist,"Pt of the Top", "Pt (GeV)", "Events/20GeV");   
    TCanvas *c2 = drawPrint( eta_top_hist,"Eta of the Top", "Eta", "Events/0.2");
    TCanvas *c3 = drawPrint( im_top_hist,"Mass of Top", "Mass (GeV)", "Events/20GeV");

    TCanvas *c4 = drawPrint( im_w_hist,"Mass of the W", "Mass (GeV)", "Events/2GeV");    
    TCanvas *c5 = drawPrint( pt_w_hist,"Pt of the W", "Pt (GeV)", "Events/20GeV");   
    TCanvas *c6 = drawPrint( eta_w_hist,"Eta of the W", "Eta", "Events/0.2");

    TCanvas *c7 = drawPrint( pt_higgs_hist, "Pt of Higgs", "Pt (GeV)", "Events/20GeV");
    TCanvas *c8 = drawPrint( eta_higgs_hist, "Eta of the Higgs", "Eta", "Events/0.2");    
    TCanvas *c9 = drawPrint( im_higgs_hist, "Mass of the Higgs", "Mass (GeV)", "Events/20GeV");

    TCanvas *c10 = drawPrint( pt_tprime_hist, "Pt of TPrime", "Pt (GeV)", "Events/20GeV");
    TCanvas *c11 = drawPrint( eta_tprime_hist, "Eta of TPrime", "Eta", "Events/0.2");
    TCanvas *c12 = drawPrint( im_tprime_hist, "Mass of TPrime", "Mass (GeV)", "Events/20GeV");

    TCanvas *c13 = drawPrint( dr_top_hist, "Delat R of Top", "DeltaR", "Events/0.2");
    TCanvas *c14 = drawPrint( dr_w_hist, "Delta R of W", "DeltaR", "Events/0.2" );
    TCanvas *c15 = drawPrint( dr_higgs_hist, "Delta R of Higgs", "DeltaR", "Events/0.2" );
    
    TCanvas *c16 = drawPrint2D(drvpt_top_hist, "DeltaR vs Pt of Top", "Pt (GeV)", "DeltaR");
    TCanvas *c17 = drawPrint2D(drvpt_higgs_hist, "DeltaR vs Pt of Higgs", "Pt (GeV)", "DeltaR");
    TCanvas *c18 = drawPrint2D(drvpt_w_hist, "DeltaR vs Pt of W", "Pt (GeV)", "DeltaR");

    TCanvas *c19 = drawPrint( im_qwt_hist, "Mass of qwt", "Mass (GeV)", "Events/0.002GeV");
    TCanvas *c20 = drawPrint( eta_qwt_hist, "Eta of qwt", "Eta", "Events/0.2");
    TCanvas *c21 = drawPrint( pt_qwt_hist, "Pt  of qwt", "Pt (GeV)", "Events/10GeV" );

    TCanvas *c22 = drawPrint( im_btT_hist, "Mass of btT", "Mass (GeV)", "Events/20GeV");
    TCanvas *c23 = drawPrint( eta_btT_hist, "Eta of btT", "Eta", "Events/0.2");
    TCanvas *c24 = drawPrint( pt_btT_hist, "Pt of btT", "Pt (GeV)", "Events/20GeV" );

    TCanvas *c25 = drawPrint( im_bhT_hist, "Mass of bht", "Mass (GeV)", "Events/20GeV");
    TCanvas *c26 = drawPrint( eta_bhT_hist, "Eta of bht", "Eta", "Events/0.2");
    TCanvas *c27 = drawPrint( pt_bhT_hist, "Pt of bht", "Pt (GeV)", "Events/0.10GeV" );

    TCanvas *c28 = drawPrint( pt_bg_hist, "Pt of b from Gluon","Pt (GeV)", "Events/20GeV" );
    TCanvas *c29 = drawPrint( eta_bg_hist, "Eta of b from Gluon","Eta", "Events/0.2" );
    TCanvas *c30 = drawPrint( im_bg_hist, "Mass of b from Gluon","Mass (GeV)", "Events/0.1GeV" );
    
    TCanvas *c31 = drawPrint( pt_qq_hist, "Pt of q from qW","Pt (GeV)", "Events/20GeV" );
    TCanvas *c32 = drawPrint( eta_qq_hist, "Eta of q from qW","Eta", "Events/0.2" );
    TCanvas *c33 = drawPrint( im_qq_hist, "Mass of q from qW","Mass (GeV)", "Events/0.02GeV" );

    TCanvas *c34 = drawPrintComb( pt_top_hist, dr_top_hist, drvpt_top_hist, "DeltaR vs Pt of Top Combo" ); 
    TCanvas *c35 = drawPrintComb( pt_higgs_hist, dr_higgs_hist, drvpt_higgs_hist, "DeltaR vs Pt of Higgs Combo" );
    TCanvas *c36 = drawPrintComb( pt_w_hist, dr_w_hist, drvpt_w_hist, "DeltaR vs Pt of W Combo" );
   
    TCanvas *c37 = drawPrint( Ht_hist, "Ht for Event Run","Ht (GeV)","Events/50GeV"); 

    cout << "Number of qw  found: " << qwcnt << endl;
    cout << "Number of T   found: " << Tcnt << endl;/////
    cout << "Number of tT  found: " << tTcnt << endl;/////
    cout << "Number of bhT found: " << bhTcnt << endl;/////
    cout << "Number of btT found: " << btTcnt << endl;/////
    cout << "Number of hT  found: " << hTcnt << endl;//////
    cout << "Number of wtT found: " << wtTcnt << endl;/////
    cout << "Number of qq  found: " << qqcnt << endl;/////
    cout << "Number of bg  found: " << bgcnt << endl;
    cout << "Number of events found: " << evcnt << endl;
    return;

}//print histogram aboves

//   In a ROOT session, you can do:
//      root> .L tprimeAnalisis.C
//      root> tprimeAnalisis t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch