#define tprimeAnalisis_cxx
#include "tpana_onerun.h"
#include <TCanvas.h>
#include "jwk_ku_tdr_style.h"
#include "CombineHistograms.h"

static const bool dodr = false;
static const int ht_cut = 1000;
static const int top_cut = 400;
static const int hig_cut = 300;

bool tprimeAnalisis::get2Vector( int numpart, int particle, int mother, int granny, TLorentzVector &vectorReturn, int &charg )
{

    	bool results = false;
//    	cout << "Looking for 2 " << particle << " "<< mother << " " << granny << endl;

	int part = Particle_PID[numpart];
	int mom = Particle_PID[Particle_Mother1[numpart]];
	int gran = Particle_PID[Particle_Mother1[Particle_Mother1[numpart]]];
	int level = Particle_Status[numpart];	

  //  	cout << "Found       " << level << " " << part << " " << mom << " " << gran << endl;
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
  //	cout << "Match? :  " << results << endl;
    	return results;
}

bool tprimeAnalisis::get2VectorFromJet( int x, vector<int> &partlist, int mother, int granny, TLorentzVector &vectorReturn, int &charg, int &part )
{

	int numparts = partlist.size();
	bool results = false;
	int jetloop = 0;
	while( (!results) && (jetloop < numparts) ){
  //		cout << "In jet Looking for : " << partlist[jetloop]<<endl;
		results = get2Vector( x, partlist[jetloop], mother, granny, vectorReturn, charg );
		if( results ){part = partlist[jetloop];} 
		jetloop++;
	}
	
//	cout << "jet Match? : " << results << endl;
	return results;

}

bool tprimeAnalisis::get2VectorFromLists( int numpart, vector<int> particle, vector<int> mother, int granny, TLorentzVector &vectorReturn, int &charg, int &part )
{

	int numparts = mother.size();
        bool results = false;
        int jetloop = 0;
        while( (!results) && (jetloop < numparts) ){
      //          cout << "In jet Looking for : " << mother[jetloop]<<endl;
                results = get2VectorFromJet( numpart, particle, mother[jetloop], granny, vectorReturn, charg, part );
                jetloop++;
        }

    //    cout << "jet Match? : " << results << endl;
        return results;

}

TCanvas* drawPrintComb( TH1D* hist1, TH1D* hist2, TH2D* hist3, string type, string title, string rfname )
{
        TCanvas *c = getTDRCanvas( title + rfname );
	CombineHistograms( hist1, hist2, hist3, c );
        string savetitle( rfname + "_" + rpsb(title) + type );
        c->SaveAs( savetitle.c_str() );
        return c;
}

TCanvas* drawPrintComb( TH1D* hist1, TH1D* hist2, TH2D* hist3, string title, string rfname )
{
        return drawPrintComb( hist1, hist2, hist3, ".png", title, rfname );
}

bool drcheck( double dr )
{
	bool result = false;
	if( dodr ){
		if( dr < 0.8 ) result = true;
	} else { result = true; }
	return result;
}

double pt30_eta5_HtBounds( TLorentzVector &vec ) 
{
	double result = 0;
	if( ( vec.Pt() > 30 ) && ( abs( vec.Eta()) < 5 ) ) result = vec.Pt();
	return result; 
}

bool ckhiggscut( TLorentzVector &hTvec, double dr_h )
{
	bool result = false;
	if( (hTvec.Pt() > hig_cut ) && ( abs( hTvec.Eta() ) < 2.4 ) && drcheck(dr_h)  ) result = true;
	return result;
}

bool cktopcut( TLorentzVector &tTvec, double dr_wb ) //&& ( dr_wb < 0.8 )
{
	bool result = false;
	if( (tTvec.Pt() > top_cut ) && ( abs( tTvec.Eta() ) < 2.4 ) && drcheck( dr_wb )  ) result = true;
	return result;
}

bool passcut( int cutrun, TLorentzVector &hTvec, double dr_h, TLorentzVector &tTvec, double dr_wb, double Ht, int &htcut, int &higgscut, int &topcut ){

/*      if( Ht > 1100 ){ htcut++;
        if( (hTvec.Pt() > 300 ) && ( abs( hTvec.Eta() ) < 2.4 ) ){ higgscut++;
        if( (tTvec.Pt() > 400 ) && ( abs( tTvec.Eta() ) < 2.4 ) ){ topcut++;
*/
	bool result = false;

	if( cutrun == 0 ){ result = true; }
	else if( ( cutrun == 1 ) && ( Ht > ht_cut ) ) { htcut++; result = true; }
	else if( ( cutrun == 2 ) && ( Ht > ht_cut ) && ckhiggscut( hTvec, dr_h ) ){ higgscut++; result = true; }
        else if( ( cutrun == 3 ) && ( Ht > ht_cut ) && ckhiggscut( hTvec, dr_h ) && cktopcut( tTvec, dr_wb) ){ topcut++; result = true; }
	
	return result;
}

void tprimeAnalisis::Loop()
{
   if (fChain == 0) return;

  string cutname[4] = { "_none", "_htcut", "_higgscut", "_topcut" };
  for( int cutrun = 0; cutrun < (4); cutrun++ ){

    Long64_t nentries = fChain->GetEntriesFast();


/////////////////////////////////////////

    int none = 0;
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
    
    int qqw_qcnt[8] = { 0,0,0,0,0,0,0,0 };
    int wqq_qcnt[8] = { 0,0,0,0,0,0,0,0 };

    double Ht = 0;
    double Heta = 0; 
    int part = 0;

    int htcut = 0;
    int higgscut = 0;
    int topcut = 0;

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

//  for( int cutrun = 0; cutrun < (number_of_cuts + 1); cutrun++ ){ 
    
   cout << "Starting Cut :" << cutname[cutrun] << endl;
//bg        
    TH1D *pt_bg_hist = new TH1D("Pt of b from Gluon", "Pt of b from Gluon", 50, 0, 1000);
    TH1D *im_bg_hist = new TH1D( "Mass of b from Gluon","Mass of b from Gluon", 100, 0, 10 );
    TH1D *eta_bg_hist = new TH1D("Eta of b from Gluon", "Eta of b from Gluon", 50, -5, 5);
//qq        
    TH1D *pt_qq_hist = new TH1D("Pt of q from qW", "Pt of q from qW", 25, 0, 500);
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
    TH1D *im_tprime_hist = new TH1D( "Mass TPrime","Mass tPrime", 200, 0, 4000 );
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
    TH1D *dr_tmin_hist = new TH1D( "Min DeltaR Top", "Min DeltaR of Top", 40, 0, 4 );
    TH1D *dr_tmax_hist = new TH1D( "Max DeltaR Top", "Max DeltaR of Top", 40, 0, 4 );
    TH1D *dr_higgs_hist = new TH1D( "DelaR Higgs", "DeltaR Higgs", 40, 0, 4 );
    TH1D *dr_w_hist = new TH1D( "DeltaR W+", "DeltaR of W+", 40, 0, 4 );
    TH1D *dr_wb_hist = new TH1D( "DeltaR of wb Jet", "DeltaR of wb Jet", 40, 0, 4 );
//dr vs Pt
    TH2D *drvpt_top_hist = new TH2D( "DeltaR vs Pt of Top", "DeltaR vs Pt of Top", 100, 0, 2000, 40, 0, 4 );
    TH2D *drvpt_higgs_hist = new TH2D( "DeltaR vs Pt of Higgs", "DeltaR vs Pt of Higgs", 100, 0, 2000, 40, 0, 4 );
    TH2D *drvpt_w_hist = new TH2D( "DeltaR vs Pt of W", "DeltaR vs Pt of W", 50, 0, 1000, 40, 0, 4 );
//Ht
    TH1D *Ht_hist = new TH1D( "Ht", "Ht", 200, 0, 4000 );
    TH1D *Heta_hist = new TH1D( "SumEta", "SumEta", 50, -5, 5 );


    Long64_t nbytes = 0, nb = 0;
    //loop through events
//    cout << "Number of events : " << nentries << endl;
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
//	cout << "Looping Particles: Number of particles in the event : " << *Event_Nparticles << endl;
        for(int partnum=0; partnum < *Event_Nparticles; partnum++){ // partnum = particle number
  //          cout << " in loop with particle number :" << partnum << endl;
	    vReturn.SetPtEtaPhiM(0,0,0,0);
	    charg = 0;
            if( get2VectorFromLists( partnum, bglist, qglist, none, vReturn, charg, part )){ bgcnt++; bgvec = vReturn;}
            if( get2VectorFromLists( partnum, qglist, qglist, none, vReturn, charg, part )){ qqcnt++; qqvec = vReturn;}
            if( get2Vector( partnum, top, tprime, none, vReturn, charg )){ tTcnt++; tTvec = vReturn;}
            if( get2Vector( partnum, w, top, tprime, vReturn, charg )) { wtTcnt++; wtTvec = vReturn; }            
	    if( get2Vector( partnum, higgs, tprime, none, vReturn, charg )) { hTcnt++; hTvec = vReturn;}
	    if( get2Vector( partnum, tprime, none, none, vReturn, charg )) { Tcnt++; Tvec = vReturn;}
            if( get2Vector( partnum, bottom, top, tprime, vReturn, charg )) {btTcnt++; btTvec = vReturn;}
            if( get2Vector( partnum, bottom, higgs, tprime, vReturn, charg )){ bhTcnt++; bhTjet.push_back( vReturn ); } 
	    if( get2VectorFromJet( partnum, qjet, w, none, vReturn, charg, part )){ qwcnt++; qwjet.push_back( vReturn );}
      

        }// find particles in above loop

// calc varibles
	cout << "Calcating values for event" << endl;
	evcnt++;
	// delta_R_top
	double dr1 = 0;
	if( qwjet.size() == 2 ){ dr1 = (qwjet[0]).DeltaR( qwjet[1] ); } else { cout << "qwjet bad" << endl;}
        cout << "stop2" << endl;
	double dr2 = (qwjet[1]).DeltaR( btTvec );
	double dr3 = (qwjet[0]).DeltaR( btTvec );
	double dr_tmax = 0;
        double dr_tmin = 0;
        cout << "stop1"<<endl;
        double dr_wb = wtTvec.DeltaR( btTvec );
	if( dr1 > dr2 ){ if( dr1 > dr3 ){ dr_tmax = dr1;} else { dr_tmax = dr3;}} else { if( dr2 > dr3 ){ dr_tmax = dr2;} else{ dr_tmax = dr3; }}	
	if( dr1 < dr2 ){ if( dr1 < dr3 ){ dr_tmin = dr1;} else { dr_tmin = dr3;}} else { if( dr2 < dr3 ){ dr_tmin = dr2;} else{ dr_tmin = dr3; }} 
	double dr_h = (bhTjet[0]).DeltaR( bhTjet[1] );
	double dr_w = dr1;
//pt>30 and eta < 5 for inclusion to Ht
	Ht = pt30_eta5_HtBounds( bgvec ) + pt30_eta5_HtBounds( qqvec ) + pt30_eta5_HtBounds( btTvec );
	Ht = Ht  + pt30_eta5_HtBounds(bhTjet[0]) + pt30_eta5_HtBounds(bhTjet[1])+ pt30_eta5_HtBounds(qwjet[0]) + pt30_eta5_HtBounds(qwjet[1]);
	Heta = bgvec.Eta() + qqvec.Eta() + btTvec.Eta() + (bhTjet[0]).Eta() + (bhTjet[1]).Eta()+ (qwjet[0]).Eta() + (qwjet[1]).Eta();

	qwjet.clear();
	bhTjet.clear();
//	Fill cutts<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	
	cout << "Checking Cut for cut level :" << cutrun << endl;
	if( passcut( cutrun, hTvec, dr_h,  tTvec, dr_wb, Ht, htcut, higgscut, topcut ) ) {
	cout << "Filling histograms" << endl;
//      fills
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

	dr_tmin_hist->Fill( dr_tmin );
	dr_tmax_hist->Fill( dr_tmax );
	dr_wb_hist->Fill( dr_wb );
        dr_w_hist->Fill( dr_w );
	dr_higgs_hist->Fill( dr_h );

        drvpt_top_hist->Fill( tTvec.Pt(), dr_tmax);
	drvpt_higgs_hist->Fill( hTvec.Pt(), dr_h );
        drvpt_w_hist->Fill( wtTvec.Pt(), dr_w );
	
	Ht_hist->Fill( Ht );
	Heta_hist->Fill( Heta );
	

 	}//<<<<<<  cuts
       
  //      cout << "Event Finished" << endl;   
    }// calc values and fill histograms
    cout << "Print and Draw Histograms" << endl;
    vector<TH1D*> histlist;
    histlist.push_back( dr_tmax_hist );
    histlist.push_back( dr_tmin_hist );
    histlist.push_back( dr_wb_hist );
     
    vector<string> bob;
    bob.push_back( "Legend" );
    bob.push_back( "DeltaR of Tmax" );
    bob.push_back( "DeltaR of Tmin" );
    bob.push_back( "DeltaR of Wb" );

    setTDRStyle();
    int stop = rootfilename.length() - namebody.length();
    string savename = rootfilename.substr(0,stop ) + "/" +  rootfilename + cutname[cutrun];

    TCanvas *c1 = drawPrint( pt_top_hist,"Pt of the Top", "Pt (GeV)", "Events/20GeV", savename );   
    TCanvas *c2 = drawPrint( eta_top_hist,"Eta of the Top", "Eta", "Events/0.2", savename);
    TCanvas *c3 = drawPrint( im_top_hist,"Mass of Top", "Mass (GeV)", "Events/20GeV", savename);

    TCanvas *c4 = drawPrint( im_w_hist,"Mass of the W", "Mass (GeV)", "Events/2GeV", savename);    
    TCanvas *c5 = drawPrint( pt_w_hist,"Pt of the W", "Pt (GeV)", "Events/20GeV", savename);   
    TCanvas *c6 = drawPrint( eta_w_hist,"Eta of the W", "Eta", "Events/0.2", savename);

    TCanvas *c7 = drawPrint( pt_higgs_hist, "Pt of Higgs", "Pt (GeV)", "Events/20GeV", savename);
    TCanvas *c8 = drawPrint( eta_higgs_hist, "Eta of the Higgs", "Eta", "Events/0.2", savename);    
    TCanvas *c9 = drawPrint( im_higgs_hist, "Mass of the Higgs", "Mass (GeV)", "Events/20GeV", savename);

    TCanvas *c10 = drawPrint( pt_tprime_hist, "Pt of TPrime", "Pt (GeV)", "Events/20GeV", savename);
    TCanvas *c11 = drawPrint( eta_tprime_hist, "Eta of TPrime", "Eta", "Events/0.2", savename);
    TCanvas *c12 = drawPrint( im_tprime_hist, "Mass of TPrime", "Mass (GeV)", "Events/20GeV", savename);

    TCanvas *c13a = drawPrint( dr_tmax_hist, "Delat R of Top Max", "DeltaR", "Events/0.2", savename);
    TCanvas *c13b = drawPrint( dr_tmin_hist, "Delat R of Top Min", "DeltaR", "Events/0.2", savename);
    TCanvas *c13c = drawPrint( dr_wb_hist, "Delat R of wb Jet", "DeltaR", "Events/0.2", savename);
    TCanvas *c13d = drawPrintMulti( histlist, "Delat R of wb Jet, Min&Max Top Jet", "DeltaR", "Events/0.2", savename,0,1200, bob);
    TCanvas *c14 = drawPrint( dr_w_hist, "Delta R of W", "DeltaR", "Events/0.2" , savename);
    TCanvas *c15 = drawPrint( dr_higgs_hist, "Delta R of Higgs", "DeltaR", "Events/0.2" , savename);
    
    TCanvas *c16 = drawPrint2D(drvpt_top_hist, "DeltaR vs Pt of Top", "Pt (GeV)", "DeltaR", savename);
    TCanvas *c17 = drawPrint2D(drvpt_higgs_hist, "DeltaR vs Pt of Higgs", "Pt (GeV)", "DeltaR", savename);
    TCanvas *c18 = drawPrint2D(drvpt_w_hist, "DeltaR vs Pt of W", "Pt (GeV)", "DeltaR", savename);

    TCanvas *c19 = drawPrintLogX( im_qwt_hist, "Mass of qwt", "Mass (GeV)", "Events/LogGeV", savename);
    TCanvas *c20 = drawPrint( eta_qwt_hist, "Eta of qwt", "Eta", "Events/0.2", savename);
    TCanvas *c21 = drawPrint( pt_qwt_hist, "Pt  of qwt", "Pt (GeV)", "Events/10GeV" , savename);

    TCanvas *c22 = drawPrint( im_btT_hist, "Mass of btT", "Mass (GeV)", "Events/20GeV", savename);
    TCanvas *c23 = drawPrint( eta_btT_hist, "Eta of btT", "Eta", "Events/0.2", savename);
    TCanvas *c24 = drawPrint( pt_btT_hist, "Pt of btT", "Pt (GeV)", "Events/20GeV" , savename);

    TCanvas *c25 = drawPrint( im_bhT_hist, "Mass of bht", "Mass (GeV)", "Events/20GeV", savename);
    TCanvas *c26 = drawPrint( eta_bhT_hist, "Eta of bht", "Eta", "Events/0.2", savename);
    TCanvas *c27 = drawPrint( pt_bhT_hist, "Pt of bht", "Pt (GeV)", "Events/0.10GeV" , savename);

    TCanvas *c28 = drawPrint( pt_bg_hist, "Pt of b from Gluon","Pt (GeV)", "Events/20GeV" , savename);
    TCanvas *c29 = drawPrint( eta_bg_hist, "Eta of b from Gluon","Eta", "Events/0.2", savename );
    TCanvas *c30 = drawPrint( im_bg_hist, "Mass of b from Gluon","Mass (GeV)", "Events/0.1GeV", savename );
    
    TCanvas *c31 = drawPrint( pt_qq_hist, "Pt of q from qW","Pt (GeV)", "Events/20GeV", savename );
    TCanvas *c32 = drawPrint( eta_qq_hist, "Eta of q from qW","Eta", "Events/0.2" , savename);
    TCanvas *c33 = drawPrintLogX( im_qq_hist, "Mass of q from qW","Mass (GeV)", "Events/LogGeV", savename );

    TCanvas *c34 = drawPrintComb( pt_top_hist, dr_tmax_hist, drvpt_top_hist, "DeltaR vs Pt of Top Combo", savename ); 
    TCanvas *c35 = drawPrintComb( pt_higgs_hist, dr_higgs_hist, drvpt_higgs_hist, "DeltaR vs Pt of Higgs Combo", savename );
    TCanvas *c36 = drawPrintComb( pt_w_hist, dr_w_hist, drvpt_w_hist, "DeltaR vs Pt of W Combo", savename );
   
    TCanvas *c37 = drawPrint( Ht_hist, "Ht for Event Run","Ht (GeV)","Events/50GeV", savename);
    TCanvas *c38 = drawPrint( Heta_hist, "Total Eta for Event Run", "Eta", "Events/0.2", savename); 
    
    cout << "Writing log file" << endl;
    ofstream  out;
    string fileName = rootfilename.substr(0,stop ) + "/" +  rootfilename + cutname[cutrun] + "_" + "log.txt";
    char* temp = new char[ fileName.length() + 1 ];
    strcpy( temp, fileName.c_str() );
    out.open( temp );
	
    out << "Number of qw  found: " << qwcnt << endl;
    out << "Number of T   found: " << Tcnt << endl;/////
    out << "Number of tT  found: " << tTcnt << endl;/////
    out << "Number of bhT found: " << bhTcnt << endl;/////
    out << "Number of btT found: " << btTcnt << endl;/////
    out << "Number of hT  found: " << hTcnt << endl;//////
    out << "Number of wtT found: " << wtTcnt << endl;/////
    out << "Number of qq  found: " << qqcnt << endl;/////
    out << "Number of bg  found: " << bgcnt << endl;
    out << "Number of events found: " << evcnt << endl;
    out << "Number of events after Ht cut found: " << htcut << endl;
    out << "Number of events after Higgs cut found: " << higgscut << endl;
    out << "Number of events after top cut found: " << topcut << endl;
    out.close();

    cout << "End of cut" << endl;

 }// cut loop

 cout << "End of file : " << rootfilename << endl;
 return;

}//print histogram aboves


void runana( ){

	vector<string> rootname;
	rootname.push_back( myrootname0 );
	rootname.push_back( myrootname1 );
	rootname.push_back( myrootname2 );
	rootname.push_back( myrootname3 );
	rootname.push_back( myrootname4 );
	rootname.push_back( myrootname5 );
	rootname.push_back( myrootname6 );
	//rootname.push_back( myrootname7 );

   	for( int filenum = 0; filenum < rootname.size(); filenum++ ){
		cout << "Making Analysis Class" << endl;
		tprimeAnalisis tana( rootname[filenum]);
		cout << "Running Analysis" << endl;
		tana.Loop();	
	   }

   cout << "Thats all Folks!" << endl;
   return;
}

void make_eff_comp_hist()
{

	double eff_had[3][4] = { {1.880123895,		1.511145623,		1.836391328,		1.935738618 },
				    {0.54668811,		0.445212964,		0.467298692,		0.431720126 },
				    {0.080810557,		0.072946987,		0.071447565,		0.049200899 } };

	TH1D* eff_ht = new TH1D( "Efficiency Ratio vd Tprime Mass", "eff_vs_T_mass", 4, 1.1, 1.9 );
        TH1D* eff_hig = new TH1D( "Efficiency Ratio vd Tprime Mass", "eff_vs_T_mass", 4, 1.1, 1.9 );
        TH1D* eff_top = new TH1D( "Efficiency Ratio vd Tprime Mass", "eff_vs_T_mass", 4, 1.1, 1.9 );

	eff_ht->SetBinContent( 1, eff_had[0][0]/eff_had[0][0] );
	eff_ht->SetBinContent( 2, eff_had[0][1]/eff_had[0][0] );
	eff_ht->SetBinContent( 3, eff_had[0][2]/eff_had[0][0] );
	eff_ht->SetBinContent( 4, eff_had[0][3]/eff_had[0][0] );

        eff_hig->SetBinContent( 1, eff_had[1][0]/eff_had[1][0] );
        eff_hig->SetBinContent( 2, eff_had[1][1]/eff_had[1][0] );
        eff_hig->SetBinContent( 3, eff_had[1][2]/eff_had[1][0] );
        eff_hig->SetBinContent( 4, eff_had[1][3]/eff_had[1][0] );

        eff_top->SetBinContent( 1, eff_had[2][0]/eff_had[2][0] );
        eff_top->SetBinContent( 2, eff_had[2][1]/eff_had[2][0] );
        eff_top->SetBinContent( 3, eff_had[2][2]/eff_had[2][0] );
        eff_top->SetBinContent( 4, eff_had[2][3]/eff_had[2][0] );

	vector<TH1D*> hists;
	hists.push_back( eff_ht );
        hists.push_back( eff_hig );
        hists.push_back( eff_top );

	setTDRStyle();

	TCanvas *can = getTDRCanvas( "Cut_vs_tag_eff" );
        
        string xtit( "T prime Mass (GeV)");
        string ytit("Norilized Ratio of Tag/Cut Eff. / 0.2 GeV");

	double miny = 0.5;
	double maxy = 1.5;

	vector<string> label;
	label.push_back( "TAG/Cut Efficieny Ratio by Mass");
	label.push_back( "Ht Cut/Tag" );
	label.push_back( "Higs Cut/Tag" );
	label.push_back( "Top Cut.Tag" );
		

	tdrHistDrawMulti( hists, can, xtit, ytit, miny, maxy, label );
	
        string savetitle( "Ratio_Tag_Cut_Eff.png" );
        can->SaveAs( savetitle.c_str() );

}

void make_cvm_plots()
{

                int htc = 0;
		int hgc = 0;
		int tpc = 0;
		string temp("");
		double mass = 0;

        	TH1D* eff_ht = new TH1D( "Efficiency Ratio vd Tprime Mass", "eff_vs_T_mass", 4, 1.1, 1.9 );
        	TH1D* eff_hig = new TH1D( "Efficiency Ratio vd Tprime Mass", "eff_vs_T_mass", 4, 1.1, 1.9 );
        	TH1D* eff_top = new TH1D( "Efficiency Ratio vd Tprime Mass", "eff_vs_T_mass", 4, 1.1, 1.9 );

		std::ifstream maskFile;
                char maskFilePath[256];
                sprintf(maskFilePath, "Cuts_by_mass.txt");
                maskFile.open(maskFilePath, std::ifstream::in);
                if (!maskFile) {
                        std::cout << "ERROR: mask file <" << maskFilePath << "> can't be opened!"<<std::endl;
                }
                std::string line;
                std::vector< std::string > tokens;
                while(getline(maskFile, line)) {
                        if (line[0] != '#') {
                                std::stringstream ss(line);
                                std::string buf;
                                tokens.clear();
                                while (ss >> buf) {
                                        tokens.push_back(buf);
                                }
                              //  std::cout << "tok0 <" << tokens[0] << "> ";
                                if (tokens[0] == "For:" && tokens.size() >= 4) {
                                        htc = atoi(tokens[6].c_str());
                                        hgc = atoi(tokens[7].c_str());
                                        tpc = atoi(tokens[8].c_str());
					std::stringstream ss2( revrpsb( tokens[1] ));
					tokens.clear();
					while (ss2 >> buf ){
						tokens.push_back( buf );
						
					}
					mass = atoi( ((tokens[3].substr(1))).c_str() );	
                                //        std::cout << "mask pixel " << roc << " " << col << " " << row << std::endl;
                                }
                        }
                }
                maskFile.close();




}


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
