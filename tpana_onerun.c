#define tprimeAnalisis_cxx
#include "tpana_onerun.h"
#include <TCanvas.h>
#include "jwk_ku_tdr_style.h"
#include "CombineHistograms.h"

/*
static const bool bkgrd = false;
static const int cutlevel = 3;
static const bool dodr = true;
static const int ht_cut = 1000;
static const int top_cut = 350;
static const int hig_cut = 300;
*/

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
        		if( ( mother < 0 ) || ( abs( mom ) == mother ) ){ 
        			if( ( granny < 0 ) || ( abs( gran ) == granny )){
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
        TCanvas *c = getTDRCanvas( title + rfname , "combo");
	CombineHistograms( hist1, hist2, hist3, c );
        string savetitle( rfname + "_" + rpsb(title) + type );
        c->SaveAs( savetitle.c_str() );
        return c;
}

TCanvas* drawPrintComb( TH1D* hist1, TH1D* hist2, TH2D* hist3, string title, string rfname )
{
        return drawPrintComb( hist1, hist2, hist3, ".png", title, rfname );
}

bool  tprimeAnalisis::drcheck( double dr )
{
	bool result = false;
	if( dodr ){
		if( dr < 0.8 ) result = true;
	} else { result = true; }
	return result;
}

double  tprimeAnalisis::pt30_eta5_HtBounds( TLorentzVector &vec ) 
{
	double result = 0;
	if( ( vec.Pt() > 30 ) && ( abs( vec.Eta()) < 5 ) ) result = vec.Pt();
	return result; 
}

bool  tprimeAnalisis::ckhiggscut( TLorentzVector &hTvec, double dr_h )
{
	bool result = false;
	if( (hTvec.Pt() > hig_cut ) && ( abs( hTvec.Eta() ) < 2.4 ) && drcheck(dr_h)  ) result = true;
	return result;
}

bool  tprimeAnalisis::cktopcut( TLorentzVector &tTvec, double dr_wb ) //&& ( dr_wb < 0.8 )
{
	bool result = false;
	if( (tTvec.Pt() > top_cut ) && ( abs( tTvec.Eta() ) < 2.4 ) && drcheck( dr_wb )  ) result = true;
	return result;
}

bool  tprimeAnalisis::passcut( int cutrun, TLorentzVector &hTvec, double dr_h, TLorentzVector &tTvec, double dr_wb, double Ht, int htcut, int higgscut, int topcut ){

/*      if( Ht > 1100 ){ htcut++;
        if( (hTvec.Pt() > 300 ) && ( abs( hTvec.Eta() ) < 2.4 ) ){ higgscut++;
        if( (tTvec.Pt() > 400 ) && ( abs( tTvec.Eta() ) < 2.4 ) ){ topcut++;
*/
	bool result = false;

	if( cutrun == 0 ){ result = true; }
	else if( ( cutrun == 1 ) && ( Ht > ht_cut ) ) { htcut++; result = true; }
	else if( ( cutrun == 2 ) && ( Ht > ht_cut ) && ckhiggscut( hTvec, dr_h ) ){ higgscut++; result = true; }
        else if( ( cutrun == 3 ) && ( Ht > ht_cut ) && ckhiggscut( hTvec, dr_h ) && cktopcut( tTvec, dr_wb) ){ topcut++; result = true; }
//	else if( ( cutrun == 4 ) && ( Ht > ht_cut ) ){ xjetcut++; result = true; }
	
	return result;
}

void tprimeAnalisis::Loop()
{
   if (fChain == 0) return;

  string cutname[5] = { "_none", "_htcut", "_higgscut", "_topcut", "_xjetcut" };
  for( int cutrun = 0; cutrun < (cutlevel + 1); cutrun++ ){       ///  loop through cuts  <<<<<<<<<<<<<<<<<<<<<<<<<<    cutrun loop

    Long64_t nentries = fChain->GetEntriesFast();


/////////////////////////////////////////

    int dnc = -1;
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
    int enutrino = 12; 
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

    int qwcnt = 0;  //  quark from W bosun
    int Tcnt = 0; //  Tprime
    int tTcnt = 0; // top from a Tprime
    int tcnt = 0;  //  top ( Parantage not checked )
    int bhTcnt = 0; // bottom from higgs
    int btTcnt = 0; //  bottom from top
    int bcnt = 0; // bottom ( paratage not checked ) 
    int hTcnt = 0;  //  higgs from Tprime
    int wtTcnt = 0; // W from a top
    int qqcnt = 0; // quark from a quark
    int bgcnt = 0; // bottom from a gluon
    int evcnt = 0; //  event count
    int ltjetcnt = 0; //  ljet with W parent
    int otherjetcnt = 0; // ljet ( paratage not checked )

    double Ht = 0;
    int part = 0;

    int htcut = 0;
    int higgscut = 0;
    int topcut = 0;
    int xjetcut = 0;

    vector<TLorentzVector> qwjet;  //  q from w
    vector<TLorentzVector> bhTjet; //  b from h
    vector<TLorentzVector> qqgjet;  //  q from q 
    vector<TLorentzVector> tjet; 
    vector<TLorentzVector> bjet;
    vector<TLorentzVector> otherjet;     
 
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

    vector<int> ljet;
    ljet.push_back(down);
    ljet.push_back(up);
    ljet.push_back(strange);
    ljet.push_back(charm);
    ljet.push_back(electron);
    ljet.push_back(bottom);
  
    int ljetpartcnt[17] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};    
    int otherpartcnt[17] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    int allcnt[25] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    cout << "Init histograms :" << cutname[cutrun] << endl;

    vector<TH1D*> thlist;
    vector<TH2D*> th2list;
	


//bg        
    TH1D *pt_bg_hist = new TH1D("Pt of b from Gluon", "Pt of b from Gluon", 40, 0, 800);
    TH1D *im_bg_hist = new TH1D( "Mass of b from Gluon","Mass of b from Gluon", 200, 0, 20 );
    TH1D *eta_bg_hist = new TH1D("Eta of b from Gluon", "Eta of b from Gluon", 50, -5, 5);
    thlist.push_back( pt_bg_hist ); thlist.push_back( im_bg_hist );thlist.push_back( eta_bg_hist );
//qq        
    TH1D *pt_qq_hist = new TH1D("Pt of q from Quark", "Pt of q from Quark", 40, 0, 800);
    TH1D *im_qq_hist = new TH1D( "Mass of q from Quark","Mass of q from Quark", 200, 0, 2 );
    TH1D *eta_qq_hist = new TH1D("Eta of q from Quark", "Eta of q from Quark", 50, -5, 5);
    thlist.push_back( pt_qq_hist ); thlist.push_back( im_qq_hist );thlist.push_back( eta_qq_hist );
//tT        
    TH1D *pt_top_hist = new TH1D("Pt Top", "Pt Top",100, 0, 2000);
    TH1D *im_top_hist = new TH1D( "Mass Top","Mass of Top", 200, 0, 200 );
    TH1D *eta_top_hist = new TH1D("Eta Top", "Eta Top", 50, -5, 5);
    thlist.push_back( pt_top_hist ); thlist.push_back( im_top_hist );thlist.push_back( eta_top_hist );
//wt
    TH1D *pt_w_hist = new TH1D("Pt W", "Pt W", 100, 0, 2000);
    TH1D *eta_w_hist = new TH1D("Eta W", "Eta W", 50, -5, 5);    
    TH1D *im_w_hist = new TH1D( "Mass W","Mass of W", 200, 0, 200 );
    thlist.push_back( pt_w_hist ); thlist.push_back( im_w_hist );thlist.push_back( eta_w_hist );
//hT
    TH1D *pt_higgs_hist = new TH1D("Pt Higgs","Pt Higgs", 80, 0, 1600 );
    TH1D *eta_higgs_hist = new TH1D("Eta of Higgs", "Eta of HIggs", 50, -5,5 );
    TH1D *im_higgs_hist = new TH1D("Mass Higgs","Mass of the Higgs", 200, 0, 200 );
    thlist.push_back( pt_higgs_hist ); thlist.push_back( im_higgs_hist );thlist.push_back( eta_higgs_hist );
//T
    TH1D *im_tprime_hist = new TH1D( "Mass TPrime","Mass tPrime", 70, 1400, 2800 );
    TH1D *pt_tprime_hist = new TH1D( "Pt TPrime","Pt of tPrime", 100, 1000, 3000 );    
    TH1D *eta_tprime_hist = new TH1D( "Eta TPrime","Eta of tPrime", 50, -5, 5 );
    thlist.push_back( pt_tprime_hist ); thlist.push_back( im_tprime_hist );thlist.push_back( eta_tprime_hist );
//qwt   <<<<<<<<<<<<
    TH1D *pt_qwt_hist = new TH1D("Pt qw Jet", "Pt qw Jet", 100, 0, 2000);
    TH1D *im_qwt_hist = new TH1D( "Mass qw Jet","Mass of qw Jet", 200, 0, 2 );
    TH1D *eta_qwt_hist = new TH1D("Eta qw Jet", "Eta qw Jet", 50, -5, 5);
    thlist.push_back( pt_qwt_hist ); thlist.push_back( im_qwt_hist );thlist.push_back( eta_qwt_hist );
//btT
    TH1D *pt_btT_hist = new TH1D("Pt btT", "Pt btT", 40, 0, 800);
    TH1D *im_btT_hist = new TH1D( "Mass btT","Mass of btT", 100, 0, 100 );
    TH1D *eta_btT_hist = new TH1D("Eta btT", "Eta btT", 50, -5, 5);
    thlist.push_back( pt_btT_hist ); thlist.push_back( im_btT_hist );thlist.push_back( eta_btT_hist );
//bhT
    TH1D *pt_bhT_hist = new TH1D("Pt bh Jet", "Pt bh Jet", 50, 0, 1000);
    TH1D *im_bhT_hist = new TH1D( "Mass bh Jet","Mass of bh Jet", 100, 0, 100 );
    TH1D *eta_bhT_hist = new TH1D("Eta bh Jet", "Eta bh Jet", 50, -5, 5);
    thlist.push_back( pt_bhT_hist ); thlist.push_back( im_bhT_hist );thlist.push_back( eta_bhT_hist );
//dr
    TH1D *dr_tmin_hist = new TH1D( "Min DeltaR Top", "Min DeltaR of Top", 50, 0, 5 );
    TH1D *dr_tmax_hist = new TH1D( "Max DeltaR Top", "Max DeltaR of Top", 50, 0, 5 );
    TH1D *dr_higgs_hist = new TH1D( "DelaR Higgs", "DeltaR Higgs", 50, 0, 5 );
    TH1D *dr_w_hist = new TH1D( "DeltaR W+", "DeltaR of W+", 50, 0, 5 );
    TH1D *dr_wb_hist = new TH1D( "DeltaR of wb Jet", "DeltaR of wb Jet", 50, 0, 5 );
    thlist.push_back( dr_tmin_hist ); thlist.push_back( dr_tmax_hist );thlist.push_back( dr_higgs_hist ); thlist.push_back( dr_w_hist ); thlist.push_back( dr_wb_hist );
//dr vs Pt
    TH2D *drvpt_top_hist = new TH2D( "DeltaR vs Pt of Top", "DeltaR vs Pt of Top", 800, 0, 1600, 50, 0, 5 );
    TH2D *drvpt_higgs_hist = new TH2D( "DeltaR vs Pt of Higgs", "DeltaR vs Pt of Higgs", 80, 0, 1600, 50, 0, 5 );
    TH2D *drvpt_w_hist = new TH2D( "DeltaR vs Pt of W", "DeltaR vs Pt of W", 80, 0, 1600, 50, 0, 5 );
    th2list.push_back( drvpt_top_hist ); th2list.push_back( drvpt_higgs_hist );th2list.push_back( drvpt_w_hist );
//Ht
    TH1D *Ht_hist = new TH1D( "Ht", "Ht", 120, 0, 2400 );
    TH1D *Heta_hist = new TH1D( "All Jet Eta", "All Jet Eta", 50, -5, 5 );
    TH1D *Hpt_hist = new TH1D( "All Jet Pt", "All Jet Pt", 50, 0, 1000 );
    thlist.push_back( Ht_hist ); thlist.push_back( Heta_hist ); thlist.push_back( Hpt_hist );
//Pt bs Eta
    TH2D *ptvseta_qqg_hist = new TH2D( "Lt Quark from Gluon Jet", "Lt Quark from Gluon Jet", 80, 0, 1600, 50, -5, 5 );
    TH2D *ptvseta_bhT_hist = new TH2D( "Bottom from Higgs Jet", "Bottom from Higgs Jet", 80, 0, 1600, 50, -5, 5 );
    TH2D *ptvseta_bg_hist = new TH2D( "Bottom from Gluon Jet", "Bottom from Gluon Jet", 80, 0, 1600, 50, -5, 5 );
    TH2D *ptvseta_qwt_hist = new TH2D( "Lt Quark from W Jet", "Lt Quark from W Jet", 80, 0, 1600, 50, -5, 5 );
    TH2D *ptvseta_btT_hist = new TH2D( "Bottom from Top Jet", "Bottom from Top Jet", 80, 0, 1600, 50, -5, 5 );
    th2list.push_back(ptvseta_qqg_hist); th2list.push_back(ptvseta_bhT_hist); th2list.push_back(ptvseta_bg_hist); 
    th2list.push_back(ptvseta_qwt_hist); th2list.push_back(ptvseta_btT_hist);


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
	//TLorentzVector qqvec;

        TLorentzVector vReturn;
        int charg = 0;

        //Loop Through Particles
//	cout << "Looping Particles: Number of particles in the event : " << *Event_Nparticles << endl;
        for(int partnum=0; partnum < *Event_Nparticles; partnum++){ // partnum = particle number
  //          cout << " in loop with particle number :" << partnum << endl;
	    vReturn.SetPtEtaPhiM(0,0,0,0);
	    charg = 0;
            if( get2VectorFromLists( partnum, bglist, qglist, dnc, vReturn, charg, part )){ bgcnt++; bgvec = vReturn;}
            if( get2VectorFromLists( partnum, qglist, qglist, dnc, vReturn, charg, part )){ qqcnt++; qqgjet.push_back( vReturn ); }
            if( get2Vector( partnum, top, tprime, dnc, vReturn, charg )){ tTcnt++; tTvec = vReturn;}
	    if( get2Vector( partnum, top, dnc, dnc, vReturn, charg )){ allcnt[top]++; tjet.push_back( vReturn ); }
            if( get2Vector( partnum, w, top, dnc, vReturn, charg )) { wtTcnt++; wtTvec = vReturn; }
            if( get2Vector( partnum, w, dnc, dnc, vReturn, charg )) { allcnt[w]++; }            
	    if( get2Vector( partnum, higgs, tprime, dnc, vReturn, charg )) { allcnt[higgs]++; hTcnt++; hTvec = vReturn;}
	    if( get2Vector( partnum, tprime, dnc, dnc, vReturn, charg )) { allcnt[8]++; Tcnt++; Tvec = vReturn;}
            if( get2Vector( partnum, bottom, top, dnc, vReturn, charg )) { btTcnt++; btTvec = vReturn;}
            if( get2Vector( partnum, bottom, higgs, dnc, vReturn, charg )){ bhTcnt++; bhTjet.push_back( vReturn ); } 
	    if( get2VectorFromJet( partnum, ljet, w, dnc, vReturn, charg, part )){ allcnt[bottom]++; qwcnt++; ljetpartcnt[part]++; qwjet.push_back( vReturn );}
	    if( get2Vector( partnum, bottom, dnc, dnc, vReturn, charg )){ allcnt[bottom]++; bcnt++; bjet.push_back( vReturn );}            
	    if( get2VectorFromLists( partnum, ljet, qglist, dnc, vReturn, charg, part )){allcnt[part]++; otherpartcnt[part]++; otherjetcnt++; otherjet.push_back( vReturn );}

        }// find particles in above loop

// calc variblesi

	if( bkgrd ){
		if( ( rand()%10 ) > 4 ) { hTvec = tjet[0]; Tvec = tjet[1]; }
		else { hTvec = tjet[1]; Tvec = tjet[0]; }
	}
	
//	cout << "Calcating values for event" << endl;
	vReturn.SetPtEtaPhiM(0,0,0,0);
	evcnt++;
	// delta_R_top
	double dr1 = 0;
	double dr2 = 0;
	double dr3 = 0;
        double dr_tmax = 0;
        double dr_tmin = 0;
	double dr_wb = 0;
	double dr_h = 0;
	double dr_w = 0;	

	if( qwjet.size() < 1 ){ qwjet.push_back( vReturn); qwjet.push_back( vReturn);cout << "qwjet empty" << endl; }
	if( qqgjet.size() < 1 ){ qqgjet.push_back( vReturn); qqgjet.push_back( vReturn);cout << "qqgjet empty" << endl; }
	if( bhTjet.size() < 1 ){ bhTjet.push_back( vReturn); bhTjet.push_back( vReturn);cout << "bhTjet empty" << endl; }

	dr1 = (qwjet[0]).DeltaR( qwjet[1] );
	dr2 = (qwjet[1]).DeltaR( btTvec );
	dr3 = (qwjet[0]).DeltaR( btTvec );
       
//    	cout << "stop1"<<endl;
        dr_wb = wtTvec.DeltaR( btTvec );
	if( dr1 > dr2 ){ if( dr1 > dr3 ){ dr_tmax = dr1;} else { dr_tmax = dr3;}} else { if( dr2 > dr3 ){ dr_tmax = dr2;} else{ dr_tmax = dr3; }}	
	if( dr1 < dr2 ){ if( dr1 < dr3 ){ dr_tmin = dr1;} else { dr_tmin = dr3;}} else { if( dr2 < dr3 ){ dr_tmin = dr2;} else{ dr_tmin = dr3; }} 
	dr_h = (bhTjet[0]).DeltaR( bhTjet[1] );
   	dr_w = dr1;
//pt>30 and eta < 5 for inclusion to Ht
	Ht = pt30_eta5_HtBounds( bgvec ) + pt30_eta5_HtBounds( qqgjet[0] ) + pt30_eta5_HtBounds( btTvec );
	Ht = Ht  + pt30_eta5_HtBounds(bhTjet[0]) + pt30_eta5_HtBounds(bhTjet[1])+ pt30_eta5_HtBounds(qwjet[0]) + pt30_eta5_HtBounds(qwjet[1]);
//	Heta = bgvec.Eta() + (qqgjet[0]).Eta() + btTvec.Eta() + (bhTjet[0]).Eta() + (bhTjet[1]).Eta()+ (qwjet[0]).Eta() + (qwjet[1]).Eta();

//	Fill cutts<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	
//	cout << "Checking Cut for cut level :" << cutrun << endl;

	if( passcut( cutrun, hTvec, dr_h,  tTvec, dr_wb, Ht, htcut, higgscut, topcut  ) ) {
//	cout << "Filling histograms" << endl;  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

//      fills
	if( bgvec.Pt() > 0 ){
        pt_bg_hist->Fill( bgvec.Pt() );
        eta_bg_hist->Fill( bgvec.Eta() );
        im_bg_hist->Fill( bgvec.M() );
	Heta_hist->Fill( bgvec.Eta() );
	Hpt_hist->Fill( bgvec.Pt() );
	}

	for( int i = 0; i < qqgjet.size(); i++ ){
        pt_qq_hist->Fill( (qqgjet[i]).Pt() );
        eta_qq_hist->Fill( (qqgjet[i]).Eta() );
        im_qq_hist->Fill( (qqgjet[i]).M() );
        Heta_hist->Fill( (qqgjet[i]).Eta() );
        Hpt_hist->Fill( (qqgjet[i]).Pt() );
	}

	if( tTvec.Pt() > 0 ){
        pt_top_hist->Fill( tTvec.Pt() );
        eta_top_hist->Fill( tTvec.Eta() );
        im_top_hist->Fill( tTvec.M() );
        }

        if(wtTvec.Pt() > 0 ){
        pt_w_hist->Fill( wtTvec.Pt() );
        eta_w_hist->Fill( wtTvec.Eta() );
        im_w_hist->Fill( wtTvec.M() );
	}

	if( hTvec.Pt() > 0 ){
	pt_higgs_hist->Fill( hTvec.Pt() );
        eta_higgs_hist->Fill( hTvec.Eta() );
  	im_higgs_hist->Fill( hTvec.M() );
	}

	if( Tvec.Pt() > 0 ){
	pt_tprime_hist->Fill( Tvec.Pt() );
        eta_tprime_hist->Fill( Tvec.Eta() );
        im_tprime_hist->Fill( Tvec.M() );
	}

	for( int c = 0; c <  qwjet.size(); c++ ){
        pt_qwt_hist->Fill( (qwjet[c]).Pt() );
        eta_qwt_hist->Fill( (qwjet[c]).Eta() );
        im_qwt_hist->Fill( (qwjet[c]).M() );
	}
	
	for( int c = 0; c < bhTjet.size(); c++ ){
        pt_bhT_hist->Fill( (bhTjet[c]).Pt() );
        eta_bhT_hist->Fill( (bhTjet[c]).Eta() );
        im_bhT_hist->Fill( (bhTjet[c]).M() );
	}

	if( btTvec.Pt() > 0 ){
        pt_btT_hist->Fill( btTvec.Pt() );
        eta_btT_hist->Fill( btTvec.Eta() );
        im_btT_hist->Fill( btTvec.M() );
	}
	
	if( dr_tmin > 0 ) dr_tmin_hist->Fill( dr_tmin );
	if( dr_tmin > 0 ) dr_tmax_hist->Fill( dr_tmax );
	if( dr_wb > 0 ) dr_wb_hist->Fill( dr_wb );
        if( dr_w > 0 ) dr_w_hist->Fill( dr_w );
	if( dr_h > 0 ) dr_higgs_hist->Fill( dr_h );

        if( (tTvec.Pt() > 0) && ( dr_tmax > 0 ) ) drvpt_top_hist->Fill( tTvec.Pt(), dr_tmax);
	if( (hTvec.Pt() > 0) && ( dr_h > 0 ) ) drvpt_higgs_hist->Fill( hTvec.Pt(), dr_h );
        if( (wtTvec.Pt() > 0) && ( dr_w > 0 ) ) drvpt_w_hist->Fill( wtTvec.Pt(), dr_w );

	for( int c = 0; c < qqgjet.size(); c++ ){ ptvseta_qqg_hist->Fill( (qqgjet[c]).Pt(), (qqgjet[c]).Eta() );}
        if( btTvec.Pt() > 0 ) ptvseta_btT_hist->Fill( btTvec.Pt(), btTvec.Eta() );
        if( bgvec.Pt() > 0 ) ptvseta_bg_hist->Fill( bgvec.Pt(),  bgvec.Eta() );
	
	for( int c = 0; c < qwjet.size(); c++ ){ ptvseta_qwt_hist->Fill( qwjet[c].Pt(),  qwjet[c].Eta() );}
        for( int c = 0; c < bhTjet.size(); c++ ){ ptvseta_bhT_hist->Fill( bhTjet[c].Pt(),  bhTjet[c].Eta() );}	

	Ht_hist->Fill( Ht );
//	Heta_hist->Fill( Heta );
	

 	}//<<<<<<  cuts  if pass cuts

        qwjet.clear();
        bhTjet.clear();
        qqgjet.clear();
        tjet.clear();
        bjet.clear();
        otherjet.clear();
       
//        cout << "Event Finished" << endl;   
    }// calc values and fill histograms
    cout << "Print and Draw Histograms" << endl;  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    vector<TH1D*> histlist;
    	histlist.push_back( dr_tmax_hist );
    	histlist.push_back( dr_tmin_hist );
    	histlist.push_back( dr_wb_hist );
    
    vector<TH1D*> histjetptlist;
    	histjetptlist.push_back(pt_bg_hist);
	histjetptlist.push_back(pt_qq_hist);
	histjetptlist.push_back(pt_btT_hist);
	histjetptlist.push_back(pt_bhT_hist);
	histjetptlist.push_back(pt_qwt_hist);

    vector<TH1D*> histjetetalist;
        histjetetalist.push_back(eta_bg_hist);
        histjetetalist.push_back(eta_qq_hist);
        histjetetalist.push_back(eta_btT_hist);
        histjetetalist.push_back(eta_bhT_hist);
        histjetetalist.push_back(eta_qwt_hist);
	
    vector<string> bob;
    	bob.push_back( "Legend" );
    	bob.push_back( "DeltaR of Tmax" );
    	bob.push_back( "DeltaR of Tmin" );
    	bob.push_back( "DeltaR of Wb" );
        bob.push_back( "DeltaR of Tmin" );
        bob.push_back( "DeltaR of Wb" );
        bob.push_back( "DeltaR of Tmin" );
        bob.push_back( "DeltaR of Wb" );

   vector<string> jim;
        jim.push_back( "Legend" );
        jim.push_back( "Bottom from Gluon" );
        jim.push_back( "Lt. Quark of Lt. Quark" );
        jim.push_back( "Bottom from Top from TPrime" );
        jim.push_back( "Bottom from Higgs from Tprime" );
        jim.push_back( "Quark from W" );

    setTDRStyle();
    string savename = rootfilename + "/" +  rootfilename + cutname[cutrun];
    cout << savename << endl;

    vector<TCanvas*> canlist;

    TCanvas *c1 = drawPrint( pt_top_hist,"Pt of the Top", "Pt (GeV)", "Events/20GeV", savename );   
    TCanvas *c2 = drawPrint( eta_top_hist,"Eta of the Top", "Eta", "Events/0.2", savename);
    TCanvas *c3 = drawPrint( im_top_hist,"Mass of Top", "Mass (GeV)", "Events/20GeV", savename);
    canlist.push_back( c1 );canlist.push_back( c2 );canlist.push_back( c3 );

    TCanvas *c4 = drawPrint( im_w_hist,"Mass of the W", "Mass (GeV)", "Events/2GeV", savename);    
    TCanvas *c5 = drawPrint( pt_w_hist,"Pt of the W", "Pt (GeV)", "Events/20GeV", savename);   
    TCanvas *c6 = drawPrint( eta_w_hist,"Eta of the W", "Eta", "Events/0.2", savename);
    canlist.push_back( c4 );canlist.push_back( c5 );canlist.push_back( c6 );

    TCanvas *c7 = drawPrint( pt_higgs_hist, "Pt of Higgs", "Pt (GeV)", "Events/20GeV", savename);
    TCanvas *c8 = drawPrint( eta_higgs_hist, "Eta of the Higgs", "Eta", "Events/0.2", savename);    
    TCanvas *c9 = drawPrint( im_higgs_hist, "Mass of the Higgs", "Mass (GeV)", "Events/20GeV", savename);
    canlist.push_back( c7 );canlist.push_back( c8 );canlist.push_back( c9 );

    TCanvas *c10 = drawPrint( pt_tprime_hist, "Pt of TPrime", "Pt (GeV)", "Events/20GeV", savename);
    TCanvas *c11 = drawPrint( eta_tprime_hist, "Eta of TPrime", "Eta", "Events/0.2", savename);
    TCanvas *c12 = drawPrint( im_tprime_hist, "Mass of TPrime", "Mass (GeV)", "Events/20GeV", savename);
    canlist.push_back( c10 );canlist.push_back( c11 );canlist.push_back( c12 );

    TCanvas *c13a = drawPrint( dr_tmax_hist, "Delat R of Top Max", "DeltaR", "Events/0.2", savename);
    TCanvas *c13b = drawPrint( dr_tmin_hist, "Delat R of Top Min", "DeltaR", "Events/0.2", savename);
    TCanvas *c13c = drawPrint( dr_wb_hist, "Delat R of wb Jet", "DeltaR", "Events/0.2", savename);
    TCanvas *c13d = drawPrintMulti( histlist, "Delat R of wb Jet, Min&Max Top Jet", "DeltaR", "Events/0.2", savename,0,1200, bob);
    TCanvas *c14 = drawPrint( dr_w_hist, "Delta R of W", "DeltaR", "Events/0.2" , savename);
    TCanvas *c15 = drawPrint( dr_higgs_hist, "Delta R of Higgs", "DeltaR", "Events/0.2" , savename);
    canlist.push_back( c13a );canlist.push_back( c13b );canlist.push_back( c13c );
    canlist.push_back( c13d );canlist.push_back( c14 );canlist.push_back( c15 );
    
    TCanvas *c16 = drawPrint2D(drvpt_top_hist, "DeltaR vs Pt of Top", "Pt (GeV)", "DeltaR", savename);
    TCanvas *c17 = drawPrint2D(drvpt_higgs_hist, "DeltaR vs Pt of Higgs", "Pt (GeV)", "DeltaR", savename);
    TCanvas *c18 = drawPrint2D(drvpt_w_hist, "DeltaR vs Pt of W", "Pt (GeV)", "DeltaR", savename);
    canlist.push_back( c16 );canlist.push_back( c17 );canlist.push_back( c18 );

    TCanvas *c19 = drawPrint( im_qwt_hist, "Mass of qw jet", "Mass (GeV)", "Events/LogGeV", savename);
    TCanvas *c20 = drawPrint( eta_qwt_hist, "Eta of qw jet", "Eta", "Events/0.2", savename);
    TCanvas *c21 = drawPrint( pt_qwt_hist, "Pt  of qw jet", "Pt (GeV)", "Events/10GeV" , savename);
    canlist.push_back( c19 );canlist.push_back( c20 );canlist.push_back( c21 );

    TCanvas *c22 = drawPrint( im_btT_hist, "Mass of btT", "Mass (GeV)", "Events/20GeV", savename);
    TCanvas *c23 = drawPrint( eta_btT_hist, "Eta of btT", "Eta", "Events/0.2", savename);
    TCanvas *c24 = drawPrint( pt_btT_hist, "Pt of btT", "Pt (GeV)", "Events/20GeV" , savename);
    canlist.push_back( c22 );canlist.push_back( c23 );canlist.push_back( c24 );

    TCanvas *c25 = drawPrint( im_bhT_hist, "Mass of bht", "Mass (GeV)", "Events/20GeV", savename);
    TCanvas *c26 = drawPrint( eta_bhT_hist, "Eta of bht", "Eta", "Events/0.2", savename);
    TCanvas *c27 = drawPrint( pt_bhT_hist, "Pt of bht", "Pt (GeV)", "Events/0.10GeV" , savename);
    canlist.push_back( c25 );canlist.push_back( c26 );canlist.push_back( c27 );

    TCanvas *c28 = drawPrint( pt_bg_hist, "Pt of b from Gluon","Pt (GeV)", "Events/20GeV" , savename);
    TCanvas *c29 = drawPrint( eta_bg_hist, "Eta of b from Gluon","Eta", "Events/0.2", savename );
    TCanvas *c30 = drawPrint( im_bg_hist, "Mass of b from Gluon","Mass (GeV)", "Events/0.1GeV", savename );
    canlist.push_back( c28 );canlist.push_back( c29 );canlist.push_back( c30 );
    
    TCanvas *c31 = drawPrint( pt_qq_hist, "Pt of q from qW","Pt (GeV)", "Events/20GeV", savename );
    TCanvas *c32 = drawPrint( eta_qq_hist, "Eta of q from qW","Eta", "Events/0.2" , savename);
    TCanvas *c33 = drawPrint( im_qq_hist, "Mass of q from qW","Mass (GeV)", "Events/LogGeV", savename );
    canlist.push_back( c31 );canlist.push_back( c32 );canlist.push_back( c33 );

    TCanvas *c34 = drawPrintComb( pt_top_hist, dr_tmax_hist, drvpt_top_hist, "DeltaR vs Pt of Top Combo", savename ); 
    TCanvas *c35 = drawPrintComb( pt_higgs_hist, dr_higgs_hist, drvpt_higgs_hist, "DeltaR vs Pt of Higgs Combo", savename );
    TCanvas *c36 = drawPrintComb( pt_w_hist, dr_w_hist, drvpt_w_hist, "DeltaR vs Pt of W Combo", savename );
    canlist.push_back( c34 );canlist.push_back( c35 );canlist.push_back( c36 );
   
    TCanvas *c37 = drawPrint( Ht_hist, "Ht for Event Run","Ht (GeV)","Events/50GeV", savename);
    TCanvas *c38 = drawPrint( Heta_hist, "Total Eta for Event Run", "Eta", "Events/0.2", savename);
    TCanvas *c39 = drawPrintMulti( histjetptlist, "Pt of Light Jets", "Pt(GeV)", "Events/20(GeV)", savename,0,5000, jim); 
    TCanvas *c40 = drawPrintMulti( histjetetalist, "Eta of Light Jets", "Eta", "Events/0.2", savename,0,5000, jim);   
    TCanvas *c46 = drawPrint( Hpt_hist, "Total Pt for Event Run", "Pt(Gev)", "Events/20(GeV)", savename);
    canlist.push_back( c37 );canlist.push_back( c38 );canlist.push_back( c39 );canlist.push_back( c40 );canlist.push_back( c46 );
    
    TCanvas *c41 = drawPrint2D(ptvseta_qqg_hist, "Pt vs Eta of Lt Quark From Gluon Jet", "Pt (GeV)", "DeltaR", savename);
    TCanvas *c42 = drawPrint2D(ptvseta_bhT_hist, "Pt vs Eta of Bottom from Higgs Jet", "Pt (GeV)", "DeltaR", savename);
    TCanvas *c43 = drawPrint2D(ptvseta_bg_hist, "Pt vs Eta of Bottom from Gluon Jet", "Pt (GeV)", "DeltaR", savename);
    TCanvas *c44 = drawPrint2D(ptvseta_qwt_hist, "Pt vs Eta of Quark from W Jet", "Pt (GeV)", "DeltaR", savename);
    TCanvas *c45 = drawPrint2D(ptvseta_btT_hist, "Pt vs Eta of Bottom from Top Jet", "Pt (GeV)", "DeltaR", savename);
    canlist.push_back( c41 );canlist.push_back( c42 );canlist.push_back( c43 );canlist.push_back( c44 );canlist.push_back( c45 );

    cout << "Writing log file" << endl;  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ofstream  out;
    string fileName = savename + "_" + "log.txt";
    char* temp = new char[ fileName.length() + 1 ];
    strcpy( temp, fileName.c_str() );
    out.open( temp );

/*	
    int none = 0;
    int down = 1;
    int up = 2;
    int strange = 3;
    int charm = 4;otherjetcnt
    int bottom = 5;
    int top = 6;
    int bprime = 7;
    int tprime = 8000001;

    int electron = 11;
    int enutrino = 12; 
    int muon = 13;
    int mnutrino = 14;
    int tau = 15;
    int tnutrino = 16;
*/
    out << "For File: " << savename << endl;
    out << "Number of quarks from W  found: " << qwcnt << endl;
    out << "Number of Tprime   found: " << Tcnt << endl;/////
    out << "Number of top from Tprime found: " << tTcnt << endl;/////
    out << "Number of bottom from higs found: " << bhTcnt << endl;/////
    out << "Number of bottom from top found: " << btTcnt << endl;/////
    out << "Number of higgs  found: " << hTcnt << endl;//////
    out << "Number of  W from top  found: " << wtTcnt << endl;/////
    out << "Number of quark from quark  found: " << qqcnt << endl;/////
    out << "Number of bottom from gluon  found: " << bgcnt << endl;
    out << "Number of lt. jets from W found: " << qwcnt << endl;
    out << "Number and Type of particles in Light Jets: Down: "<< ljetpartcnt[1] <<  " Up: " << ljetpartcnt[2] << " Strange: " 
		<< ljetpartcnt[3] <<  " Charm: " << ljetpartcnt[4] <<  " Bottom: " << ljetpartcnt[5] << " Electron: " << ljetpartcnt[11] << endl;
    out << "Number of bottom overall : " << bcnt << endl;
    out << "Number of tops overall : " << allcnt[top] << endl;
    out << "Number of Other Lt. jets : " << otherjetcnt << endl;
    out << "Number and Type of particles in Other Jets: Down: "<< otherpartcnt[1] <<  " Up: " << otherpartcnt[2] << " Strange: " 
		<< otherpartcnt[3] <<  " Charm: " << otherpartcnt[4] <<  " Bottom: " << otherpartcnt[5] << " Electron: " << otherpartcnt[11] << endl;
    out << "Number of events found: " << evcnt << endl;    
    out << "Number and Type of particles in Event: Down: "<< allcnt[1] << " Up: " << allcnt[2] << " Strange: " << allcnt[3] <<  " Charm: " 
		<< allcnt[4] <<  " Bottom: " << allcnt[5] << " Top : " << allcnt[6] << " Tprime : " << allcnt[8] << " W : " << allcnt[w] << " Electron: " << allcnt[11] << endl;
    out << "Number of events after Ht cut found: " << htcut << endl;
    out << "Number of events after Higgs cut found: " << higgscut << endl;
    out << "Number of events after top cut found: " << topcut << endl;

    out.close();

    cout << "End of cut" << endl;

    cout << "Cleaning Up" << endl;
    cout << "Deleteing Canvas" << endl;
    for( int cl = 0; cl < canlist.size(); cl++ ){ delete canlist[cl]; }	
    cout << "Deleting TH1D" << endl;
    for( int cl = 0; cl < thlist.size(); cl++ ){ delete thlist[cl]; }    
    cout << "Deleteing TH2Ds" << endl;
    for( int cl = 0; cl < th2list.size(); cl++ ){ delete th2list[cl]; } 
    cout << "Clearing lsits" << endl;
    canlist.clear();
    thlist.clear();
    th2list.clear();

    cout << "Clean Up Finished" << endl;

 }// cut loop
   
 cout << "End of file : " << rootfilename << endl;
 return;

}//print histogram aboves

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  loops tprime analysis over masses and cuts  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

/*
static const bool bkgrd = false;
static const int cutlevel = 3;
static const bool dodr = true;
static const int ht_cut = 1000;
static const int top_cut = 350;
static const int hig_cut = 300;
*/

void runana( ){

	vector<string> rootname;
	rootname.push_back( myrootname0 );
	rootname.push_back( myrootname1 );
	rootname.push_back( myrootname2 );
	rootname.push_back( myrootname3 );
	rootname.push_back( myrootname4 );
	rootname.push_back( myrootname5 );
	rootname.push_back( myrootname6 );
	rootname.push_back( myrootname7 );

//   loop on cut types - add internal ref and cut varibles


   	for( int filenum = 0; filenum < rootname.size(); filenum++ ){
		cout << "Making Analysis Class " << rootname[filenum] << endl;
		tprimeAnalisis *tana = new tprimeAnalisis( rootname[filenum] );
		cout << "Running Analysis" << endl;
		tana->Loop();
		delete tana;	
	   }

//    end cut tyoes loop ( 16 flavors )

   cout << "Thats all Folks!" << endl;
   return;
}

//  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  make efficiency comparisons plot  <<<<<<<<<<<<<<<<<<<<<<<<<

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

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  creates Figure of Signficance plots  <<<<<<<<<<<<<<<<<<<<<<<<<<

void make_cvm_plots()
{

                int htc = 0;
		int hgc = 0;
		int tpc = 0;
		string temp("");
		double mass = 0;

//        	TH1D* eff_ht = new TH1D( "Efficiency Ratio vs Tprime Mass", "eff_vs_T_mass", 4, 1.1, 1.9 );
//        	TH1D* eff_hig = new TH1D( "Efficiency Ratio vs Tprime Mass", "eff_vs_T_mass", 4, 1.1, 1.9 );
//        	TH1D* eff_top = new TH1D( "Efficiency Ratio vs Tprime Mass", "eff_vs_T_mass", 4, 1.1, 1.9 );

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
