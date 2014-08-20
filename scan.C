#include "TTreeCache.h"
#include "TChain.h"
#include "TH1.h"
#include "TString.h"
#include "Math/VectorUtil.h"
#include "TDatabasePDG.h"

#include <iostream>
#include <math.h>
#include <utility>

#include "./CORE/CMS2.h"
// #include "./CORE/jetSelections.h"
#include "/home/users/namin/macros/utils.C"
#include "./particles.C"

using namespace std;
using namespace tas;


struct JetStruct {
    LorentzVector jet;
    float mindRToQuark;
    unsigned int idx;
};

void scan(float jetPtCut=40.0, TString tag = "") {

    TChain *ch = new TChain("Events");

    // ch->Add("/hadoop/cms/store/user/namin/mStop825to900_mLSP1/*.root");

    ch->Add("/hadoop/cms/store/user/namin/mStop825to900_mLSP1/ntuple_309_4_16a.root");
    ch->Add("/hadoop/cms/store/user/namin/mStop825to900_mLSP1/ntuple_253_1_f2L.root");
    ch->Add("/hadoop/cms/store/user/namin/mStop825to900_mLSP1/ntuple_343_3_l8P.root");
    ch->Add("/hadoop/cms/store/user/namin/mStop825to900_mLSP1/ntuple_329_3_4dB.root");
    ch->Add("/hadoop/cms/store/user/namin/mStop825to900_mLSP1/ntuple_336_1_yfv.root");
    ch->Add("/hadoop/cms/store/user/namin/mStop825to900_mLSP1/ntuple_358_1_IPX.root");
    ch->Add("/hadoop/cms/store/user/namin/mStop825to900_mLSP1/ntuple_346_1_ptQ.root");
    ch->Add("/hadoop/cms/store/user/namin/mStop825to900_mLSP1/ntuple_281_1_Qqp.root");
    ch->Add("/hadoop/cms/store/user/namin/mStop825to900_mLSP1/ntuple_362_1_NNy.root");
    ch->Add("/hadoop/cms/store/user/namin/mStop825to900_mLSP1/ntuple_361_1_Yhs.root");
    ch->Add("/hadoop/cms/store/user/namin/mStop825to900_mLSP1/ntuple_359_1_YrF.root");
    ch->Add("/hadoop/cms/store/user/namin/mStop825to900_mLSP1/ntuple_360_1_sfh.root");


    ch->Add("/hadoop/cms/store/user/namin/mStop150to475_mLSP1/*.root");
    ch->Add("/hadoop/cms/store/user/namin/mStop500to800_mLSP1/*.root");

    int nEventsTotal = ch->GetEntries();
    int nEventsSoFar = 0;
    int nGoodEvents = 0;

    TFile *currentFile = 0;
    TObjArray *listOfFiles = ch->GetListOfFiles();
    TIter fileIter(listOfFiles);




    vector<TH1F*> h1D_njets_vec;
    TH1F* h1D_njets_reco = new TH1F("", ";;Entries", 15, 0, 15); h1D_njets_vec.push_back(h1D_njets_reco);
    TH1F* h1D_njets_gen = new TH1F("", ";;Entries", 15, 0, 15); h1D_njets_vec.push_back(h1D_njets_gen);

    vector<TH1F*> h1D_mindR_vec;
    TH1F* h1D_mindR = new TH1F("", "min #DeltaR between reco jet and quark;;Entries", 100, 0.0, 2.0); h1D_mindR_vec.push_back(h1D_mindR);

    // all
    vector<TH1F*> h1D_dRqq_vec;
    TH1F* h1D_dRqq1j = new TH1F("", "#DeltaR between qq if only 1 matched jet;;Entries", 100, 0.0, 5.0); h1D_dRqq_vec.push_back(h1D_dRqq1j);
    TH1F* h1D_dRqq2j = new TH1F("", "#DeltaR between qq if 2 matched jets;;Entries", 100, 0.0, 5.0); h1D_dRqq_vec.push_back(h1D_dRqq2j);

    // 200
    vector<TH1F*> h1D_dRqq200_vec;
    TH1F* h1D_dRqq1j200 = new TH1F("", "#DeltaR between qq if only 1 matched jet (stop 200);;Entries", 100, 0.0, 5.0); h1D_dRqq200_vec.push_back(h1D_dRqq1j200);
    TH1F* h1D_dRqq2j200 = new TH1F("", "#DeltaR between qq if 2 matched jets (stop 200);;Entries", 100, 0.0, 5.0); h1D_dRqq200_vec.push_back(h1D_dRqq2j200);

    // 500
    vector<TH1F*> h1D_dRqq500_vec;
    TH1F* h1D_dRqq1j500 = new TH1F("", "#DeltaR between qq if only 1 matched jet (stop 500);;Entries", 100, 0.0, 5.0); h1D_dRqq500_vec.push_back(h1D_dRqq1j500);
    TH1F* h1D_dRqq2j500 = new TH1F("", "#DeltaR between qq if 2 matched jets (stop 500);;Entries", 100, 0.0, 5.0); h1D_dRqq500_vec.push_back(h1D_dRqq2j500);

    // 900
    vector<TH1F*> h1D_dRqq900_vec;
    TH1F* h1D_dRqq1j900 = new TH1F("", "#DeltaR between qq if only 1 matched jet (stop 900);;Entries", 100, 0.0, 5.0); h1D_dRqq900_vec.push_back(h1D_dRqq1j900);
    TH1F* h1D_dRqq2j900 = new TH1F("", "#DeltaR between qq if 2 matched jets (stop 900);;Entries", 100, 0.0, 5.0); h1D_dRqq900_vec.push_back(h1D_dRqq2j900);

    vector<TH1F*> h1D_stop_mass_vec;
    TH1F* h1D_stop_mass = new TH1F("", "Stop mass (events have 2 matched jets);;Entries", 100, 100, 900); h1D_stop_mass_vec.push_back(h1D_stop_mass); 

    vector<TH1F*> h1D_stop_mass_lt2match_vec;
    TH1F* h1D_stop_mass_lt2match = new TH1F("", "Stop mass (events have <2 matched jets);;Entries", 100, 100, 900); h1D_stop_mass_lt2match_vec.push_back(h1D_stop_mass_lt2match); 

    vector<TH1F*> h1D_stop_mass_all_vec;
    TH1F* h1D_stop_mass_all = new TH1F("", "Stop mass (events with at least 1 W->qq);;Entries", 100, 100, 900); h1D_stop_mass_all_vec.push_back(h1D_stop_mass_all); 

    vector<TH1F*> h1D_lonelyjet_vec;
    TH1F* h1D_lonelyjet = new TH1F("", "Just 1 jet;;Entries", 100, 100, 900); h1D_lonelyjet_vec.push_back(h1D_lonelyjet); 

    vector<TH1F*> h1D_nonlonelyjet_vec;
    TH1F* h1D_nonlonelyjet = new TH1F("", "2 jets;;Entries", 100, 100, 900); h1D_nonlonelyjet_vec.push_back(h1D_nonlonelyjet); 

    vector<TH1F*> h1D_dRjj_masses_vec;
    TH1F* h1D_dRjj200 = new TH1F("", ";;Entries", 100, 0.0, 5.0); h1D_dRjj_masses_vec.push_back(h1D_dRjj200);
    TH1F* h1D_dRjj300 = new TH1F("", ";;Entries", 100, 0.0, 5.0); h1D_dRjj_masses_vec.push_back(h1D_dRjj300);
    TH1F* h1D_dRjj400 = new TH1F("", ";;Entries", 100, 0.0, 5.0); h1D_dRjj_masses_vec.push_back(h1D_dRjj400);
    TH1F* h1D_dRjj500 = new TH1F("", ";;Entries", 100, 0.0, 5.0); h1D_dRjj_masses_vec.push_back(h1D_dRjj500);
    TH1F* h1D_dRjj600 = new TH1F("", ";;Entries", 100, 0.0, 5.0); h1D_dRjj_masses_vec.push_back(h1D_dRjj600);
    TH1F* h1D_dRjj700 = new TH1F("", ";;Entries", 100, 0.0, 5.0); h1D_dRjj_masses_vec.push_back(h1D_dRjj700);
    TH1F* h1D_dRjj800 = new TH1F("", ";;Entries", 100, 0.0, 5.0); h1D_dRjj_masses_vec.push_back(h1D_dRjj800);
    TH1F* h1D_dRjj900 = new TH1F("", ";;Entries", 100, 0.0, 5.0); h1D_dRjj_masses_vec.push_back(h1D_dRjj900);

    vector<TH1F*> h1D_dRqq_nomatch_vec;
    TH1F* h1D_dRqq200_nomatch = new TH1F("", ";;Entries", 100, 0.0, 5.0); h1D_dRqq_nomatch_vec.push_back(h1D_dRqq200_nomatch);
    TH1F* h1D_dRqq300_nomatch = new TH1F("", ";;Entries", 100, 0.0, 5.0); h1D_dRqq_nomatch_vec.push_back(h1D_dRqq300_nomatch);
    TH1F* h1D_dRqq400_nomatch = new TH1F("", ";;Entries", 100, 0.0, 5.0); h1D_dRqq_nomatch_vec.push_back(h1D_dRqq400_nomatch);
    TH1F* h1D_dRqq500_nomatch = new TH1F("", ";;Entries", 100, 0.0, 5.0); h1D_dRqq_nomatch_vec.push_back(h1D_dRqq500_nomatch);
    TH1F* h1D_dRqq600_nomatch = new TH1F("", ";;Entries", 100, 0.0, 5.0); h1D_dRqq_nomatch_vec.push_back(h1D_dRqq600_nomatch);
    TH1F* h1D_dRqq700_nomatch = new TH1F("", ";;Entries", 100, 0.0, 5.0); h1D_dRqq_nomatch_vec.push_back(h1D_dRqq700_nomatch);
    TH1F* h1D_dRqq800_nomatch = new TH1F("", ";;Entries", 100, 0.0, 5.0); h1D_dRqq_nomatch_vec.push_back(h1D_dRqq800_nomatch);
    TH1F* h1D_dRqq900_nomatch = new TH1F("", ";;Entries", 100, 0.0, 5.0); h1D_dRqq_nomatch_vec.push_back(h1D_dRqq900_nomatch);

    vector<TH1F*> h1D_leadjmass_masses_vec;
    TH1F* h1D_leadjmass200 = new TH1F("", ";;Entries", 30, 0.0, 100.0); h1D_leadjmass_masses_vec.push_back(h1D_leadjmass200);
    TH1F* h1D_leadjmass300 = new TH1F("", ";;Entries", 30, 0.0, 100.0); h1D_leadjmass_masses_vec.push_back(h1D_leadjmass300);
    TH1F* h1D_leadjmass400 = new TH1F("", ";;Entries", 30, 0.0, 100.0); h1D_leadjmass_masses_vec.push_back(h1D_leadjmass400);
    TH1F* h1D_leadjmass500 = new TH1F("", ";;Entries", 30, 0.0, 100.0); h1D_leadjmass_masses_vec.push_back(h1D_leadjmass500);
    TH1F* h1D_leadjmass600 = new TH1F("", ";;Entries", 30, 0.0, 100.0); h1D_leadjmass_masses_vec.push_back(h1D_leadjmass600);
    TH1F* h1D_leadjmass700 = new TH1F("", ";;Entries", 30, 0.0, 100.0); h1D_leadjmass_masses_vec.push_back(h1D_leadjmass700);
    TH1F* h1D_leadjmass800 = new TH1F("", ";;Entries", 30, 0.0, 100.0); h1D_leadjmass_masses_vec.push_back(h1D_leadjmass800);
    TH1F* h1D_leadjmass900 = new TH1F("", ";;Entries", 30, 0.0, 100.0); h1D_leadjmass_masses_vec.push_back(h1D_leadjmass900);





    initCounter();
    // File Loop
    while ( (currentFile = (TFile*)fileIter.Next()) ) { 

        // Get File Content
        TFile *file = new TFile( currentFile->GetTitle() );
        TTree *tree = (TTree*)file->Get("Events");
        cms2.Init(tree);

        // Loop over Events in current file
        unsigned int nEventsTree = tree->GetEntriesFast();
        for( unsigned int event = 0; event < nEventsTree; ++event) {


            // if(nEventsSoFar > 300) break;
            // if(event > 300) break;

            addToCounter("total");

            // Get Event Content
            cms2.GetEntry(event);
            nEventsSoFar++;

            // each element of this vector is a pair of genps_ indices for the 2 quarks from a W
            vector<pair<int, int> > WtoqqIdxs; 
            // each element of this vector is a pair of goodJet indices for the 2 matched jets for the quarks
            vector<pair<int, int> > qqMatchedJets; 

            int stopMass = (int)sparm_values().at(0);

            vector<int> WpDaughterIdxs, WmDaughterIdxs;
            for(unsigned int iPart = 0; iPart < genps_id().size(); iPart++) {

                // same requirements on quarks as jets
                if (genps_p4().at(iPart).pt() < jetPtCut) continue;
                if (fabs(genps_p4().at(iPart).eta()) > 2.4) continue;

                if(genps_id_mother().at(iPart) == 24) WpDaughterIdxs.push_back(iPart);
                else if(genps_id_mother().at(iPart) == -24) WmDaughterIdxs.push_back(iPart);
                else continue;
            }

            if(WpDaughterIdxs.size() == 2) {
                if(genps_id().at( WpDaughterIdxs[0] ) < 7 && genps_id().at( WpDaughterIdxs[1] ) < 7) 
                    WtoqqIdxs.push_back( std::make_pair(WpDaughterIdxs[0], WpDaughterIdxs[1]) );
            }
            if(WmDaughterIdxs.size() == 2) {
                if(genps_id().at( WmDaughterIdxs[0] ) < 7 && genps_id().at( WmDaughterIdxs[1] ) < 7) 
                    WtoqqIdxs.push_back( std::make_pair(WmDaughterIdxs[0], WmDaughterIdxs[1]) );
            }

            if(WtoqqIdxs.size() < 1) continue; // we don't even have 1 W->qq :(
            addToCounter("at least 1 Wtoqq");

            fill(h1D_stop_mass_all, stopMass);

            std::vector<JetStruct> goodJets;
            for (unsigned int iJet = 0; iJet < pfjets_p4().size(); iJet++){
                if (pfjets_p4().at(iJet).pt()*pfjets_corL1FastL2L3().at(iJet) < jetPtCut) continue;
                if (fabs(pfjets_p4().at(iJet).eta()) > 2.4) continue;
                // if (!passesLoosePFJetID(iJet)) continue;


                ///// FIXME passesLoosePFJetID implementation
                float pfjet_chf_  = cms2.pfjets_chargedHadronE()[iJet] / cms2.pfjets_p4()[iJet].energy();
                float pfjet_nhf_  = cms2.pfjets_neutralHadronE()[iJet] / cms2.pfjets_p4()[iJet].energy();
                float pfjet_cef_  = cms2.pfjets_chargedEmE()[iJet] / cms2.pfjets_p4()[iJet].energy();
                float pfjet_nef_  = cms2.pfjets_neutralEmE()[iJet] / cms2.pfjets_p4()[iJet].energy();
                int   pfjet_cm_   = cms2.pfjets_chargedMultiplicity()[iJet];
                if (cms2.pfjets_pfcandIndicies().size() < 2) continue;
                if (pfjet_nef_ >= 0.99) continue;
                if (pfjet_nhf_ >= 0.99) continue;
                if (pfjet_cm_ < 1) continue;
                if (pfjet_chf_ < 1e-6) continue;
                if (pfjet_cef_ >= 0.99) continue;
                //// end passesLoosePFJetID implementation


                JetStruct myJet = {*(new LorentzVector()), 0.0, 999};
                myJet.jet = pfjets_p4().at(iJet);
                myJet.mindRToQuark = 9999;
                myJet.idx = iJet;

                goodJets.push_back(myJet);
            }


            for(unsigned int iPair = 0; iPair < WtoqqIdxs.size(); iPair++) {

                int q1idx = WtoqqIdxs[iPair].first;
                int q2idx = WtoqqIdxs[iPair].second;

                float mindR1 = 9999, mindR2 = 9999;
                int j1idx = -1, j2idx = -1;
                for(unsigned int iJet = 0; iJet < goodJets.size(); iJet++) {
                    float dR1 = deltaR( goodJets[iJet].jet, genps_p4().at(q1idx) );
                    float dR2 = deltaR( goodJets[iJet].jet, genps_p4().at(q2idx) );
                    if(dR1 < mindR1) {
                        mindR1 = dR1;
                        j1idx = iJet;
                    }
                    if(dR2 < mindR2) {
                        mindR2 = dR2;
                        j2idx = iJet;
                    }
                }
                float dRqq = deltaR( genps_p4().at(q1idx), genps_p4().at(q2idx) );

                if(        stopMass == 200 ) {
                    fill(h1D_dRqq200_nomatch, dRqq);
                } else if( stopMass == 300 ) {
                    fill(h1D_dRqq300_nomatch, dRqq);
                } else if( stopMass == 400 ) {
                    fill(h1D_dRqq400_nomatch, dRqq);
                } else if( stopMass == 500 ) {
                    fill(h1D_dRqq500_nomatch, dRqq);
                } else if( stopMass == 600 ) {
                    fill(h1D_dRqq600_nomatch, dRqq);
                } else if( stopMass == 700 ) {
                    fill(h1D_dRqq700_nomatch, dRqq);
                } else if( stopMass == 800 ) {
                    fill(h1D_dRqq800_nomatch, dRqq);
                } else if( stopMass == 900 ) {
                    fill(h1D_dRqq900_nomatch, dRqq);
                }

                if(j1idx > -1 && j2idx > -1) {

                    if(j1idx != j2idx) {
                        goodJets[j1idx].mindRToQuark = mindR1;
                        goodJets[j2idx].mindRToQuark = mindR2;

                        qqMatchedJets.push_back( std::make_pair( j1idx, j2idx ) );

                        fill(h1D_nonlonelyjet, stopMass);

                        fill(h1D_dRqq2j, dRqq);
                        if( stopMass == 200 ) fill(h1D_dRqq2j200, dRqq);
                        if( stopMass == 500 ) fill(h1D_dRqq2j500, dRqq);
                        if( stopMass == 900 ) fill(h1D_dRqq2j900, dRqq);
                    } else {
                        // cout << "just 1 jet that matches both quarks :(" << endl;
                        fill(h1D_lonelyjet, stopMass);

                        fill(h1D_dRqq1j, dRqq);
                        if( stopMass == 200 ) fill(h1D_dRqq1j200, dRqq);
                        if( stopMass == 500 ) fill(h1D_dRqq1j500, dRqq);
                        if( stopMass == 900 ) fill(h1D_dRqq1j900, dRqq);
                    }
                } else {
                    // cout << "no good jets?!" << goodJets.size() << " " << pfjets_p4().size() <<  endl;
                }

                fill(h1D_mindR, mindR1);
                fill(h1D_mindR, mindR2);

            }

            for(unsigned int iPair = 0; iPair < qqMatchedJets.size(); iPair++) {
                int j1idx = qqMatchedJets[iPair].first;
                int j2idx = qqMatchedJets[iPair].second;

                int nMatchedJets = 0;
                // if(goodJets[j1idx].mindRToQuark > 0.1) continue;
                // if(goodJets[j2idx].mindRToQuark > 0.1) continue;
                if(goodJets[j1idx].mindRToQuark < 0.1) nMatchedJets++;
                if(goodJets[j2idx].mindRToQuark < 0.1) nMatchedJets++;

                if(nMatchedJets < 2) {
                    fill(h1D_stop_mass_lt2match, stopMass);
                    continue;
                }

                // at this point, both quarks for this pair have 2 matched jets

                float dRjj = deltaR( goodJets[j1idx].jet, goodJets[j2idx].jet );

                float leadingJetInvMass = goodJets[j1idx].jet.pt() > goodJets[j2idx].jet.pt() ? 
                    goodJets[j1idx].jet.mass() : goodJets[j2idx].jet.mass();

                if(        stopMass == 200 ) {
                    fill(h1D_dRjj200, dRjj);
                    fill(h1D_leadjmass200, leadingJetInvMass);
                } else if( stopMass == 300 ) {
                    fill(h1D_dRjj300, dRjj);
                    fill(h1D_leadjmass300, leadingJetInvMass);
                } else if( stopMass == 400 ) {
                    fill(h1D_dRjj400, dRjj);
                    fill(h1D_leadjmass400, leadingJetInvMass);
                } else if( stopMass == 500 ) {
                    fill(h1D_dRjj500, dRjj);
                    fill(h1D_leadjmass500, leadingJetInvMass);
                } else if( stopMass == 600 ) {
                    fill(h1D_dRjj600, dRjj);
                    fill(h1D_leadjmass600, leadingJetInvMass);
                } else if( stopMass == 700 ) {
                    fill(h1D_dRjj700, dRjj);
                    fill(h1D_leadjmass700, leadingJetInvMass);
                } else if( stopMass == 800 ) {
                    fill(h1D_dRjj800, dRjj);
                    fill(h1D_leadjmass800, leadingJetInvMass);
                } else if( stopMass == 900 ) {
                    fill(h1D_dRjj900, dRjj);
                    fill(h1D_leadjmass900, leadingJetInvMass);
                }

                fill(h1D_stop_mass, stopMass);
            } // outside of (beyond) this loop, only condition is that we have at least 1 W->qq

            CMS2::progress( nEventsSoFar, nEventsTotal );

            fill(h1D_njets_gen, genjets_p4().size());
            fill(h1D_njets_reco, goodJets.size());


            addToCounter("good events");

            nGoodEvents++;

        }//event loop

    }//file loop MCMC

    printCounter();

    std::cout << " nGoodEvents: " << nGoodEvents << " nEventsTotal: " << nEventsTotal << std::endl;


    TString prefix("plots");
    prefix += tag;
    prefix += "/";

    TString common = "";
    TString massTitles = " --titles m_{#tilde{t}}=200| \
                          m_{#tilde{t}}=300| \
                          m_{#tilde{t}}=400| \
                          m_{#tilde{t}}=500| \
                          m_{#tilde{t}}=600| \
                          m_{#tilde{t}}=700| \
                          m_{#tilde{t}}=800| \
                          m_{#tilde{t}}=900| \
                          ";

    TH1F* dog = new TH1F("dog","",1,0,1);



    drawStacked(dog, h1D_njets_vec,prefix+"h1D_njets.pdf","--nostack --title Njets --titles reco|gen"+common);
    drawStacked(dog, h1D_stop_mass_vec,prefix+"h1D_stop_mass.pdf",""+common);
    drawStacked(dog, h1D_stop_mass_lt2match_vec,prefix+"h1D_stop_mass_lt2match.pdf",""+common);
    drawStacked(dog, h1D_stop_mass_all_vec,prefix+"h1D_stop_mass_all.pdf",""+common);
    drawStacked(dog, h1D_lonelyjet_vec,prefix+"h1D_lonelyjet.pdf",""+common);

    h1D_lonelyjet->Divide(h1D_nonlonelyjet);
    drawStacked(dog, h1D_lonelyjet_vec,prefix+"h1D_lonelyjetfraction.pdf","--title lonely jet ratio (1 jet/2 jets)"+common);

    h1D_stop_mass_lt2match->Divide(h1D_stop_mass);
    drawStacked(dog, h1D_stop_mass_lt2match_vec,prefix+"h1D_jetmatchfraction.pdf","--title ratio of evts with <2 matched jets to 2 matched jets"+common);

    drawStacked(dog, h1D_dRqq_vec,prefix+"h1D_dRqq.pdf","--centerlabel --nostack --title dR between qq --titles if 1 common mindR jet|if 2 distinct jets"+common);
    drawStacked(dog, h1D_dRqq200_vec,prefix+"h1D_dRqq200.pdf","--centerlabel --nostack --title dR between qq (stop 200) --titles if 1 common mindR jet|if 2 distinct jets"+common);
    drawStacked(dog, h1D_dRqq500_vec,prefix+"h1D_dRqq500.pdf","--centerlabel --nostack --title dR between qq (stop 500) --titles if 1 common mindR jet|if 2 distinct jets"+common);
    drawStacked(dog, h1D_dRqq900_vec,prefix+"h1D_dRqq900.pdf","--centerlabel --nostack --title dR between qq (stop 900) --titles if 1 common mindR jet|if 2 distinct jets"+common);

    drawStacked(dog, h1D_dRjj_masses_vec,prefix+"h1D_dRjj_masses.pdf","--keeporder --centerlabel --nostack --title dR between jj matched to qq"+massTitles+common);
    drawStacked(dog, h1D_dRqq_nomatch_vec,prefix+"h1D_dRqq_nomatch.pdf","--keeporder --centerlabel --nostack --title dR between qq (no matching at this point)"+massTitles+common);
    drawStacked(dog, h1D_leadjmass_masses_vec,prefix+"h1D_leadjmass_masses.pdf","--keeporder --centerlabel --nostack --title inv mass of leading j (of jj matched to qq) --nofill --normalize "+massTitles+common);

    drawStacked(dog, h1D_mindR_vec,prefix+"h1D_mindR.pdf","--centerlabel --title min #DeltaR between q-j "+common);

    vector<vector<float> > xvals;
    vector<vector<float> > yvals_jj, yvals_qq;

    vector<float> thresholds;
    TString kTitles = " --titles ";
    thresholds.push_back(0.6); kTitles += "k = 0.6|";
    thresholds.push_back(0.7); kTitles += "k = 0.7|";
    thresholds.push_back(0.8); kTitles += "k = 0.8|";
    thresholds.push_back(0.9); kTitles += "k = 0.9|";
    thresholds.push_back(1.0); kTitles += "k = 1.0|";

    for(unsigned int i = 0; i < thresholds.size(); i++) {
        vector<float> masses, fatfractions, qqfractions;
        masses.push_back(200); 
        masses.push_back(300); 
        masses.push_back(400); 
        masses.push_back(500); 
        masses.push_back(600); 
        masses.push_back(700); 
        masses.push_back(800); 
        masses.push_back(900); 

        fatfractions.push_back( getFractionBetween(h1D_dRjj200, 0.0, thresholds.at(i)) );
        fatfractions.push_back( getFractionBetween(h1D_dRjj300, 0.0, thresholds.at(i)) );
        fatfractions.push_back( getFractionBetween(h1D_dRjj400, 0.0, thresholds.at(i)) );
        fatfractions.push_back( getFractionBetween(h1D_dRjj500, 0.0, thresholds.at(i)) );
        fatfractions.push_back( getFractionBetween(h1D_dRjj600, 0.0, thresholds.at(i)) );
        fatfractions.push_back( getFractionBetween(h1D_dRjj700, 0.0, thresholds.at(i)) );
        fatfractions.push_back( getFractionBetween(h1D_dRjj800, 0.0, thresholds.at(i)) );
        fatfractions.push_back( getFractionBetween(h1D_dRjj900, 0.0, thresholds.at(i)) );

        qqfractions.push_back( getFractionBetween(h1D_dRqq200_nomatch, 0.0, thresholds.at(i)) );
        qqfractions.push_back( getFractionBetween(h1D_dRqq300_nomatch, 0.0, thresholds.at(i)) );
        qqfractions.push_back( getFractionBetween(h1D_dRqq400_nomatch, 0.0, thresholds.at(i)) );
        qqfractions.push_back( getFractionBetween(h1D_dRqq500_nomatch, 0.0, thresholds.at(i)) );
        qqfractions.push_back( getFractionBetween(h1D_dRqq600_nomatch, 0.0, thresholds.at(i)) );
        qqfractions.push_back( getFractionBetween(h1D_dRqq700_nomatch, 0.0, thresholds.at(i)) );
        qqfractions.push_back( getFractionBetween(h1D_dRqq800_nomatch, 0.0, thresholds.at(i)) );
        qqfractions.push_back( getFractionBetween(h1D_dRqq900_nomatch, 0.0, thresholds.at(i)) );

        xvals.push_back(masses);
        yvals_jj.push_back(fatfractions);
        yvals_qq.push_back(qqfractions);
    }


    drawGraph(xvals, yvals_jj, prefix+"g1D_fatfraction.pdf", "--title Fraction of jj with #DeltaR_{jj}<k --xlabel m_{#tilde{t}}-m_{LSP} --ylabel fat fraction --legendposition bottom"+kTitles);
    drawGraph(xvals, yvals_qq, prefix+"g1D_qqfraction.pdf", "--title Fraction of qq (no match) with #DeltaR_{jj}<k --xlabel m_{#tilde{t}}-m_{LSP} --ylabel fraction --legendposition bottom"+kTitles);


    // return 0;
}
