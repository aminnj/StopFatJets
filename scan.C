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
#include "./CORE/mcSelections.cc"
// #include "./CORE/jetSelections.h"
#include "/home/users/namin/macros/utils.C"
#include "./particles.C"

#define NSTOPMASSES 7

using namespace std;
using namespace tas;


struct JetStruct {
    LorentzVector jet;
    float mindRToQuark;
    unsigned int idx;
};

int stopMassToIndex(int stopMass) {
    if(      stopMass == 200 )  return 0;
    else if( stopMass == 300 )  return 1;
    else if( stopMass == 400 )  return 2;
    else if( stopMass == 500 )  return 3;
    else if( stopMass == 600 )  return 4;
    else if( stopMass == 700 )  return 5;
    else if( stopMass == 800 )  return 6;
    else if( stopMass == 900 )  return 7;

    cout << "ERROR in stopMassToIndex!!!!" << endl;
    return 0;
}

// void scan(float jetPtCut=40.0, TString tag = "", bool requireLeptonic = false, float mtLow = -1.0, float mtHigh = -1.0) {
void scan(float jetPtCut=40.0, TString tag = "", bool requireLeptonic = false) {

    TChain *ch = new TChain("Events");

    // ch->Add("/hadoop/cms/store/user/namin/mStop500to800_mLSP1/ntuple_78_2_RJ3.root");

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

    vector<TH1F*> h1D_mindR_vec;
    TH1F* h1D_mindR = new TH1F("", "min #DeltaR between reco jet and quark;;Entries", 100, 0.0, 2.0); h1D_mindR_vec.push_back(h1D_mindR);

    vector<TH1F*> h1D_stop_mass_vec;
    TH1F* h1D_stop_mass = new TH1F("", "Stop mass (events have 2 matched jets);;Entries", 100, 100, 900); h1D_stop_mass_vec.push_back(h1D_stop_mass); 

    vector<TH1F*> h1D_stop_mass_lt2match_vec;
    TH1F* h1D_stop_mass_lt2match = new TH1F("", "Stop mass (events have <2 matched jets);;Entries", 100, 100, 900); h1D_stop_mass_lt2match_vec.push_back(h1D_stop_mass_lt2match); 

    vector<TH1F*> h1D_stop_mass_all_vec;
    TH1F* h1D_stop_mass_all = new TH1F("", "Stop mass (events passing W daughter reqs);;Entries", 100, 100, 900); h1D_stop_mass_all_vec.push_back(h1D_stop_mass_all); 

    vector<TH1F*> h1D_lonelyjet_vec;
    TH1F* h1D_lonelyjet = new TH1F("", "Just 1 jet;;Entries", 100, 100, 900); h1D_lonelyjet_vec.push_back(h1D_lonelyjet); 

    vector<TH1F*> h1D_nonlonelyjet_vec;
    TH1F* h1D_nonlonelyjet = new TH1F("", "2 jets;;Entries", 100, 100, 900); h1D_nonlonelyjet_vec.push_back(h1D_nonlonelyjet); 

    TH2F* h2D_dRqq_Wpt = new TH2F("", "#DeltaR_{qq} (no matching);W_{p_{T}};#DeltaR_{qq}", 50,0.0,900, 50,0.0,5.0);

    vector<TH1F*> h1D_njets_gen_vec;
    vector<TH1F*> h1D_njets_reco_vec;
    vector<TH1F*> h1D_dRjj_masses_vec;
    vector<TH1F*> h1D_dRqq_nomatch_vec;
    vector<TH1F*> h1D_dRcone_nomatch_vec;
    vector<TH1F*> h1D_leadjmass_masses_vec;
    vector<TH1F*> h1D_mt_masses_vec;
    for(int iStop = 0; iStop <= NSTOPMASSES; iStop++) {
        TH1F* h1D_njets_gen = new TH1F("", ";;Entries", 15, 0,15); h1D_njets_gen_vec.push_back(h1D_njets_gen);
        TH1F* h1D_njets_reco = new TH1F("", ";;Entries", 15, 0,15); h1D_njets_reco_vec.push_back(h1D_njets_reco);
        TH1F* h1D_dRjj_masses = new TH1F("", ";;Entries", 100, 0.0, 5.0); h1D_dRjj_masses_vec.push_back(h1D_dRjj_masses);
        TH1F* h1D_dRqq_nomatch = new TH1F("", ";;Entries", 100, 0.0, 5.0); h1D_dRqq_nomatch_vec.push_back(h1D_dRqq_nomatch);
        TH1F* h1D_dRcone_nomatch = new TH1F("", ";;Entries", 100, 0.0, 5.0); h1D_dRcone_nomatch_vec.push_back(h1D_dRcone_nomatch);
        TH1F* h1D_leadjmass_masses = new TH1F("", ";;Entries", 30, 0.0, 100.0); h1D_leadjmass_masses_vec.push_back(h1D_leadjmass_masses);
        TH1F* h1D_mt_masses = new TH1F("", ";;Entries", 45, 0.0, 900.0); h1D_mt_masses_vec.push_back(h1D_mt_masses);
    }

    initCounter();
    // File Loop
    while ( (currentFile = (TFile*)fileIter.Next()) ) { 

        // Get File Content
        TFile *file = new TFile( currentFile->GetTitle() );
        TTree *tree = (TTree*)file->Get("Events");
        cms2.Init(tree);

        TString filename(currentFile->GetTitle());

        // Loop over Events in current file
        unsigned int nEventsTree = tree->GetEntriesFast();
        for( unsigned int event = 0; event < nEventsTree; ++event) {


            // if(event > 300) break;

            addToCounter("total");

            // Get Event Content
            cms2.GetEntry(event);
            nEventsSoFar++;

            // each element of this vector is a pair of genps_ indices for the 2 quarks from a W
            vector<pair<int, int> > WtoqqIdxs; 
            // each element of this vector is a pair of goodJet indices for the 2 matched jets for the quarks
            vector<pair<int, int> > qqMatchedJets; 
            // indices of b quarks from same tops as Ws; this is in same order as WtoqqIdxs
            vector<int> bIdxs; 
            // indices of Ws
            vector<int> WIdxs;

            int stopMass = (int)sparm_values().at(0);
            int stopMassIdx = stopMassToIndex(stopMass);

            vector<int> WpDaughterIdxs, WmDaughterIdxs;
            for(unsigned int iPart = 0; iPart < genps_id().size(); iPart++) {
                //store Ws
                // if(genps_id().at(iPart) == 24) WIdxs[0] = iPart;
                // else if(genps_id().at(iPart) == -24) WIdxs[1] = iPart;

                if(genps_status().at(iPart) != 3) continue; // only look at status=3

                // same requirements on quarks as jets
                if (genps_p4().at(iPart).pt() < jetPtCut) continue;
                if (fabs(genps_p4().at(iPart).eta()) > 2.5) continue;

                if(genps_id_mother().at(iPart) == 24) {
                    WpDaughterIdxs.push_back(iPart);
                }
                else if(genps_id_mother().at(iPart) == -24) {
                    WmDaughterIdxs.push_back(iPart);
                }
                else continue;
            }


            int nLeptonicDecay = 0;
            int leptonIdx = -1;
            if(WpDaughterIdxs.size() == 2) {
                if(abs(genps_id().at( WpDaughterIdxs[0] )) < 7 && abs(genps_id().at( WpDaughterIdxs[1] )) < 7) 
                    WtoqqIdxs.push_back( std::make_pair(WpDaughterIdxs[0], WpDaughterIdxs[1]) );
                else if(genps_id().at( WpDaughterIdxs[0] ) * genps_id().at( WpDaughterIdxs[1] ) < 0) {  // opposite signs for pdg ids
                    if( abs( abs(genps_id().at( WpDaughterIdxs[0] )) - abs(genps_id().at( WpDaughterIdxs[1] )) ) == 1 ) {  // lep + nu are same flavor
                        nLeptonicDecay++;
                        // lepton pdg codes are odd
                        if( abs( genps_id().at( WpDaughterIdxs[0] ) ) % 2 != 0) leptonIdx = WpDaughterIdxs[0];
                        else leptonIdx = WpDaughterIdxs[1];
                    }
                }
            }
            if(WmDaughterIdxs.size() == 2) {
                if(abs(genps_id().at( WmDaughterIdxs[0] )) < 7 && abs(genps_id().at( WmDaughterIdxs[1] )) < 7) 
                    WtoqqIdxs.push_back( std::make_pair(WmDaughterIdxs[0], WmDaughterIdxs[1]) );
                else if(genps_id().at( WmDaughterIdxs[0] ) * genps_id().at( WmDaughterIdxs[1] ) < 0) {  // opposite signs for pdg ids
                    if( abs( abs(genps_id().at( WmDaughterIdxs[0] )) - abs(genps_id().at( WmDaughterIdxs[1] )) ) == 1 ) {  // lep + nu are same flavor
                        nLeptonicDecay++;
                        // lepton pdg codes are odd
                        if( abs( genps_id().at( WmDaughterIdxs[0] ) ) % 2 != 0) leptonIdx = WmDaughterIdxs[0];
                        else leptonIdx = WmDaughterIdxs[1];
                    }
                }
            }

            if(requireLeptonic) {
                // want 1 W to go hadronically and 1 leptonically
                if(WtoqqIdxs.size() != 1) continue;
                if(nLeptonicDecay != 1) continue;

                if(leptonIdx == -1) cout << "woah dude, error, leptonIdx is -1" << endl;

            } else {

                if(WtoqqIdxs.size() < 1) continue; // we don't even have 1 W->qq :(
            }

            // cout << endl;
            // // cout << "_______________" << endl;

            // find bs (and Ws) from top that Ws came from. then we have the whole familia: b,q,q
            for(unsigned int iPair = 0; iPair < WtoqqIdxs.size(); iPair++) {

                // take first quark for example and get W id 
                int Wid = genps_id_mother().at( WtoqqIdxs[iPair].first );
                // std::cout << " Wid: " << Wid << std::endl;


                // std::cout << " genps_id().size(): " << genps_id().size() << std::endl;
                
                for(unsigned int iPart = 0; iPart < genps_id().size(); iPart++) {

                    // iPart genps_id().at(iPart) Wid genps_id().at(iPart)==Wid
                    // std::cout << " iPart: " << iPart << " genps_id().at(iPart): " << genps_id().at(iPart) << " Wid: " << Wid  << std::endl;
                    // std::cout << " genps_id().at(iPart)-Wid: " << genps_id().at(iPart) - Wid << std::endl;
                    // std::cout << " iPart: " << iPart << " genps_id().at(iPart): " << genps_id().at(iPart) << " Wid: " << Wid << " genps_id().at(iPart)==Wid: " << genps_id().at(iPart)==Wid << std::endl;

                    if(genps_id().at(iPart) == Wid) WIdxs.push_back(iPart);
                    // std::cout << " iPart: " << iPart << std::endl;

                    // only consider b quarks
                    if(abs(genps_id().at(iPart)) != 5) continue;

                    int bid = genps_id().at(iPart);

                    // W+ goes with b; W- goes with bbar => product of W,b ID should be positive
                    if(Wid * bid > 0) {
                        bIdxs.push_back(iPart);
                        // break;
                    }
                }
                // if(WIdxs.size() < 1) { // XXX
                //     std::cout << " Wid: " << Wid << " WtoqqIdxs[iPair].first: " << WtoqqIdxs[iPair].first  << std::endl;
                //     for(unsigned int iPart = 0; iPart < genps_id().size(); iPart++) {
                //         std::cout << " iPart: " << iPart << " genps_id().at(iPart): " << genps_id().at(iPart) << std::endl;
                //     }
                // }
            }

            // if(WIdxs.size() < 1) { // XXX
            //     std::cout << " filename: " << filename << std::endl;
            //     dumpDocLines();
            // }
            // for(unsigned int iB = 0; iB < bIdxs.size(); iB++) {
            // std::cout << " genps_p4().at(bIxs[iB]).pt(): " << genps_p4().at(bIdxs[iB]).pt() << std::endl;
            // }
            if(bIdxs.size() != WtoqqIdxs.size() || WIdxs.size() < 1) {
                // t->qW, but q is mainly b. if it is not b, skip the event
                // also, for some reason, sometimes we can't find the W?! very rare, don't care to debug
                // dumpDocLines();
                continue;
            }

            // std::cout << " genps_p4().at(WIdxs[iPair]).pt(): " << genps_p4().at(WIdxs[0]).pt() << std::endl;


            // std::cout << " WIdxs.size(): " << WIdxs.size() << " bIdxs.size(): " << bIdxs.size() << " WtoqqIdxs.size(): " << WtoqqIdxs.size()  << std::endl;


            // mt stuff
            // sweet. gen met considers only lsps+neutrinos! 
            float mt = MT( genps_p4().at(leptonIdx), gen_met(), gen_metPhi() );
            fill(h1D_mt_masses_vec.at(stopMassIdx), mt);

            //             if(mtLow > 0 && mtHigh > 0) {
            //                 // std::cout << " mt: " << mt << " mtLow: " << mtLow << " mtHigh: " << mtHigh << std::endl;
            //                 if(mt > mtHigh || mt < mtLow) {
            //                     continue;
            //                 }
            //             }


            // end mt stuff

            fill(h1D_stop_mass_all, stopMass);

            std::vector<JetStruct> goodJets;
            std::vector<JetStruct> goodGenJets;
            for (unsigned int iJet = 0; iJet < pfjets_p4().size(); iJet++){
                if (pfjets_p4().at(iJet).pt()*pfjets_corL1FastL2L3().at(iJet) < jetPtCut) continue;
                if (fabs(pfjets_p4().at(iJet).eta()) > 2.5) continue;
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
            for (unsigned int iJet = 0; iJet < genjets_p4().size(); iJet++) {
                if (genjets_p4().at(iJet).pt() < jetPtCut) continue;
                if (fabs(genjets_p4().at(iJet).eta()) > 2.5) continue;

                JetStruct myJet = {*(new LorentzVector()), 0.0, 999};
                myJet.jet = pfjets_p4().at(iJet);
                myJet.mindRToQuark = 9999;
                myJet.idx = iJet;

                goodGenJets.push_back(myJet);
            }

            fill(h1D_njets_gen_vec.at(stopMassIdx), goodGenJets.size());
            fill(h1D_njets_reco_vec.at(stopMassIdx), goodJets.size());

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
                fill(h1D_dRqq_nomatch_vec.at(stopMassIdx), dRqq);
                fill2D(h2D_dRqq_Wpt, genps_p4().at( WIdxs[iPair] ).pt(), dRqq);
                // cout << "before stuff" << endl; 
                // std::cout << " iPair: " << iPair << " WtoqqIdxs.size(): " << WtoqqIdxs.size() << std::endl;
                float dRq1b = deltaR( genps_p4().at( q1idx ), genps_p4().at( bIdxs[iPair] ) );
                float dRq2b = deltaR( genps_p4().at( q2idx ), genps_p4().at( bIdxs[iPair] ) );
                // cout << "after stuff" << endl; 
                float dRcone = max( max( dRq1b, dRq2b ), dRqq );
                fill(h1D_dRcone_nomatch_vec.at(stopMassIdx), dRcone);

                if(j1idx > -1 && j2idx > -1) {

                    if(j1idx != j2idx) {
                        goodJets[j1idx].mindRToQuark = mindR1;
                        goodJets[j2idx].mindRToQuark = mindR2;

                        qqMatchedJets.push_back( std::make_pair( j1idx, j2idx ) );

                        fill(h1D_nonlonelyjet, stopMass);
                    } else {
                        // cout << "just 1 jet that matches both quarks :(" << endl;
                        fill(h1D_lonelyjet, stopMass);

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
                if(goodJets[j1idx].mindRToQuark < 0.2) nMatchedJets++;
                if(goodJets[j2idx].mindRToQuark < 0.2) nMatchedJets++;

                if(nMatchedJets < 2) {
                    fill(h1D_stop_mass_lt2match, stopMass);
                    continue;
                }

                // at this point, both quarks for this pair have 2 matched jets

                float dRjj = deltaR( goodJets[j1idx].jet, goodJets[j2idx].jet );

                float leadingJetInvMass = goodJets[j1idx].jet.pt() > goodJets[j2idx].jet.pt() ? 
                    goodJets[j1idx].jet.mass() : goodJets[j2idx].jet.mass();

                fill(h1D_dRjj_masses_vec.at(stopMassIdx), dRjj);
                fill(h1D_leadjmass_masses_vec.at(stopMassIdx), leadingJetInvMass);
                // fill(h1D_njets_gen_vec.at(stopMassIdx), goodGenJets.size());
                // fill(h1D_njets_reco_vec.at(stopMassIdx), goodJets.size());

                fill(h1D_stop_mass, stopMass);
            } // outside of (beyond) this loop, only condition is that we have at least 1 W->qq

            CMS2::progress( nEventsSoFar, nEventsTotal );

            addToCounter("good events");

            // dumpDocLines();

            nGoodEvents++;

        }//event loop

    }//file loop MCMC

    printCounter();

    std::cout << " nGoodEvents: " << nGoodEvents << " nEventsTotal: " << nEventsTotal << std::endl;


    TString prefix("plots");
    prefix += tag;
    prefix += "/";

    TString common = "";
    // if(mtLow > 0 && mtHigh > 0) {
    //     common += " --label ";
    //     common += mtLow;
    //     common += " < M_{T} < ";
    //     common += mtHigh;
    // }

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


    drawStacked(dog, h1D_lonelyjet_vec,prefix+"h1D_lonelyjet.pdf",""+common);
    h1D_lonelyjet->Divide(h1D_nonlonelyjet);
    drawStacked(dog, h1D_lonelyjet_vec,prefix+"h1D_lonelyjetfraction.pdf","--title lonely jet ratio (=1 jet/2 jets) (mindR, but no dR cut yet)"+common);

    drawStacked(dog, h1D_stop_mass_lt2match_vec,prefix+"h1D_stop_mass_lt2match.pdf",""+common);
    h1D_stop_mass_lt2match->Divide(h1D_stop_mass);
    drawStacked(dog, h1D_stop_mass_lt2match_vec,prefix+"h1D_jetmatchfraction.pdf","--title ratio of evts with <2 matched jets to 2 matched jets (after dR cut)"+common);

    drawStacked(dog, h1D_njets_gen_vec,prefix+"h1D_njets_gen.pdf","--keeporder --nostack --title Ngenjets --nofill --normalize --label normalized "+massTitles+common);
    drawStacked(dog, h1D_njets_reco_vec,prefix+"h1D_njets_reco.pdf","--keeporder --nostack --title Nrecojets --nofill --normalize --label normalized "+massTitles+common);
    drawStacked(dog, h1D_stop_mass_vec,prefix+"h1D_stop_mass.pdf",""+common);
    drawStacked(dog, h1D_stop_mass_all_vec,prefix+"h1D_stop_mass_all.pdf",""+common);
    drawStacked(dog, h1D_dRjj_masses_vec,prefix+"h1D_dRjj_masses.pdf","--keeporder --centerlabel --nostack --title dR between jj matched to qq --nofill "+massTitles+common);
    drawStacked(dog, h1D_dRqq_nomatch_vec,prefix+"h1D_dRqq_nomatch.pdf","--keeporder --centerlabel --nostack --title dR between qq (no matching at this point) --nofill "+massTitles+common);
    drawStacked(dog, h1D_dRcone_nomatch_vec,prefix+"h1D_dRcone_nomatch.pdf","--keeporder --centerlabel --nostack --title dR of cone with bqq (no matching at this point) --nofill "+massTitles+common);
    drawStacked(dog, h1D_leadjmass_masses_vec,prefix+"h1D_leadjmass_masses.pdf","--keeporder --centerlabel --nostack --title inv mass of leading j (of jj matched to qq) --nofill --normalize --label normalized "+massTitles+common);
    drawStacked(dog, h1D_mt_masses_vec,prefix+"h1D_mt_masses.pdf","--keeporder --centerlabel --nostack --title M_{T} --nofill --normalize --label normalized "+massTitles+common);
    drawStacked(dog, h1D_mindR_vec,prefix+"h1D_mindR.pdf","--centerlabel --title min #DeltaR between q-j --normalize --label normalized "+common);

    vector<float> thresholds;
    TString kTitles = " --titles ";
    thresholds.push_back(0.4); kTitles += "k = 0.4|";
    thresholds.push_back(0.5); kTitles += "k = 0.5|";
    thresholds.push_back(0.6); kTitles += "k = 0.6|";
    thresholds.push_back(0.7); kTitles += "k = 0.7|";
    thresholds.push_back(0.8); kTitles += "k = 0.8|";
    thresholds.push_back(0.9); kTitles += "k = 0.9|";
    thresholds.push_back(1.0); kTitles += "k = 1.0|";

    vector<float> masses;
    for(int iStop = 0; iStop <= NSTOPMASSES; iStop++) {
        masses.push_back(iStop*100+200);
    }

    vector<vector<float> > xvals;
    vector<vector<float> > yvals_jj, yvals_qq, yvals_bqq;
    for(unsigned int i = 0; i < thresholds.size(); i++) {
        vector<float> fatfractions, qqfractions, bqqfractions;

        for(int iStop = 0; iStop <= NSTOPMASSES; iStop++) {
            fatfractions.push_back( getFractionBetween(h1D_dRjj_masses_vec.at(iStop), 0.0, thresholds.at(i)) );
            qqfractions.push_back( getFractionBetween(h1D_dRqq_nomatch_vec.at(iStop), 0.0, thresholds.at(i)) );
            bqqfractions.push_back( getFractionBetween(h1D_dRcone_nomatch_vec.at(iStop), 0.0, thresholds.at(i)) );
        }

        xvals.push_back(masses);
        yvals_jj.push_back(fatfractions);
        yvals_qq.push_back(qqfractions);
        yvals_bqq.push_back(bqqfractions);
    }


    drawGraph(xvals, yvals_jj, prefix+"g1D_fatfraction.pdf", "--title Fraction of jj with #DeltaR_{jj}<k --xlabel m_{#tilde{t}}-m_{LSP} --ylabel fat fraction --legendposition bottom "+kTitles+common);
    drawGraph(xvals, yvals_qq, prefix+"g1D_qqfraction.pdf", "--title Fraction of qq (no match) with #DeltaR_{qq}<k --xlabel m_{#tilde{t}}-m_{LSP} --ylabel fraction --legendposition bottom "+kTitles+common);
    drawGraph(xvals, yvals_bqq, prefix+"g1D_bqqfraction.pdf", "--title Fraction of bqq cones (no match) with #DeltaR_{bqqcone}<k --xlabel m_{#tilde{t}}-m_{LSP} --ylabel fraction --legendposition bottom "+kTitles+common);

    vector<vector<float> > xval_masses;
    vector<vector<float> > yval_genjets, yval_recojets;
    for(int njets = 3; njets <= 4; njets++) {
        vector<float> genjetfractions, recojetfractions;

        for(int iStop = 0; iStop <= NSTOPMASSES; iStop++) {
            genjetfractions.push_back( getFractionBetween(h1D_njets_gen_vec.at(iStop), njets, 100) );
            recojetfractions.push_back( getFractionBetween(h1D_njets_reco_vec.at(iStop), njets, 100) );
        }

        yval_genjets.push_back(genjetfractions); 
        yval_recojets.push_back(recojetfractions); 
        xval_masses.push_back(masses); 
    }

    drawGraph(xval_masses, yval_genjets, prefix+"g1D_genjetfraction.pdf", "--title Fraction of ngenjets #geq k --xlabel m_{#tilde{t}}-m_{LSP} --ylabel fraction --legendposition bottom --titles k=3|k=4 "+common);
    drawGraph(xval_masses, yval_recojets, prefix+"g1D_recojetfraction.pdf", "--title Fraction of nrecojets #geq k --xlabel m_{#tilde{t}}-m_{LSP} --ylabel fraction --lerecodposition bottom --titles k=3|k=4 "+common);

    // 2D hists
    drawHist2D(h2D_dRqq_Wpt, prefix+"h2D_dRqq_Wpt.pdf");

    // return 0;
}
