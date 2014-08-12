#include "TString.h"
#include "TH1.h"
#include "THStack.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TString.h"
#include "TMath.h"
#include "TRegexp.h"
#include "TPRegexp.h"
#include "TPavesText.h"

#include <stdlib.h>

using namespace std;
TString pdgToStr(int id) {
    TString out = "";
    bool isSusy = false;
    if(abs(id) > 1000000) {
        isSusy = true;
        if(id > 0) id -= 1000000;
        else id += 1000000;
        // id = ( (id>0)-(id<0) ) * (id - 1000000);
    }
    switch(abs(id)) {
        case 1: out += (id > 0 ? "d" : "dbar"); break;
        case 2: out += (id > 0 ? "u" : "ubar"); break;
        case 3: out += (id > 0 ? "s" : "sbar"); break;
        case 4: out += (id > 0 ? "c" : "cbar"); break;
        case 5: out += (id > 0 ? "b" : "bbar"); break;
        case 6: out += (id > 0 ? "t" : "tbar"); break;

        case 11: out += (id > 0 ? "e-" : "e+"); break;
        case 12: out += (id > 0 ? "ve" : "vebar"); break;
        case 13: out += (id > 0 ? "mu-" : "mu+"); break;
        case 14: out += (id > 0 ? "vmu" : "vmubar"); break;
        case 15: out += (id > 0 ? "tau-" : "tau+"); break;
        case 16: out += (id > 0 ? "vtau" : "vtaubar"); break;

        case 21: out += (id > 0 ? "g" : "g"); break;
        case 22:
                 if(isSusy) {
                     out += (id > 0 ? "LSP" : "LSP"); break;
                 } else {
                     out += (id > 0 ? "gamma" : "gamma"); break;
                 }
                 
        case 23: out += (id > 0 ? "Z0" : "Z0"); break;
        case 24: out += (id > 0 ? "W+" : "W-"); break;

        case 2212: out += (id > 0 ? "p" : "p"); break;

    }

    if(out.Length() == 0) { // we didn't find it
        // couldn't find it
        return Form("%i", id);
    } else {
        if(isSusy) out += "~";
        return out;
    }

}
