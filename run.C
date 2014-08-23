{

TStopwatch * timer = new TStopwatch();
timer->Start();



gROOT->ProcessLine(".L ./CORE/CMS2.cc+");
gROOT->ProcessLine(".L /home/users/namin/macros/utils.C+");
gROOT->ProcessLine(".L scan.C+");
// gROOT->ProcessLine(".L fraction.C+");

// scan(40.0, "", false);
//     jetPt   tag               lept     mtLow    mtHigh
// scan(  40.0,           "_lept",  true,    -1.0,    -1.0);
scan(  40.0,           "_lept",  true);
// scan(  40.0,   "_lept100to200",  true,    -1.0,    -1.0);

timer->Stop();

timer->Print();
}
