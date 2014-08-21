{

TStopwatch * timer = new TStopwatch();
timer->Start();



gROOT->ProcessLine(".L ./CORE/CMS2.cc+");
gROOT->ProcessLine(".L scan.C+");
// gROOT->ProcessLine(".L fraction.C+");

scan(40.0, "", false);
scan(40.0, "_lept", true);
// scan(40.0, "_pt40");
// scan(30.0, "_pt30");
// scan(20.0, "_pt20");
// scan(20.0, "_pt20");
// scan(30.0, "_pt30");

timer->Stop();

timer->Print();
}
