{

TStopwatch * timer = new TStopwatch();
timer->Start();



gROOT->ProcessLine(".L ./CORE/CMS2.cc+");
gROOT->ProcessLine(".L scan.C+");

scan();
// scan(40.0, "_pt40");
// scan(30.0, "_pt30");
// scan(20.0, "_pt20");
// scan(20.0, "_pt20");
// scan(30.0, "_pt30");

timer->Stop();

timer->Print();
}
