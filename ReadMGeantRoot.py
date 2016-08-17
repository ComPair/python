import ROOT

a = 2.0

ROOT.gROOT.ProcessLine(
"struct h10 {\
   Float_t         Runevt;\
   Float_t       Xcor;\
   Float_t       Ycor;\
   Float_t       Zcor;\
   Float_t       Estep;\
   Float_t       Itra;\
}" )


train = ROOT.TFile('100MeV_8cm_0deg.root','read')
t = train.Get("h10")
h = ROOT.h10()

t.SetBranchAddress("Runevt", ROOT.AddressOf(h, 'Runevt'))
t.SetBranchAddress("Xcor", ROOT.AddressOf(h, 'Xcor'))
t.SetBranchAddress("Ycor", ROOT.AddressOf(h, 'Ycor'))
t.SetBranchAddress("Zcor", ROOT.AddressOf(h, 'Zcor'))
t.SetBranchAddress("Estep", ROOT.AddressOf(h, 'Estep'))
t.SetBranchAddress("Itra", ROOT.AddressOf(h, 'Itra'))

a = t.GetEntries()

for i in range(0,a)[1000:1100]:
    t.GetEntry(i)
    print h.Runevt, h.Xcor
