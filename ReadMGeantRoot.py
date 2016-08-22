'''This is just an example script on how to read the MGeant 
files from Alex.  Not really useful for much else.'''

from ROOT import TFile, AddressOf
from Sim2Root import MGSimulation

MGS = MGSimulation()

train = TFile('100MeV_8cm_0deg.root','read')
t = train.Get("h10")

t.SetBranchAddress("Runevt", AddressOf(MGS.MGStruct, 'Runevt'))
t.SetBranchAddress("Xcor", AddressOf(MGS.MGStruct, 'Xcor'))
t.SetBranchAddress("Ycor", AddressOf(MGS.MGStruct, 'Ycor'))
t.SetBranchAddress("Zcor", AddressOf(MGS.MGStruct, 'Zcor'))
t.SetBranchAddress("Estep", AddressOf(MGS.MGStruct, 'Estep'))
t.SetBranchAddress("Itra", AddressOf(MGS.MGStruct, 'Itra'))

a = t.GetEntries()

for i in range(0,a)[1000:1100]:
    t.GetEntry(i)
    print MGS.MGStruct.Runevt, MGS.MGStruct.Xcor
