from EventViewer import parse
from ROOT import TFile, TTree, AddressOf

class MGSimulation():

    def __init__(self):

        from ROOT import gROOT

        #Structure in the tree
        gROOT.ProcessLine(
        "struct h10 {\
        Float_t         Runevt;\
        Float_t       Xcor;\
        Float_t       Ycor;\
        Float_t       Zcor;\
        Float_t       Estep;\
        Float_t       Itra;\
        }" )

        #Root is dumb
        from ROOT import h10

        self.MGStruct = h10()


filename = 'FarFieldPointSource_1MeV.inc1.id1.sim'

sim = parse(filename)

print "There are {} events in this file.".format(len(sim.events))

#for event in sim.events[:10]:
#    print event.hits

#Make the tree
f = TFile('test.root', 'RECREATE')
t = TTree('h10','Simulations')

MGS = MGSimulation()

t.Branch('Runevt', AddressOf(MGS.MGStruct, 'Runevt'), 'Runevt/F')
t.Branch('Xcor', AddressOf(MGS.MGStruct, 'Xcor'), 'Xcor/F')
t.Branch('Ycor', AddressOf(MGS.MGStruct, 'Ycor'), 'Ycor/F')
t.Branch('Zcor', AddressOf(MGS.MGStruct, 'Zcor'), 'Zcor/F')
t.Branch('Estep', AddressOf(MGS.MGStruct, 'Estep'), 'Estep/F')

for event in sim.events:

    hits = zip(event.hits.detector,
               event.hits.x,
               event.hits.y,
               event.hits.z,
               event.hits.energy)
    
    for hit in hits:
        MGS.MGStruct.Runevt = float(event.id_trigger)
        MGS.MGStruct.Xcor = hit[1]
        MGS.MGStruct.Ycor = hit[2]
        MGS.MGStruct.Zcor = hit[3]
        MGS.MGStruct.Estep = hit[4]
        t.Fill()


f.Write()
f.Close()
