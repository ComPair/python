from EventViewer import parse
from ROOT import TFile, TTree, gROOT, AddressOf


filename = 'FarFieldPointSource_1MeV.inc1.id1.sim'

sim = parse(filename)

print "There are {} events in this file.".format(len(sim.events))

#for event in sim.events[:10]:
#    print event.hits

#Make the tree
f = TFile('test.root', 'RECREATE')
t = TTree('h10','Simulations')


#Structure in the tree
gROOT.ProcessLine(
"struct h10 {\
   Float_t         Runevt;\
   Float_t       Xcor;\
   Float_t       Ycor;\
   Float_t       Zcor;\
   Float_t       Estep;\
}" )

#Root is stupid.
from ROOT import h10

h = h10()
t.Branch('Runevt', AddressOf(h, 'Runevt'), 'Runevt/F')
t.Branch('Xcor', AddressOf(h, 'Xcor'), 'Xcor/F')
t.Branch('Ycor', AddressOf(h, 'Ycor'), 'Ycor/F')
t.Branch('Zcor', AddressOf(h, 'Zcor'), 'Zcor/F')
t.Branch('Estep', AddressOf(h, 'Estep'), 'Estep/F')

for event in sim.events:

    hits = zip(event.hits.detector,
               event.hits.x,
               event.hits.y,
               event.hits.z,
               event.hits.energy)
    
    for hit in hits:
        h.Runevt = float(event.id_trigger)
        h.Xcor = hit[1]
        h.Ycor = hit[2]
        h.Zcor = hit[3]
        h.Estep = hit[4]
        t.Fill()


f.Write()
f.Close()
