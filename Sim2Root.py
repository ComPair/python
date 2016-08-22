class MGSimulation():

    ''' This class is the description of the MGeant 
    ntuple in Root.  Root is dumb.'''
    
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


def cli():
    
    helpString = "This script takes an input Sim file from Cosima and spits out a Root\
    file which can be used in Alex MGeant analysis scripts.  Takes as\
    input the Simulation File name and the Root File name."

    import argparse

    parser = argparse.ArgumentParser(description=helpString)
    parser.add_argument("SimFile", help="Input Simulation File (from Cosima)")
    parser.add_argument("RootFile", help="Output Root File")

    args = parser.parse_args()
    
    from EventViewer import parse
    from ROOT import TFile, TTree, AddressOf

    sim = parse(args.SimFile)
    print "There are {} events in this file.".format(len(sim.events))

    #for event in sim.events[:10]:
    #    print event.hits

    #Make the tree
    f = TFile(args.RootFile, 'RECREATE')
    t = TTree('h10','Simulations')

    MGS = MGSimulation()

    t.Branch('Runevt', AddressOf(MGS.MGStruct, 'Runevt'), 'Runevt/F')
    t.Branch('Xcor', AddressOf(MGS.MGStruct, 'Xcor'), 'Xcor/F')
    t.Branch('Ycor', AddressOf(MGS.MGStruct, 'Ycor'), 'Ycor/F')
    t.Branch('Zcor', AddressOf(MGS.MGStruct, 'Zcor'), 'Zcor/F')
    t.Branch('Estep', AddressOf(MGS.MGStruct, 'Estep'), 'Estep/F')
    #t.Branch('Det', AddressOf(MGS.MGStruct, 'Det'), 'Det/F')
    
    for event in sim.events:

        hits = zip(event.hits.detector,
                   event.hits.x,
                   event.hits.y,
                   event.hits.z,
                   event.hits.energy)
    
        for hit in hits:
            MGS.MGStruct.Runevt = float(event.id_trigger)
            #MGS.MGStruct.Det = hit[0]
            MGS.MGStruct.Xcor = hit[1]
            MGS.MGStruct.Ycor = hit[2]
            MGS.MGStruct.Zcor = hit[3]
            MGS.MGStruct.Estep = hit[4]
            t.Fill()

    f.Write()
    f.Close()

if __name__ == '__main__': cli()
