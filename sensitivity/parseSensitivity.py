import glob
import pandas 

#Script to parse through sensitivity output log files
#Output is a table with BG rates, effecitve area, sensitivity for each event type,
#ordered by energy.
#Author: Henrike Fleischhack


#Set your input files here
files = glob.glob("Sensitivity/pixel_rightOrbit/*log")

data = []

for f in files:

    res = {}

    #File names should look like Sensitivity_Continuum_90-270_80keV_R1_bash.log
    #where the energy bin goes from 90 to 270 keV.
    E=int(f.split("_")[3].split("-")[0])
    print(f, E)
    res["E"] = E

    with open(f, "r") as fi:
        
        for line in fi.readlines():

            if "Best achievable sensitivity - untracked compton" in line:
                current = "UC"

            if "Best achievable sensitivity - tracked compton" in line:
                current = "TC"

            if "Best achievable sensitivity - pairs" in line:
                current = "P"

            if " -> " not in line: 
                continue
                
            if "EffectiveArea" in line:
                EA = line.split()[2]
                res[f"Aeff_{current}"] = EA

            if "Sensitivity" in line:
                sens = line.split()[2]
                res[f"Sens_{current}"] = sens
            if "background" in line:
                rate = line.split()[4]
                res[f"BG_{current}"] = rate

    data.append(res)

#dataframe magic
cols = [ f"{i}_{j}" for j in ["UC", "TC", "P" ] for i in ["BG", "Aeff", "Sens"] ]
df = pandas.DataFrame( data ).set_index("E").sort_index()[cols].fillna("")

print(df)

#Set your output filename here
df.to_csv("pixel_rightOrbit_sens.csv")
