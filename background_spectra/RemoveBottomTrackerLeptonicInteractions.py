import gzip
import re
import os
import glob

#trafile = gzip.open('TraFiles_100420/R1_Files/TrappedLeptonicBackground/TrappedLeptonicBackground.p2.inc30.id1.tra.gz', 'rb')
#newtrafile = gzip.open('TrappedLeptonicBackground_mod.p2.inc30.id1.tra.gz', 'wb')
			
def RemoveEvents(filename, newtrafilename):

	trafile = gzip.open(filename, 'rb')
	print('Reading in file ', filename) 
	
	newtrafile = gzip.open(newtrafilename, 'wb')


	#Write the first 7 lines of the tra file into the new condensed file
	for x in range(7):
		line = trafile.readline()
		newtrafile.write(line)

	#Read the whole file, and append the event details to new_event. If the position is not consistent with the last tracker layer, save the event into the new trafile
	while True:
		line = trafile.readline()
		if not line:
			break
		#Find the start of an event with the SE flag
		if re.match('SE', line.decode('utf-8')):
			new_event = []
			#append SE to be the first element in the new_event list
			new_event.append(line)

			#Read the next line and see if it matches with a CO event (since these are the only ones in the TrappedLeptonic file other than unknowns)
			line = trafile.readline()
			if re.match('ET CO', line.decode('utf-8')):
				new_event.append(line)
			

	#			line = trafile.readline()
	#			if re.match('ID 5349', line.decode('utf-8')):
	#				print('found ID 5349')
	#				new_event.append(line)

				#If it's a compton event, append the next 9 lines of the event details to the new_event list to get to the position info
				for x in range(7):
					line = trafile.readline()
					new_event.append(line)
				#Some events have the TQ flag which makes them one line longer, so we have to check for that
				if re.match('TQ', line.decode('utf-8')):
					for x in range(3):
						line = trafile.readline()
						new_event.append(line)
				else: 
					for x in range(2):
						line = trafile.readline()
						new_event.append(line)

							

				#Read the CH line and see if it's consistent with the first event being in the bottom layer of the tracker
				line = trafile.readline()
				if re.match('CH 0 -?\d+\.\d+ -?\d+\.\d+ 0.75', line.decode('utf-')):
					#print('Got a bad one!')
					pass

				#If the event isn't consistent, then write it to the new tra file
				else:	
					new_event.append(line)
					line = trafile.readline()
					while re.match('CH \d ', line.decode('utf-8')):
						new_event.append(line)
						line = trafile.readline()
					#print(new_event)
					for i in new_event:
						newtrafile.write(i)
					
		#If you've reached the last event, copy over the final lines of the file to include the FT info
		if re.match('EN', line.decode('utf-8')):
			newtrafile.write(line)
			while True:
				line = trafile.readline()
				newtrafile.write(line)
				if not line:
					break

	

	trafile.close()
	newtrafile.close()

	return
			
###########################################333

	

	
#Change the path here!
#path = 'TraFiles_040520/R1_Files/TrappedLeptonicBackground/'
path = '/data/slag2/hfleisc1/amego_sims/simfiles/pixel/background/TrappedLeptonicBackground/'
path = '/data/slag2/hfleisc1/amego_sims/simfiles/base/background/TrappedLeptonicBackground/'
path = '/data/slag2/hfleisc1/amego_sims/simfiles/pixel/BG_above_0.01MeV/TrappedLeptonicBackground/'
#for filename in glob.glob(os.path.join(path, 'TrappedLeptonicBackground.p*.inc*.tra.gz')):
for filename in os.listdir(path):
        if os.path.isfile(os.path.join(path, filename)):
                if re.match('TrappedLeptonicBackground.p\d.inc\d+.id1.tra.gz', filename):
                        newtrafilename = filename.replace('TrappedLeptonicBackground','TrappedLeptonicBackground_old')
                        print (filename, newtrafilename)
                        os.rename( os.path.join(path, filename), os.path.join(path, newtrafilename))
                        RemoveEvents(os.path.join(path, newtrafilename), os.path.join(path,filename))



