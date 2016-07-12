import sys
import glob
import os

htcutstring = "Number of events after Ht cut found:"
higgscutstring = "Number of events after Higgs cut found:"
topcutstring = "Number of events after top cut found:"
eventstring = "Number of events found:"

htcutnum = 0
higgscutnum = 0
topcutnum = 0
eventnum = 0

thisdir = os.getcwd()

outfile = open( 'Cuts_by_mass.txt', 'w' )

for root, dirs, files in os.walk( thisdir ):
	for cdir in dirs:	
		for troot, tdir, tfiles in os.walk( cdir ):
	    		for tfile in tfiles:
        			if tfile.endswith(".txt"):
					filepath = cdir + '/' + tfile 
             				search = open( filepath )
					for line in search:
 						if htcutstring in line:
  							if htcutnum == 0: htcutnum = int( line.split()[-1] )
						if higgscutstring in line:
							if higgscutnum == 0: higgscutnum = int( line.split()[-1] )
						if topcutstring in line:
							if topcutnum == 0: topcutnum = int( line.split()[-1] ) 
						if eventstring in line:
                                                        if eventnum == 0: eventnum = int( line.split()[-1] )
    		out = 'For: ' + cdir + ' cut effs are : ' + str(htcutnum) + ' ' + str(higgscutnum) + ' ' + str(topcutnum) + " " + str(eventnum) + "\n"
    		print out
		outfile.write( out )
    		htcutnum = eventnum = higgscutnum = topcutnum = 0
outfile.close()


