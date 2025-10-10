#==============================================================================
import argparse
import sys
from collections import defaultdict
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("infile", type=str,
                    help="A blast result file. Tabular output")
parser.add_argument("outfile", type=str,
                    help="The name of the file")
# This checks if the user supplied any arguments. If not, help is printed.
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()
#==============================================================================
#Functions=====================================================================
#==============================================================================
def read_blastout(source):
    countdeleted = 0
    try:
        with open(source, "r") as infile:
            count =0
            blasthits = {}
            blastlist = []
            currentscaffold = None
            for line in infile:
                count +=1
                sys.stdout.write('%d\r' % (count))
                sys.stdout.flush()
                qseqid = line.rstrip().split()[0]
                if currentscaffold is None:
                    currentscaffold = qseqid
                    blastlist.append(line)
                else:
                    if currentscaffold == qseqid:
                        blastlist.append(line)
                    else:
                        blasthits, countdel = read_blastout_bitscore(blastlist, blasthits)
                        countdeleted = countdeleted + countdel
                        blastlist = []
                        currentscaffold = qseqid
                        blastlist.append(line)
            #final loop for final scaffold
            blasthits, countdel = read_blastout_bitscore(blastlist, blasthits)
            countdeleted = countdeleted + countdel

            print("Number of lines =", count)
            print("Number of blast hits filtered because no one tophit =", countdeleted)
        return blasthits
    except IOError:
        print("!----ERROR----!")
        print("File %s does not exit!" % source)
        sys.exit(1)
    except KeyboardInterrupt:
        sys.exit(1)

def read_blastout_bitscore(blastlist, blasthits):
    if len(blastlist) == 0:
        line = line.rstrip().split()
        qseqid = line[0]
        sseqid = line[1]
        pidentity = float(line[2])
        sstart = float(line[8])
        send = float(line[9])
        bitscore = float(line[11])
        blasthits[qseqid] = (sseqid, bitscore, pidentity, sstart, send)
    else:
        countdel = 0
        check = []
        for line in blastlist:
            line = line.rstrip().split()
            qseqid = line[0]
            sseqid = line[1]
            pidentity = float(line[2])
            sstart = float(line[8])
            send = float(line[9])
            bitscore = float(line[11])

            # Check for min 30% identity
            if pidentity > 30:
                if qseqid in blasthits:
                    if blasthits[qseqid][1] < bitscore:
                        blasthits[qseqid] = (sseqid, bitscore, pidentity, sstart, send)
                    elif blasthits[qseqid][1] == bitscore:
                        if blasthits[qseqid][2] < pidentity:
                            blasthits[qseqid] = (sseqid, bitscore, pidentity, sstart, send)
                        elif blasthits[qseqid][2] == pidentity:
                            info = [bitscore, pidentity]
                            check.append(info)
                else:
                    blasthits[qseqid] = (sseqid, bitscore, pidentity, sstart, send)
        #check that there arent two hits which have the highest bitscore and pidentity but are identical
        if float(len(check)) > 0:
            sortcheck = sorted(check, key=lambda t: t[0], reverse=True)
            bitscore = float(sortcheck[0][0])
            pidentity = float(sortcheck[0][1])
            if bitscore == blasthits[qseqid][1]:
                if pidentity == blasthits[qseqid][2]:
                    countdel += 1
                    del blasthits[qseqid]
    return blasthits, countdel

#==============================================================================
#Main==========================================================================
#==============================================================================
def main():
    blastoutput = read_blastout(args.infile)
    print("No. of genes with top blasthits on scaffolds =", len(blastoutput))

    print("Starting to run...")
    with open(args.outfile, "w") as outfile:
        for blast in blastoutput:
            # print(blastoutput[blast])
            outfile.write(blast)
            outfile.write(",")
            outfile.write(blastoutput[blast][0])
            outfile.write(",")
            outfile.write(str(blastoutput[blast][1]))
            outfile.write(",")
            outfile.write(str(blastoutput[blast][2]))
            outfile.write(",")
            outfile.write(str(blastoutput[blast][3]))
            outfile.write(",")
            outfile.write(str(blastoutput[blast][4]))
            outfile.write("\n")

if __name__ == '__main__':
    main()