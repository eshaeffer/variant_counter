#! /bin/python3


# I started with a simple CommandLine interface where the user would enter the location of the vcfs, but this does
# not work so well with argon job submission so you always had to run it on a login node. To fix this, I'm going to
# use an argument parser that will read in the file and store it in a global variable.
# add file location


import os, argparse

parser = argparse.ArgumentParser()
parser.add_argument("-input",
                    help="text file containing the paths to the vcf files with each one being on seperate lines")
parser.add_argument("-output", help="location and name of textfile to be created containing the counts")
parser.add_argument("-option", default="rawCalls",
                    help="the option for running the script: rawCalls, mafCounting, vepCounting")
parser.add_argument("-k", default=False, action='store_true',
                    help="flag to activate countin of MAF in an alternative way (uses VEP annotations instead of bcftools)")
parser.add_argument("-outputDir", default="",
                    help="input directory to be used to output of rare and deletorious vars to a vcf file")

#test git
args = parser.parse_args()

# open file and get paths, storing them in filePaths
file = open(args.input, "r")

files = []

line = file.readline()

while (line != ""):
    files.append(line.rstrip())
    line = file.readline()
file.close()

# store other arguments as global variables
outputLoc = args.output
option = args.option
altCount = args.k
outputRare = args.outputDir


class rawVariantCounter:

    def __init__(self):
        self.sampleIds = []
        self.numberOfVars = []

    def countRawCalls(self):
        for x in files:
            numVar = 0
            f = open(x, 'r')
            line = f.readline()
            passedHeaders = False
            while (line != ""):
                if (passedHeaders):
                    numVar += 1
                else:
                    if "#CHROM" in line:
                        passedHeaders = True
                        temp = line.split()
                        self.sampleIds.append(temp[9])
                        print(temp[9])
                line = f.readline()

            self.numberOfVars.append(numVar)

    def outputRawCountTextFile(self):
        IDs = self.sampleIds
        counts = self.numberOfVars
        f = open(outputLoc, "w+")
        f.write("no.\tSample_ID\tnumVars\n")
        for num in range(0, len(IDs)):
            count = str(counts[num])
            number = str(num + 1)
            data = number + "\t" + IDs[num] + "\t" + count + "\n"
            data = data.format()
            f.write(data)
        f.close()


class afOnlyVariantCounter:
    def __init__(self):
        self.sampleIds = []
        self.numberOfVars = []
        # Cut-off: 0.01
        self.numberOfVarsGoodAf = []
        self.numberOfNoAfVars = []

    def countAfOnlyVariants(self):
        for x in files:
            numVarSansAf = 0
            numVarFilt = 0
            numVar = 0
            f = open(x, 'r')
            line = f.readline()
            passedHeaders = False
            # The reason why passedHeaders is neccessary is because a VCF file looks like:
            #	<headers that contain options that describe the file>
            #	<headers that contain options that describe the file>
            #	                          .
            #	                          .
            #                       #CHROM ....... SampleID
            #                        Variant information
            #                                 .
            #                                 .
            #                        Variant information
            # None of the information in headers is pertinant and must be ignored
            while (line != ""):
                # Some variants don't have GnomAD scores and need to be included in the filtered file. This is a
                # a boolean that will switch back and for each line written, and if while searching for a AF
                # score -1 there happensto be no GnomAD score, this boolean will simply remain false and cause
                # line to be written.
                if (passedHeaders):
                    GnomAdScorePresent = 0
                    # Each line past header is a different variant (variant line)
                    # so numVar is incremented for each line.
                    numVar += 1
                    # The filter information is a mess of different scores
                    # and values seperated by semicolons, and is the 8th column
                    # in each variant line.
                    temp = line.split()
                    filterInfo = temp[7].split(";")
                    # This checks each member filter info looking for the keyword tagging the
                    # MAF scores. If it finds it, it replaces it with empty string and
                    # casts the remainder to a float. If the float is less that 0.01
                    # with 9 degrees of precision then it increments the number of
                    # variants with MAF <= 0.01. The keyword can easily be replaced.
                    # simple change the if statement and 1st argument of replace method
                    # the the files' tag for MAF scores.
                    for t in filterInfo:
                        if "bcfGnomAD_AF=" in t:
                            GnomAdScorePresent = 1
                            t = t.replace("bcfGnomAD_AF=", "")
                            try:
                                temp = float(t)
                                if (temp < 0.01000000001):
                                    numVarFilt += 1
                            except:
                                t = t.split(",")
                                try:
                                    temp = float(t[0])
                                    if (temp < 0.010000001):
                                        numVarFilt += 1
                                except:
                                    GnomAdScorePresent = 0

                    if not GnomAdScorePresent:
                        numVarSansAf += 1

                else:
                    # The first line after the headers are finished
                    # always contains #CHROM to my knowlage. This
                    # could very easily be modified if this is not
                    # the case.
                    if ("#CHROM" in line):
                        passedHeaders = True
                        # convieniently, the first line after the headers
                        # also contains the sample id as the last column
                        # in a tab seperated line.
                        temp = line.split()
                        sampleId = temp[9]
                line = f.readline()

            self.sampleIds.append(sampleId)
            self.numberOfVars.append(numVar)
            self.numberOfVarsGoodAf.append(numVarFilt)
            self.numberOfNoAfVars.append(numVarSansAf)
            f.close()

    def outputVariantNumberAfOnly(self):
        f = open(outputLoc, "w+")
        # Writes data in a formatted fashion to the file location specified by user.
        f.write("no.\tSampleID\tNumber_of_Variants\tNum_Variants_MAF<=0.01\tNum_Variants_noMAF\n")
        for id in range(0, len(self.sampleIds)):
            sampleID = self.sampleIds[id]
            numberOfVarsForSample = str(self.numberOfVars[id])
            numberOfPassedVars = str(self.numberOfVarsGoodAf[id])
            numberOfNoAfVars = str(self.numberOfNoAfVars[id])
            number = str(id + 1)
            data = number + "\t" + sampleID + "\t" + \
                   numberOfVarsForSample + "\t" + numberOfPassedVars + \
                   "\t" + numberOfNoAfVars + "\n"
            data = data.format()
            f.write(data)
        f.close()


class VepCounter:
    def __init__(self):
        self.sampleIds = []
        self.numVars = []
        self.varsNoAF = []
        # MAF <=0.01
        self.rareVars = []
        # CADD >= 20
        self.deletoriousVars = []
        # CADD >= 15
        self.maybeDeletoriousVars = []
        # CADD >= 20 & MAF <= 0.01
        self.rareDeletoriousVars = []
        # CADD >= 15 & MAF <= 0.01
        self.rarePossibleDelVars = []
        # no MAF cadd >= 20
        self.deletoriousNoMaf = []
        # no MAF cadd >= 15
        self.possiblyDelNoMaf = []
        # no CADD
        self.noCADD = []

        self.outputPaths20L = []
        self.outputPaths15L = []
        self.outputPaths20N = []
        self.outputPaths15N = []
        self.outputPathsMaf = []

        if os.path.isdir(args.outputDir):

            outputRare = args.outputDir
            inputNoTxt = args.input.replace(".txt", "")
            outputRare = os.path.join(outputRare, inputNoTxt)

            outpt20L = os.path.join(outputRare, "20CADD_lowMAf")
            outpt15L = os.path.join(outputRare, "15CADD_lowMAf")
            outpt20N = os.path.join(outputRare, "CADD20")
            outpt15N = os.path.join(outputRare, "cADD15")
            outptMaf = os.path.join(outputRare, "lowMAF")

            if not os.path.isdir(outpt20L):
                os.makedirs(outpt20L)
            if not os.path.isdir(outpt20N):
                os.makedirs(outpt20N)
            if not os.path.isdir(outpt15L):
                os.makedirs(outpt15L)
            if not os.path.isdir(outpt15N):
                os.makedirs(outpt15N)
            if not os.path.isdir(outptMaf):
                os.makedirs(outptMaf)

            for filePath in files:
                tmp = os.path.basename(filePath)
                tmp = tmp.replace(".vcf", "")

                outpt20Ltmp = os.path.join(outpt20L, tmp)
                outpt15Ltmp = os.path.join(outpt15L, tmp)
                outpt20Ntmp = os.path.join(outpt20N, tmp)
                outpt15Ntmp = os.path.join(outpt15N, tmp)
                outptMaftmp = os.path.join(outptMaf, tmp)

                tmp = os.path.basename(filePath)
                tmp = tmp.replace(".vcf", "")

                outpt20Ltmp += ".20Cadd.vcf"
                outpt20Ntmp += ".20CaddNoMAF.vcf"
                outpt15Ltmp += ".15CaddLowMAF.vcf"
                outpt15Ntmp += ".15Cadd.vcf"
                outptMaftmp += ".lowMAF.vcf"

                self.outputPaths20L.append(outpt20Ltmp)
                self.outputPaths15L.append(outpt15Ltmp)
                self.outputPaths20N.append(outpt20Ntmp)
                self.outputPaths15N.append(outpt15Ntmp)
                self.outputPathsMaf.append(outptMaftmp)


        else:
            print("Invalid path specified or path was not entered. Defaulting to count only mode.")

    def countVarsVEP(self):
        fileCount = 0

        pth20L = self.outputPaths20L
        pth20N = self.outputPaths20N
        pth15L = self.outputPaths15L
        pth15N = self.outputPaths15N
        pthMaf = self.outputPathsMaf

        for x in files:
            # set/reset local variables to 0
            numVarSansAf = 0
            numVarFilt = 0
            numG20 = 0
            numG15 = 0
            numG20L = 0
            numG15L = 0
            numG20N = 0
            numG15N = 0
            numNoCadd = 0
            numVar = 0

            try:
                f = open(x, 'r')
            except:
                print(x + " could not be opened")
            try:
                output20L = open(pth20L[fileCount], 'w+')
            except:
                print(pth20L[fileCount] + " could not be created")
            try:
                output20N = open(pth20N[fileCount], 'w+')
            except:
                print(pth20N[fileCount] + " could not be created")
            try:
                output15N = open(pth15N[fileCount], 'w+')
            except:
                print(pth15N[fileCount] + " could not be created")
            try:
                output15L = open(pth15L[fileCount], 'w+')
            except:
                print(pth15L[fileCount] + " could not be created")
            try:
                outputMaf = open(pthMaf[fileCount], 'w+')
            except:
                print(pthMaf[fileCount] + " could not be created")
            try:
                line = f.readline()
            except:
                line = ""

            passedHeaders = False

            while (line != ""):
                if (passedHeaders):
                    numVar += 1
                    temp = line.split()
                    try:
                        filterInfo = temp[7].split(";")
                        cadd = filterInfo[len(filterInfo) - 1].split("|")
                    except:
                        filterInfo = "junk"
                        cadd = "junk"
                    G = False

                    GnomAdScorePresent = 0
                    if altCount:
                        t = cadd[len(cadd) - 1]
                        try:
                            t = float(t)
                            if (t < 0.010000001):
                                G = True
                            GnomAdScorePresent = 1
                        except:
                            try:
                                t = t.replace(",", ".")
                                t = float(t)
                                if (t < 0.010000001):
                                    G = True
                                    GnomAdScorePresent = 1
                            except:
                                pass
                    else:
                        for t in filterInfo:
                            if "bcfGnomAD_AF=" in t:
                                GnomAdScorePresent = 1
                                t = t.replace("bcfGnomAD_AF=", "")
                                try:
                                    temp = float(t)
                                    if (temp < 0.01000000001):
                                        G = True
                                except:
                                    t = t.split(",")
                                    try:
                                        temp = float(t[0])
                                        if (temp < 0.010000001):
                                            G = True
                                    except:
                                        GnomAdScorePresent = 0

                    caddPresent = 0
                    try:
                        caddScore = float(cadd[len(cadd) - 6])
                        caddPresent = 1
                    except:
                        pass

                    if (caddPresent):
                        if (caddScore > 19.9999):
                            numG20 += 1
                            try:
                                output20N.write(line)
                            except:
                                pass
                        elif (caddScore > 14.999999):
                            numG15 += 1
                            try:
                                output15N.write(line)
                            except:
                                pass
                        if (caddScore > 19.9999 and GnomAdScorePresent and G):
                            numG20L += 1
                            try:
                                output20L.write(line)
                            except:
                                pass
                        elif (caddScore > 14.9999 and GnomAdScorePresent and G):
                            numG15L += 1
                            try:
                                output15L.write(line)
                            except:
                                pass
                        if (caddScore > 19.99999 and not GnomAdScorePresent):
                            numG20N += 1
                        elif (caddScore > 14.9999 and not GnomAdScorePresent):
                            numG15N += 1
                    if G:
                        numVarFilt += 1
                        try:
                            outputMaf.write(line)
                        except:
                            pass
                    if (not GnomAdScorePresent):
                        numVarSansAf += 1
                    if (not caddPresent):
                        numNoCadd += 1
                else:
                    output20N.write(line)
                    output15N.write(line)
                    output20L.write(line)
                    output15L.write(line)
                    outputMaf.write(line)
                    if ("#CHROM" in line):
                        passedHeaders = True
                        temp = line.split()
                        sampleId = temp[9]
                        print("Reading vepped variants of the current sample: " + sampleId)
                line = f.readline()

            self.sampleIds.append(sampleId)
            self.numVars.append(numVar)
            self.varsNoAF.append(numVarSansAf)
            self.rareVars.append(numVarFilt)
            self.deletoriousVars.append(numG20)
            self.maybeDeletoriousVars.append(numG15)
            self.rareDeletoriousVars.append(numG20L)
            self.rarePossibleDelVars.append(numG15L)
            self.deletoriousNoMaf.append(numG20N)
            self.possiblyDelNoMaf.append(numG15N)
            self.noCADD.append(numNoCadd)
            fileCount += 1
            f.close()
            output20N.close()
            output20L.close()
            output15L.close()
            output15N.close()
            outputMaf.close()

    def outputVepCount(self):
        f = open(outputLoc, "w+")
        f.write(
            "no.\tSampleID\tNumVar_Gatk\tNumVarGatk_MAF<=0.01\tNumVarGatk_noMAF\t20Cadd\t20CaddL\t20CaddN\t15Cadd\t15CaddL\t15CaddN\tnoCadd\n")
        for id in range(0, len(self.sampleIds)):
            sampleID = self.sampleIds[id]
            numberOfVarsForSampleGatk = str(self.numVars[id])
            numberOfPassedVarsGatk = str(self.rareVars[id])
            numberOfNoAfVarsGatk = str(self.varsNoAF[id])
            Cadd20 = str(self.deletoriousVars[id])
            Cadd20L = str(self.rareDeletoriousVars[id])
            Cadd20N = str(self.deletoriousNoMaf[id])
            Cadd15 = str(self.maybeDeletoriousVars[id])
            Cadd15L = str(self.rarePossibleDelVars[id])
            Cadd15N = str(self.possiblyDelNoMaf[id])
            noCadd = str(self.noCADD[id])
            number = str(id + 1)
            data = number + "\t" + sampleID + "\t" + \
                   numberOfVarsForSampleGatk + "\t" + numberOfPassedVarsGatk + \
                   "\t" + numberOfNoAfVarsGatk + "\t" + \
                   Cadd20 + "\t" + Cadd20L + "\t" + Cadd20N + "\t" + \
                   Cadd15 + "\t" + Cadd15L + "\t" + Cadd15N + "\t" + noCadd + "\n"
            data = data.format()
            f.write(data)
        f.close()


if __name__ == "__main__":
    if option == "rawCalls":
        variants = rawVariantCounter()
        variants.countRawCalls()
        variants.outputRawCountTextFile()
    elif option == "mafCounting":
        variants = afOnlyVariantCounter()
        variants.countAfOnlyVariants()
        variants.outputVariantNumberAfOnly()
    elif option == "vepCounting":
        variants = VepCounter()
        variants.countVarsVEP()
        variants.outputVepCount()
