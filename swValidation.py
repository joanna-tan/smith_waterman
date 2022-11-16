import sys
import matplotlib.pyplot as plt
from termcolor import colored

#Keep this code in the "bigp/lab9" directory and keep the CSV files in the "bigp" directory
#Run command:
#from "bigp" directory: python3 lab9/swValidation.py 'swCLEAN.txt' ['Changed SW #1 file name','Changed SW #2 file name'...]


def getData(score_CSVfile):
    file_in = open(score_CSVfile, "r")
    data = {"score":[],"start_col":[],"start_row":[],"end_col":[],"end_row":[],"pIDtot":[],"pSIMtot":[],
            "pBOTHtot":[],"pGAPtot":[],"pIDseed":[],"pSIMseed":[],"pBOTHseed":[],"pGAPseed":[],"pIDmid":[],
            "pSIMmid":[],"pBOTHmid":[],"pGAPmid":[],"ntGAPtot":[],"nqGAPtot":[],"ntGAPseed":[],"nqGAPseed":[],
            "ntGAPmid":[],"nqGAPmid":[],"t_str": [],"q_str": [],"align_str": [],"matchLenQ":[],"matchLenT":[]}

    for line in file_in:
        l = line.strip().split(",")
        data["start_col"].append(int(l[0]))
        data["score"].append(float(l[1]))
        data["start_row"].append(int(l[4]))
        data["end_col"].append(int(l[6]))
        data["end_row"].append(int(l[7]))
        data["pIDtot"].append(float(l[8]))
        data["pSIMtot"].append(float(l[9]))
        data["pBOTHtot"].append(float(l[10]))
        data["pGAPtot"].append(float(l[11]))
        data["pIDseed"].append(float(l[12]))
        data["pSIMseed"].append(float(l[13]))
        data["pBOTHseed"].append(float(l[14]))
        data["pGAPseed"].append(float(l[15]))
        data["pIDmid"].append(float(l[16]))
        data["pSIMmid"].append(float(l[17]))
        data["pBOTHmid"].append(float(l[18]))
        data["pGAPmid"].append(float(l[19]))
        data["ntGAPtot"].append(int(l[20]))
        data["nqGAPtot"].append(int(l[21]))
        data["ntGAPseed"].append(int(l[22]))
        data["nqGAPseed"].append(int(l[23]))
        data["ntGAPmid"].append(int(l[24]))
        data["nqGAPmid"].append(int(l[25]))
        data["t_str"].append(l[2])
        data["q_str"].append(l[5])
        data["align_str"].append(l[3])
        data["matchLenQ"].append(data["end_row"][-1]-data["start_row"][-1])
        data["matchLenT"].append(data["end_col"][-1] - data["start_col"][-1])
    return data

def validation(swCLEAN, swUSlist):
    swCleanData = getData(swCLEAN)
    swUsData = []
    maxMatches = len(swCleanData["score"])
    for swUS in swUSlist:
        swUsData.append(getData(swUS))
        if len(swUsData[-1]["score"]) < maxMatches:     #Due to the possibility of tied scores, there is not the same number of matches
            maxMatches = len(swUsData[-1]["score"])     #in each data set so this picks out the number of matches in the smallest set
    nDataSet = len(swUsData)+1                          #Used for opacities of the histograms (+1 for swCleanData)

############################################### PRINT ALIGNMENT DIAGRAMS ###############################################
    print("Comparing the top "+str(maxMatches)+" matches generated by "+str(nDataSet)+" different SW algorithms...")
    print(colored("Matches from the unaltered SW algorithm ("+swCLEAN+"):","green", attrs=['bold','underline']))
    scoreNumber = 1
    for i in range(maxMatches):
        print(colored("Match #{}, score={:3.2f}".format(scoreNumber, swCleanData["score"][i]),"cyan", attrs=['bold']))
        print("{:4.0f}: {} :{:<4.0f} {}".format(swCleanData["start_col"][i], swCleanData["t_str"][i], swCleanData["end_col"][i],
                                            ">>> target"))
        print("      {}".format(swCleanData["align_str"][i]))
        print("{:4.0f}: {} :{:<4.0f} {}\n".format(swCleanData["start_row"][i], swCleanData["q_str"][i], swCleanData["end_row"][i],
                                              "<<< query"))
        scoreNumber = scoreNumber + 1
    for j in range(nDataSet-1):
        print(colored("-------------------------------------------------------------------------------------------",'red',attrs=['bold']))
        print(colored("Matches from the altered SW algorithm #"+str(j+1)+" ("+swUSlist[j]+"):", "red", attrs=['bold','underline']))
        scoreNumber = 1
        for i in range(maxMatches):
            print(colored("Match #{}, score={:3.2f}".format(scoreNumber, swUsData[j]["score"][i]), "cyan", attrs=['bold']))
            print("{:4.0f}: {} :{:<4.0f} {}".format(swUsData[j]["start_col"][i], swUsData[j]["t_str"][i],
                                                swUsData[j]["end_col"][i],
                                                ">>> target"))
            print("      {}".format(swUsData[j]["align_str"][i]))
            print("{:4.0f}: {} :{:<4.0f} {}\n".format(swUsData[j]["start_row"][i], swUsData[j]["q_str"][i],
                                                  swUsData[j]["end_row"][i],
                                                  "<<< query"))
            scoreNumber = scoreNumber + 1

    percentbins = [i for i in range(102) if i % 2 == 0]
################################ TOTAL %ID, %SIM, %BOTH, %GAP (and ntGAP, nqGAP) PLOTS #################################
    #Histogram for %ID and %Sim
    fig, axs = plt.subplots(2)
    fig.suptitle("Total Percent Identity and Total Percent Similarity")
    ax1 = axs[0]
    ax2 = axs[1]
    ax1.hist(swCleanData["pIDtot"][0:maxMatches], bins=percentbins, color='k', alpha=1 / nDataSet,
             label='Original SW %ID')
    iter = 0
    for swRun in swUsData:
        iter += 1
        ax1.hist(swRun["pIDtot"][0:maxMatches], bins=percentbins, alpha=1 / nDataSet,
                 label='Altered SW' + str(iter) + ' %ID')
    ax1.legend()
    ax1.set_xlabel("% Identity")
    ax1.set_ylabel("Frequency")
    ax1.set_xlim(0, 102)
    ax2.hist(swCleanData["pSIMtot"][0:maxMatches], bins=percentbins, color='k', alpha=1 / nDataSet, label='Original SW %Sim')
    iter = 0
    for swRun in swUsData:
        iter += 1
        ax2.hist(swRun["pSIMtot"][0:maxMatches], bins=percentbins, alpha=1 / nDataSet, label='Altered SW' + str(iter) + ' %Sim')
    ax2.legend()
    ax2.set_xlabel("% Similarity")
    ax2.set_ylabel("Frequency")
    ax2.set_xlim(0, 102)
    plt.show()

    #Histogram for %Both and %Gap
    fig, axs = plt.subplots(2)
    fig.suptitle("Total Percent Identity and Similar and Total Percent Gap")
    ax1 = axs[0]
    ax2 = axs[1]
    ax1.hist(swCleanData["pBOTHtot"][0:maxMatches], bins=percentbins, color='k', alpha=1 / nDataSet,
             label='Original SW %Both')
    iter = 0
    for swRun in swUsData:
        iter += 1
        ax1.hist(swRun["pBOTHtot"][0:maxMatches], bins=percentbins, alpha=1 / nDataSet,
                 label='Altered SW' + str(iter) + ' %Both')
    ax1.legend()
    ax1.set_xlabel("% Identity and Similar")
    ax1.set_ylabel("Frequency")
    ax1.set_xlim(0, 102)
    ax2.hist(swCleanData["pGAPtot"][0:maxMatches], bins=percentbins, color='k', alpha=1 / nDataSet, label='Original SW %Gap')
    iter = 0
    for swRun in swUsData:
        iter += 1
        ax2.hist(swRun["pGAPtot"][0:maxMatches], bins=percentbins, alpha=1 / nDataSet, label='Altered SW' + str(iter) + ' %Gap')
    ax2.legend()
    ax2.set_xlabel("% Gap")
    ax2.set_ylabel("Frequency")
    ax2.set_xlim(0, 102)
    plt.show()

    # Histogram for ntGAP and nqGAP
    nGAPbins = [i/2 for i in range(52)]
    fig, axs = plt.subplots(2)
    fig.suptitle("Total Number of Gaps in either miRNA or 3'UTR sequence")
    ax1 = axs[0]
    ax2 = axs[1]
    ax1.hist(swCleanData["ntGAPtot"][0:maxMatches], bins=nGAPbins, color='k', alpha=1 / nDataSet,
             label="Original SW 3'UTR Gaps")
    iter = 0
    for swRun in swUsData:
        iter += 1
        ax1.hist(swRun["ntGAPtot"][0:maxMatches], bins=nGAPbins, alpha=1 / nDataSet,
                 label="Altered SW" + str(iter) + " 3'UTR Gaps")
    ax1.legend()
    ax1.set_xlabel("Number of 3'UTR Gaps")
    ax1.set_ylabel("Frequency")
    ax1.set_xlim(0, 25)
    ax2.hist(swCleanData["nqGAPtot"][0:maxMatches], bins=nGAPbins, color='k', alpha=1 / nDataSet,
             label="Original SW miRNA Gaps")
    iter = 0
    for swRun in swUsData:
        iter += 1
        ax2.hist(swRun["nqGAPtot"][0:maxMatches], bins=nGAPbins, alpha=1 / nDataSet,
                 label='Altered SW' + str(iter) + ' miRNA Gaps')
    ax2.legend()
    ax2.set_xlabel("Number of miRNA Gaps")
    ax2.set_ylabel("Frequency")
    ax2.set_xlim(0, 25)
    plt.show()

################################## SEED %ID, %SIM, %BOTH, %GAP PLOTS ##################################
    # Histogram for %ID and %Sim
    fig, axs = plt.subplots(2)
    fig.suptitle("Seed Percent Identity and Seed Percent Similarity")
    ax1 = axs[0]
    ax2 = axs[1]
    ax1.hist(swCleanData["pIDseed"][0:maxMatches], bins=percentbins, color='k', alpha=1 / nDataSet,
             label='Original SW %ID')
    iter = 0
    for swRun in swUsData:
        iter += 1
        ax1.hist(swRun["pIDseed"][0:maxMatches], bins=percentbins, alpha=1 / nDataSet,
                 label='Altered SW' + str(iter) + ' %ID')
    ax1.legend()
    ax1.set_xlabel("% Identity")
    ax1.set_ylabel("Frequency")
    ax1.set_xlim(0, 102)
    ax2.hist(swCleanData["pSIMseed"][0:maxMatches], bins=percentbins, color='k', alpha=1 / nDataSet, label='Original SW %Sim')
    iter = 0
    for swRun in swUsData:
        iter += 1
        ax2.hist(swRun["pSIMseed"][0:maxMatches], bins=percentbins, alpha=1 / nDataSet, label='Altered SW' + str(iter) + ' %Sim')
    ax2.legend()
    ax2.set_xlabel("% Similarity")
    ax2.set_ylabel("Frequency")
    ax2.set_xlim(0, 102)
    plt.show()

    # Histogram for %Both and %Gap
    fig, axs = plt.subplots(2)
    fig.suptitle("Seed Percent Identity and Similar and Seed Percent Gap")
    ax1 = axs[0]
    ax2 = axs[1]
    ax1.hist(swCleanData["pBOTHseed"][0:maxMatches], bins=percentbins, color='k', alpha=1 / nDataSet,
             label='Original SW %Both')
    iter = 0
    for swRun in swUsData:
        iter += 1
        ax1.hist(swRun["pBOTHseed"][0:maxMatches], bins=percentbins, alpha=1 / nDataSet,
                 label='Altered SW' + str(iter) + ' %Both')
    ax1.legend()
    ax1.set_xlabel("% Identity and Similar")
    ax1.set_ylabel("Frequency")
    ax1.set_xlim(0, 102)
    ax2.hist(swCleanData["pGAPseed"][0:maxMatches], bins=percentbins, color='k', alpha=1 / nDataSet, label='Original SW %Gap')
    iter = 0
    for swRun in swUsData:
        iter += 1
        ax2.hist(swRun["pGAPseed"][0:maxMatches], bins=percentbins, alpha=1 / nDataSet, label='Altered SW' + str(iter) + ' %Gap')
    ax2.legend()
    ax2.set_xlabel("% Gap")
    ax2.set_ylabel("Frequency")
    ax2.set_xlim(0, 102)
    plt.show()

    # Histogram for ntGAP and nqGAP
    fig, axs = plt.subplots(2)
    fig.suptitle("Number of Gaps in either miRNA or 3'UTR sequence in Seed Matching Region")
    ax1 = axs[0]
    ax2 = axs[1]
    ax1.hist(swCleanData["ntGAPseed"][0:maxMatches], bins=nGAPbins, color='k', alpha=1 / nDataSet,
             label="Original SW 3'UTR Gaps")
    iter = 0
    for swRun in swUsData:
        iter += 1
        ax1.hist(swRun["ntGAPseed"][0:maxMatches], bins=nGAPbins, alpha=1 / nDataSet,
                 label="Altered SW" + str(iter) + " 3'UTR Gaps")
    ax1.legend()
    ax1.set_xlabel("Number of 3'UTR Gaps")
    ax1.set_ylabel("Frequency")
    ax1.set_xlim(0, 25)
    ax2.hist(swCleanData["nqGAPseed"][0:maxMatches], bins=nGAPbins, color='k', alpha=1 / nDataSet,
             label="Original SW miRNA Gaps")
    iter = 0
    for swRun in swUsData:
        iter += 1
        ax2.hist(swRun["nqGAPseed"][0:maxMatches], bins=nGAPbins, alpha=1 / nDataSet,
                 label='Altered SW' + str(iter) + ' miRNA Gaps')
    ax2.legend()
    ax2.set_xlabel("Number of miRNA Gaps")
    ax2.set_ylabel("Frequency")
    ax2.set_xlim(0, 25)
    plt.show()

################################## 12-17 %ID, %SIM, %BOTH, %GAP PLOTS ##################################
    # Histogram for %ID and %Sim
    fig, axs = plt.subplots(2)
    fig.suptitle("12-17 Percent Identity and 12-17 Percent Similarity")
    ax1 = axs[0]
    ax2 = axs[1]
    ax1.hist(swCleanData["pIDmid"][0:maxMatches], bins=percentbins, color='k', alpha=1 / nDataSet,
             label='Original SW %ID')
    iter = 0
    for swRun in swUsData:
        iter += 1
        ax1.hist(swRun["pIDmid"][0:maxMatches], bins=percentbins, alpha=1 / nDataSet,
                 label='Altered SW' + str(iter) + ' %ID')
    ax1.legend()
    ax1.set_xlabel("% Identity")
    ax1.set_ylabel("Frequency")
    ax1.set_xlim(0, 102)
    ax2.hist(swCleanData["pSIMmid"][0:maxMatches], bins=percentbins, color='k', alpha=1 / nDataSet, label='Original SW %Sim')
    iter = 0
    for swRun in swUsData:
        iter += 1
        ax2.hist(swRun["pSIMmid"][0:maxMatches], bins=percentbins, alpha=1 / nDataSet, label='Altered SW' + str(iter) + ' %Sim')
    ax2.legend()
    ax2.set_xlabel("% Similarity")
    ax2.set_ylabel("Frequency")
    ax2.set_xlim(0, 102)
    plt.show()

    # Histogram for %Both and %Gap
    fig, axs = plt.subplots(2)
    fig.suptitle("12-17 Percent Identity and Similar and 12-17 Percent Gap")
    ax1 = axs[0]
    ax2 = axs[1]
    ax1.hist(swCleanData["pBOTHmid"][0:maxMatches], bins=percentbins, color='k', alpha=1 / nDataSet,
             label='Original SW %Both')
    iter = 0
    for swRun in swUsData:
        iter += 1
        ax1.hist(swRun["pBOTHmid"][0:maxMatches], bins=percentbins, alpha=1 / nDataSet,
                 label='Altered SW' + str(iter) + ' %Both')
    ax1.legend()
    ax1.set_xlabel("% Identity and Similar")
    ax1.set_ylabel("Frequency")
    ax1.set_xlim(0, 102)
    ax2.hist(swCleanData["pGAPmid"][0:maxMatches], bins=percentbins, color='k', alpha=1 / nDataSet, label='Original SW %Gap')
    iter = 0
    for swRun in swUsData:
        iter += 1
        ax2.hist(swRun["pGAPmid"][0:maxMatches], bins=percentbins, alpha=1 / nDataSet, label='Altered SW' + str(iter) + ' %Gap')
    ax2.legend()
    ax2.set_xlabel("% Gap")
    ax2.set_ylabel("Frequency")
    ax2.set_xlim(0, 102)
    plt.show()

    # Histogram for ntGAP and nqGAP
    fig, axs = plt.subplots(2)
    fig.suptitle("Number of Gaps in either miRNA or 3'UTR sequence in miRNA 12-17 Region")
    ax1 = axs[0]
    ax2 = axs[1]
    ax1.hist(swCleanData["ntGAPmid"][0:maxMatches], bins=nGAPbins, color='k', alpha=1 / nDataSet,
             label="Original SW 3'UTR Gaps")
    iter = 0
    for swRun in swUsData:
        iter += 1
        ax1.hist(swRun["ntGAPmid"][0:maxMatches], bins=nGAPbins, alpha=1 / nDataSet,
                 label="Altered SW" + str(iter) + " 3'UTR Gaps")
    ax1.legend()
    ax1.set_xlabel("Number of 3'UTR Gaps")
    ax1.set_ylabel("Frequency")
    ax1.set_xlim(0, 25)
    ax2.hist(swCleanData["nqGAPmid"][0:maxMatches], bins=nGAPbins, color='k', alpha=1 / nDataSet,
             label="Original SW miRNA Gaps")
    iter = 0
    for swRun in swUsData:
        iter += 1
        ax2.hist(swRun["nqGAPmid"][0:maxMatches], bins=nGAPbins, alpha=1 / nDataSet,
                 label='Altered SW' + str(iter) + ' miRNA Gaps')
    ax2.legend()
    ax2.set_xlabel("Number of miRNA Gaps")
    ax2.set_ylabel("Frequency")
    ax2.set_xlim(0, 25)
    plt.show()

################################## START AND STOP POSITION PLOTS ##################################
    #Scatter plots of start and stop positions
    plt.scatter(swCleanData["start_col"][0:maxMatches], swCleanData["start_row"][0:maxMatches], alpha=.1,label="Original SW Start",c="blue")
    plt.scatter(swCleanData["end_col"][0:maxMatches], swCleanData["end_row"][0:maxMatches], alpha=.1, label="Original SW End", c="k")
    iter = 0
    for swRun in swUsData:
        iter += 1
        plt.scatter(swRun["start_col"][0:maxMatches], swRun["start_row"][0:maxMatches],alpha=.1,label="Adapted SW"+str(iter)+" Start")
        plt.scatter(swRun["end_col"][0:maxMatches], swRun["end_row"][0:maxMatches], alpha=.1, label="Adapted SW"+str(iter)+" End")
    plt.legend()
    plt.xlabel("3' UTR Position (bp)")
    plt.ylabel("miRNA Position (bp)")
    plt.xlim(0, 4000)
    plt.ylim(0, 25)
    plt.show()
    iter = nDataSet
    for swRun in swUsData[::-1]:
        iter -= 1
        plt.scatter(swRun["start_col"][0:maxMatches], swRun["start_row"][0:maxMatches], alpha=.1, label="Adapted SW"+str(iter)+" Start")
        plt.scatter(swRun["end_col"][0:maxMatches], swRun["end_row"][0:maxMatches], alpha=.1, label="Adapted SW"+str(iter)+" End")
    plt.scatter(swCleanData["start_col"][0:maxMatches], swCleanData["start_row"][0:maxMatches], alpha=.1, label="Original SW Start", c="blue")
    plt.scatter(swCleanData["end_col"][0:maxMatches], swCleanData["end_row"][0:maxMatches], alpha=.1, label="Original SW End", c="k")
    plt.legend()
    plt.xlabel("3' UTR Position (bp)")
    plt.ylabel("miRNA Position (bp)")
    plt.xlim(0, 4000)
    plt.ylim(0, 25)
    plt.title("Same as last plot but with the stack order flipped")
    plt.show()

    #histograms for start and stop position
    UTRstartstopbins = [i for i in range(4000) if i%20==0]
    fig, axs = plt.subplots(2)
    fig.suptitle("Start and Stop Postions on 3'UTR")
    ax1 = axs[0]
    ax2 = axs[1]
    ax1.hist(swCleanData["start_col"][0:maxMatches], bins=UTRstartstopbins, color='k', alpha=1/nDataSet, label='Original SW Start')
    iter = 0
    for swRun in swUsData:
        iter += 1
        ax1.hist(swRun["start_col"][0:maxMatches], bins=UTRstartstopbins, alpha=1/nDataSet, label='Altered SW'+str(iter)+' Start')
    ax1.legend()
    ax1.set_xlabel("3' UTR Position (bp)")
    ax1.set_ylabel("Frequency")
    ax1.set_xlim(0, 4000)
    ax2.hist(swCleanData["end_col"][0:maxMatches], bins=UTRstartstopbins, color='k', alpha=1/nDataSet, label='Original SW End')
    iter = 0
    for swRun in swUsData:
        iter += 1
        ax2.hist(swRun["end_col"][0:maxMatches], bins=UTRstartstopbins, alpha=1/nDataSet, label='Altered SW'+str(iter)+' End')
    ax2.legend()
    ax2.set_xlabel("3' UTR Position (bp)")
    ax2.set_ylabel("Frequency")
    ax2.set_xlim(0, 4000)
    plt.show()

    miRNAstartstopbins = [i/2 for i in range(50)]
    fig, axs = plt.subplots(2)
    fig.suptitle("Start and Stop Postions on miRNA")
    ax1 = axs[0]
    ax2 = axs[1]
    ax1.hist(swCleanData["start_row"][0:maxMatches], bins=miRNAstartstopbins, color='k', alpha=1/nDataSet, label='Original SW Start')
    iter = 0
    for swRun in swUsData:
        iter += 1
        ax1.hist(swRun["start_row"][0:maxMatches], bins=miRNAstartstopbins, alpha=1/nDataSet, label='Altered SW'+str(iter)+' Start')
    ax1.legend()
    ax1.set_xlabel("miRNA Position (bp)")
    ax1.set_ylabel("Frequency")
    ax1.set_xlim(0, 25)
    ax2.hist(swCleanData["end_row"][0:maxMatches], bins=miRNAstartstopbins, color='k', alpha=1/nDataSet, label='Original SW End')
    iter = 0
    for swRun in swUsData:
        iter += 1
        ax2.hist(swRun["end_row"][0:maxMatches], bins=miRNAstartstopbins, alpha=1/nDataSet, label='Altered SW'+str(iter)+' End')
    ax2.legend()
    ax2.set_xlabel("miRNA Position (bp)")
    ax2.set_ylabel("Frequency")
    ax2.set_xlim(0, 25)
    plt.show()

################################## LENGTH OF MATCHES PLOT ##################################
    lenPlotbins = [i/2 for i in range(2*40)]
    fig, axs = plt.subplots(2)
    fig.suptitle("Length of the Matches in the miRNA and the 3'UTR")
    ax1 = axs[0]
    ax2 = axs[1]
    ax1.hist(swCleanData["matchLenQ"][0:maxMatches], bins=lenPlotbins, color='k', alpha=1 / nDataSet, label='Original SW miRNA')
    iter = 0
    for swRun in swUsData:
        iter += 1
        ax1.hist(swRun["matchLenQ"][0:maxMatches], bins=lenPlotbins, alpha=1 / nDataSet, label='Altered SW' + str(iter) + ' miRAN')
    ax1.legend()
    ax1.set_xlabel("Matched miRNA Length (bp)")
    ax1.set_ylabel("Frequency")
    ax1.set_xlim(0, 40)
    ax2.hist(swCleanData["matchLenT"][0:maxMatches], bins=lenPlotbins, color='k', alpha=1 / nDataSet, label='Original SW End')
    iter = 0
    for swRun in swUsData:
        iter += 1
        ax2.hist(swRun["matchLenT"][0:maxMatches], bins=lenPlotbins, alpha=1 / nDataSet, label='Altered SW' + str(iter) + ' End')
    ax2.legend()
    ax2.set_xlabel("Matched 3'UTR Length (bp)")
    ax2.set_ylabel("Frequency")
    ax2.set_xlim(0, 40)
    plt.show()


#######################################################################################################################
n = len(sys.argv[-1])
swUSlist = sys.argv[-1][1:n-1]
swUSlist = swUSlist.split(',')
validation(sys.argv[-2],swUSlist)