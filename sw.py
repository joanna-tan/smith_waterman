#!/usr/bin/env python3
# Smith-Waterman: dynamic-programming-based matching. -*- mode: python -*-
# (c) 2011-2016 duane a. bailey and <your name here>

#python3 lab4/fastad.py < lab9/utr.fa | python3 lab9/sw.py `python3 lab4/fastad.py< lab9/hsa-let-7a-1.fa`

# first character in match counts as a "match," remaining are "similar"
# everything else is a mismatch
rnaMatch = {'A': 'U', 'C': 'G', 'G': 'CU', 'U': 'AG'}
    #notice the similar matches for G and U (second letter)
        #Allows for a range of "wobble" matches


# these parameters determine how to score match
# if seeding is used, the scores can be multiplied in the seed region
rnaValues = {"match": 1.0,
             "targetGap": -1, "queryGap": -.75,        #t-1,q-.75 slightly preferred query gaps as query gaps still allow for the miRNA to have a match at every location
             "seedqGap" : .25,                         #.25
             "mismatch": -1.0,
             "similar": 1.0,
             "seedLength": 7.0,
             "seed6merWeight": 2.0,       #2 Want perfect matches of this 6mer
             "seed7merm8Weight": 0.6,     #.6 Want this 7mer to be prefered to the 7mer-A1
             "seed7merA1Weight": 0.5,     #.5 See above point
             "seed8merWeight": 0.75,      #.75 Want to be same scoring (maybe a bit higher) as 7mer-m8
             "suppression": -1.5}         #-1.5 For suppression of similar matches in the 12-17 range

def sw(query, target, scoreCount=30,
       match=rnaMatch, scoreValue=rnaValues,
       qcomment="<<< query", tcomment=">>> target",
       csv=True):
    """Smith-Waterman matcher.

    Prints the scoreCount best matches based on the settings of the
    matching function and the parameters that control the scoring.
    """

    # table has an extra row and column
    height = len(query) + 1
    width = len(target) + 1

    # initialize the edges of the table.
    # each table entry has a score, a match type, and a
    table = {}
    for row in range(0, height):
        table[row, 0] = {"score": 0.0, "start": (row, 0), "type": "begin", "merBoost": 0}
    for col in range(1, width):
        table[0, col] = {"score": 0.0, "start": (0, col), "type": "begin", "merBoost": 0}

    # grab parameters for later use
    seedLength = scoreValue["seedLength"]
    seed6merWeight = scoreValue["seed6merWeight"]
    boost = 0
    seed7merA1Weight = scoreValue["seed7merA1Weight"]
    seed7merm8Weight = scoreValue["seed7merm8Weight"]
    seed8merWeight = scoreValue["seed8merWeight"]
    suppression = scoreValue["suppression"]
    seedqGap = scoreValue["seedqGap"]

    # fill out the Smith Waterman matrix with values
    for row in range(1, height):
        for col in range(1, width):

            # grab the query and target characters
            # strings start at 0, but bases count from 1
            qc = query[row - 1]
            tc = target[col - 1]

            # assume a match or similar extension
            newScore = {"score": table[row - 1, col - 1]["score"],
                        "start": table[row - 1, col - 1]["start"],
                        "merBoost": 0}

            # "row" is the position of the base currently being matched
            # get rid of one of these two terms if you want
            # seed region at one end
            in_seed = (row <= seedLength) #or ((height - row) <= seedLength)
                # We got rid of ((height - row) <= seedLength) because the seeds occur at the begging of the miRNA
            in_12to17 = row >= 12 and row <= 17
                # Used to boost Watson-Crick pairing in this region

            # types of matching: strong/match, weak/similar, bad/mismatch
            if match[qc][0] == tc:                      #match[] = allows us to check the compliment
                newScore["type"] = "match"
            elif (tc in match[qc]):
                newScore["type"] = "similar"
            else:
                newScore["type"] = "mismatch"

            if row == 1 and len(target) - col >= 8:                                # If at first position of miRNA
                if tc == "A":                           # Checks if there is an A at this position in the 3'UTR
                    newScore["merBoost"] = seed7merA1Weight
                    if match[query[row + 7 - 1]] == target[col + 7 - 1]:    # Checks if there is also a match at position 8 of the miRNA
                        newScore["merBoost"] = seed8merWeight
                    # Add score in for position 1 with boost from being in 7mer-A1 or 8mer
                    newScore["score"] += (1 + newScore["merBoost"]) * scoreValue[newScore["type"]]
                elif match[query[row + 7 - 1]][0] == target[col + 7 - 1]:      # Checks if there is just a match at miRNA position 8
                    newScore["merBoost"] = seed7merm8Weight
                    # Add score in for position 1 with no boost
                    newScore["score"] += scoreValue[newScore["type"]]
                else:                                                       # If just a 6mer seed
                    newScore["merBoost"] = 0
                    newScore["score"] += scoreValue[newScore["type"]]
            else:                                       # For postion 2+ of the miRNA (or if less than an 8mer left in UTR)
                if in_seed:
                    if "merBoost" in table[row-1,col-1]:
                        newScore["merBoost"] = table[row-1,col-1]["merBoost"]
                        boost = newScore["merBoost"]
                    newScore["score"] += (seed6merWeight + boost) * scoreValue[newScore["type"]]
                        # scoreValue gets scaled up by value determined by which "mer" it is part of
                elif row == 8:
                    if len(target) - col >= 8:
                        if table[row-1,col-1]["merBoost"] == seed7merm8Weight or table[row-1,col-1]["merBoost"] == seed8merWeight:
                            boost = table[row-1,col-1]["merBoost"]
                            # If in 8mer or 7mer-m8 and with more than 8mer left in UTR, scores scaled up accordingly
                            newScore["score"] += (seed6merWeight + boost) * scoreValue[newScore["type"]]
                        else:
                            newScore["score"] += scoreValue[newScore["type"]]
                    else:
                        newScore["score"] += scoreValue[newScore["type"]]

                elif in_12to17:                                 # Boost sites with Watson-Crick pairing in 12-17
                    supVal = 0
                    if newScore["type"] == "similar":
                        supVal = suppression                         #changes similars to mismatches scoring-wise
                    newScore["score"] += scoreValue[newScore["type"]] + supVal
                else:
                    newScore["score"] += scoreValue[newScore["type"]]

            # two other ways these strings may be aligned:
            # current query character is dropped (a targetGap)
            # or the current target character is dropped (a queryGap)
            tgboost = 0
            # Score target gaps more harshly in the seed since that means no base paring for the miRNA
            if in_seed:
                if "merBoost" in newScore:
                    tgboost = newScore["merBoost"]
                tgap_score = table[row - 1, col]["score"] + ((seed6merWeight+tgboost)*scoreValue["targetGap"])
            elif row == 8:
                if table[row - 1, col]["merBoost"] == seed7merm8Weight or table[row - 1, col]["merBoost"] == seed8merWeight:
                    tgboost = table[row - 1, col]["merBoost"]
                    tgap_score = table[row - 1, col]["score"] + ((seed6merWeight+tgboost)*scoreValue["targetGap"])
                else:
                    tgap_score = table[row - 1, col]["score"] + scoreValue["targetGap"]
            else:
                tgap_score = table[row - 1, col]["score"] + scoreValue["targetGap"]

            qgboost = 0
            # Score query gaps less harshly in the seed since they may allow for better base pairing for the miRNA
            if in_seed:
                if "merBoost" in newScore:
                    qgboost = newScore["merBoost"]
                qgap_score = table[row, col - 1]["score"] + scoreValue["queryGap"] + seedqGap
            elif row == 8:
                if table[row, col - 1]["merBoost"] == seed7merm8Weight or table[row, col - 1]["merBoost"] == seed8merWeight:
                    tgboost = table[row, col - 1]["merBoost"]
                    qgap_score = table[row, col - 1]["score"] + scoreValue["queryGap"] + seedqGap
                else:
                    qgap_score = table[row, col - 1]["score"] + scoreValue["queryGap"]
            else:
                qgap_score = table[row, col - 1]["score"] + scoreValue["queryGap"]

            # check if dropping target is a better approach
            if tgap_score > newScore["score"]:
                newScore = {"score": tgap_score,
                            "start": table[row - 1, col]["start"],
                            "merBoost": tgboost,
                            "type": "targetGap"}

            # perhaps dropping a query character is better
            if qgap_score > newScore["score"]:
                newScore = {"score": qgap_score,
                            "start": table[row, col - 1]["start"],
                            "merBoost": qgboost,
                            "type": "queryGap"}

            # if the best score is not positive,
            # we terminate any matching here
            if newScore["score"] <= 0:
                newScore = {"score": 0,
                            "start": (row, col),
                            "merBoost": 0,
                            "type": "begin"}

            # write the entry in the table
            table[row, col] = newScore

    # IF you have a score adjustment pass, put that here, making adjustments
    # to table[row,col]["score"].

    # Finally, compute the top scoreCount highest scores from the table
    # if there are ties, we may get a few more
    scores = [(table[row, col]["score"], (row, col))
              for row in range(1, height) for col in range(1, width)]
    # sorts based on (1) score, (2) position in target
    scores = sorted(scores, key=lambda s: (s[0], width - s[1][1]), reverse=True)

    # print highest scores, plus low-end ties if needed.
    # I changed this because the seed boosting part of the algorithm meant that the same 1-2 starting points with the
    # same final score were returned by the old high score method while this method attempts to show different/more
    # starting points and final scores.  The information you lose is alternative matchings for the last 14 or so
    # nucleotides given the same starting seed.  This method also tends to highlight (show more often) matches with more
    # gaps (to get the traceback start position sufficiently different) but that does not mean the rest of the algorithm
    # is actually finding more of these matches.  I tried upping the gap penalties a bit to counteract this but it had
    # limited effectiveness.  I also think that gaps in the miRNA are more favorable as they can still allow for perfect
    # base pairing of the miRNA just with some bulges in the UTR so I have changed the gap panalties to reflect this idea.
    highScores = []
    lastPos = (0,0)
    lastScore = 0
    numScores = 0
    lowScore = 100
    if len(scores) <= scoreCount:
        highScores = scores
    else:
        for ind,score in enumerate(scores):
            if abs(score[1][0]-lastPos[0]) > 8 or abs(score[1][1]-lastPos[1]) > 8 or score[0] != lastScore:
                highScores += [score]
                lastPos = score[1]
                lastScore = score[0]
                numScores += 1
                if numScores == scoreCount:
                    lowScore = lastScore
                if scores[ind+1][0] != lowScore and numScores >= scoreCount:        # First bit allows for low-end ties
                    if scores[ind+1][0] != lowScore:                    # Keeps the code from printing 1 extra match
                        del highScores[-1]
                    break


    # print the high scores
    printScores(highScores, table, query, qcomment, target, tcomment, csv=csv)


def printScores(high_scores, table, query, qcomment, target, tcomment, csv=True):
    """Print matches, possibly in comma-separated values format."""
    if csv:
        fout = open("newestSW.txt", 'w')

    scoreNumber = 1
    for entry in high_scores:
        # get the score
        (s, end) = entry
        # the location (in table) of end of match
        (row, col) = end
        (end_row, end_col) = end
        # the beginning of the match
        (start_row, start_col) = table[end]["start"]

        # these strings grow to form picture of alignment
        t_str = ""
        align_str = ""
        q_str = ""

        # Counters used to determine %identity, %similarity, etc.
        seedCounter = 0     # need seedCounter because you might end up in the seed for longer than the seed length (8 for us) due to gaps
        midCounter = 0      # need seedCounter because you might end up in the 12-17 postions for longer than 5 due to gaps
        nIDtot = 0
        nSIMtot = 0
        nBOTHtot = 0
        nGAPtot = 0
        ntGAPtot = 0
        nqGAPtot = 0
        nIDseed = 0
        nSIMseed = 0
        nBOTHseed = 0
        nGAPseed = 0
        ntGAPseed = 0
        nqGAPseed = 0
        nIDmid = 0
        nSIMmid = 0
        nBOTHmid = 0
        nGAPmid = 0
        ntGAPmid = 0
        nqGAPmid = 0

        # traverse the table from end to beginning of match,
        # building a picture of the query, target, and their
        # relative alignments
        while table[row, col]["type"] != "begin":
            type = table[row, col]["type"]
            # based on the type of match, we introduce a match,
            # mismatch, or gap.  We assemble (from right to left)
            # the picture of the alignment
            if type == "match":
                nIDtot += 1
                nBOTHtot += 1
                if row <= 8:
                    nIDseed += 1
                    nBOTHseed += 1
                    seedCounter += 1
                elif 12 <= row <= 17:
                    nIDmid += 1
                    nBOTHmid += 1
                    midCounter += 1
                (row, col) = (row - 1, col - 1)
                t_str = target[col] + t_str
                align_str = "|" + align_str
                q_str = query[row] + q_str
            elif type == "similar":
                nSIMtot += 1
                nBOTHtot += 1
                if row <= 8:
                    nSIMseed += 1
                    nBOTHseed += 1
                    seedCounter += 1
                elif 12 <= row <= 17:
                    nSIMmid += 1
                    nBOTHmid += 1
                    midCounter += 1
                (row, col) = (row - 1, col - 1)
                t_str = target[col] + t_str
                align_str = ":" + align_str
                q_str = query[row] + q_str
            elif type == "mismatch":
                if row <= 8:
                    seedCounter += 1
                elif 12 <= row <= 17:
                    midCounter += 1
                (row, col) = (row - 1, col - 1)
                t_str = target[col] + t_str
                align_str = " " + align_str
                q_str = query[row] + q_str
            elif type == "queryGap":
                nGAPtot += 1
                nqGAPtot += 1
                if row <= 8:
                    nGAPseed += 1
                    nqGAPseed += 1
                    seedCounter += 1
                elif 12 <= row <= 17:
                    nGAPmid += 1
                    nqGAPmid += 1
                    midCounter += 1
                (row, col) = (row, col - 1)
                t_str = target[col] + t_str
                align_str = " " + align_str
                q_str = "-" + q_str  # the gap
            elif type == "targetGap":
                nGAPtot += 1
                ntGAPtot += 1
                if row <= 8:
                    nGAPseed += 1
                    ntGAPseed += 1
                    seedCounter += 1
                elif 12 <= row <= 17:
                    nGAPmid += 1
                    ntGAPmid += 1
                    midCounter += 1
                (row, col) = (row - 1, col)
                t_str = "-" + t_str  # the gap
                align_str = " " + align_str
                q_str = query[row] + q_str

        pIDtot = 100*(nIDtot/len(align_str))
        pSIMtot = 100*(nSIMtot/len(align_str))
        pBOTHtot = 100*(nBOTHtot/len(align_str))
        pGAPtot = 100*(nGAPtot/len(align_str))
        if seedCounter != 0:
            pIDseed = 100 * (nIDseed / seedCounter)
            pSIMseed = 100 * (nSIMseed / seedCounter)
            pBOTHseed = 100 * (nBOTHseed / seedCounter)
            pGAPseed = 100 * (nGAPseed / seedCounter)
        else:
            pIDseed = 0
            pSIMseed = 0
            pBOTHseed = 0
            pGAPseed = 0
        if midCounter != 0:
            pIDmid = 100 * (nIDmid / midCounter)
            pSIMmid = 100 * (nSIMmid / midCounter)
            pBOTHmid = 100 * (nBOTHmid / midCounter)
            pGAPmid = 100 * (nGAPmid / midCounter)
        else:
            pIDmid = 0
            pSIMmid = 0
            pBOTHmid = 0
            pGAPmid = 0

        # actually print the information
        if csv:
            # for use in Excel and Numbers
            fout.write('{},{},"{}","{}",{},"{}",{},{},{:.2f},{:.2f},{:.2f},{:.2f},{:.2f},{:.2f},'
                       '{:.2f},{:.2f},{:.2f},{:.2f},{:.2f},{:.2f},{},{},{},{},{},{}\n'.format(start_col + 1,
                                                    table[end]["score"], t_str,
                                                    align_str,
                                                    start_row + 1, q_str, end_col, end_row,
                                                    pIDtot,
                                                    pSIMtot,
                                                    pBOTHtot,
                                                    pGAPtot,
                                                    pIDseed,
                                                    pSIMseed,
                                                    pBOTHseed,
                                                    pGAPseed,
                                                    pIDmid,
                                                    pSIMmid,
                                                    pBOTHmid,
                                                    pGAPmid,
                                                    ntGAPtot,
                                                    nqGAPtot,
                                                    ntGAPseed,
                                                    nqGAPseed,
                                                    ntGAPmid,
                                                    nqGAPmid))
        else:
            # for use by humans
            print("Match #{}, score={}".format(scoreNumber, table[end]["score"]))
            print("{:4d}: {} :{:<4d} {}".format(start_col + 1, t_str, end_col,
                                                tcomment))
            print("      {}".format(align_str))
            print("{:4d}: {} :{:<4d} {}\n".format(start_row + 1, q_str, end_row,
                                                  qcomment))
        scoreNumber = scoreNumber + 1
    if csv:
        fout.close()

if __name__ == "__main__":
    import sys
    q = sys.argv[1]  # first argument is the query
    for t in sys.stdin.readlines():
        # preprocess these in whatever way seems useful
        # this is a UTR region
        t = t.strip()  # remove leading and trailing whitespace
        t = t.upper()  # make upper case
        t = t.replace('T', 'U')  # ensure RNA alphabet

        # this is the miRNA
        q = q.upper()  # make upper case
        q = q.replace('T', 'U')  # ensure RNA alphabet
        q = q[::-1]    # reverse

        # perform the matching
        sw(q, t)
