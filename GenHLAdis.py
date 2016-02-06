import sys
import time
from GenHLA import *
import random

Tb=time.perf_counter()

def loadPatHLAdata(file):
    baseHLAlist=list()
    hlafile=open(file, 'r')
    for line in hlafile:
        newPatHLA=PatHLA(line.split())
        baseHLAlist.append(newPatHLA)
    hlafile.close()
    '''newfile=open("testgenhladis2.txt", 'w')
    for patHLA in baseHLAlist:
        newfile.write(patHLA.ID)
        for HLA in patHLA.HLAlist:
            newfile.write("\t"+HLA)
        newfile.write("\n")
    newfile.close() '''   
    return baseHLAlist

def genHLAGroups(baseHLAlist):
    def AssignGroup(PatHLA, HLAgroups, N):
        reschance=(27)/N
        nreschance=(73)/N
        lschance=(10)/N
        chance=random.random()
        #print(reschance, nreschance, lschance, chance)
        if chance<=lschance and (len(HLAgroups[1]))<=9:
            HLAgroups[1].append(PatHLA)
            return HLAgroups
        if chance<=(reschance+lschance) and (len(HLAgroups[2]))<=26:
            HLAgroups[2].append(PatHLA)
            return HLAgroups
        if chance<=(reschance+lschance+nreschance) and (len(HLAgroups[0]))<=72:
            HLAgroups[0].append(PatHLA)
            return HLAgroups
        else:
            #AssignGroup(PatHLA, HLAgroups, N)
            return HLAgroups        
    newHLAlist=list()
    for pathla in baseHLAlist:
        newHLAlist.append(pathla)
    HLAgroups=list()
    HLAgroups.append(list())
    HLAgroups.append(list())
    HLAgroups.append(list())
    while len(newHLAlist)>0:
        random.seed()
        number=random.randint(0,len(newHLAlist)-1)
        PatHLA=newHLAlist.pop(number)
        HLAgroups=AssignGroup(PatHLA, HLAgroups, len(newHLAlist)+1)   
    #print("length baselist",len(newHLAlist))
    return HLAgroups
    

    
def genHLADis(HLAgroups, newfile, count):
    newfile.write("\n")
    for HLAgroup in HLAgroups:
        CountA=list()
        CountB=list()
        CountC=list()
        for PatHLA in HLAgroup:
            for HLA in PatHLA.HLAlist:
                if "HLA-A" in HLA:
                    if HLA not in CountA:
                        CountA.append(HLA)
                if "HLA-B" in HLA:
                    if HLA not in CountB:
                        CountB.append(HLA)
                if "HLA-C" in HLA:
                    if HLA not in CountC:
                        CountC.append(HLA)
                  
        newfile.write(str(len(CountA)-1)+"\t"+str(len(CountB)-1)+"\t"+str(len(CountC)-1)+"\t")
        
    return #HLAdis    
    
def Generation(runcount):    
    baseHLAlist=loadPatHLAdata("testgenhladis.txt")
    returnfile="HLAdis.txt"
    newfile=open("HLAdis.txt",'w')
    #newfile.write("Non Response\t\t\tLong-Survival\t\t\tResponse\n")
    #newfile.write("\tHLA-A\tHLA-B\tHLA-C\tHLA-A\tHLA-B\tHLA-C\tHLA-A\tHLA-B\tHLA-C\n")
    count=0
    while count<runcount:
        count+=1
        HLAgroups=genHLAGroups(baseHLAlist)
        genHLADis(HLAgroups, newfile, count)
        #for HLAgroup in HLAgroups:
    
            #for PatHLA in HLAgroup:
            #    print(PatHLA.HLAlist)
            #print("nextgroup")    
        '''HLAdis=genHLADis(HLAgroups)'''
        print("generationround= ", count)
    
    newfile.close()
    return returnfile

def Curvefit(HLA,N):
    baseHLAlist=loadPatHLAdata("testgenhladis.txt")
    #print(len(baseHLAlist))
    #newfile=open("curvedataHLC-A10000.txt",'w')
    #newfile.write("Runcount")
    count=0
    def PatNumUniqueHLA(baseHLAlist, count, HLA):
        
        newHLAlist=list()
        #newfile.write("\n"+str(count))
        for pathla in baseHLAlist:
            newHLAlist.append(pathla)
        HLAgroups=list()
        runcount=list()
        while len(newHLAlist)>0:
            random.seed()
            number=random.randint(0,len(newHLAlist)-1)
            PatHLA=newHLAlist.pop(number)
            for HLAs in PatHLA.HLAlist:
                if HLA in HLAs:
                    
                    if HLAs in HLAgroups:
                        #print(HLAs)    
                        continue
                        
                    else:
                        HLAgroups.append(HLAs)
            runcount.append(len(HLAgroups))
        return runcount
            #newfile.write("\t"+str(len(HLAgroups)))  
    curvefit=list()        
    while count<N:
        count+=1
        #print(count)
        curvefit.append(PatNumUniqueHLA(baseHLAlist,count,HLA))
    #newfile.close()
    #print(len(curvefit[0]))
    #print(curvefit[0])
    #print(len(curvefit))
    return curvefit

def Nantload(option="false"):
    def loadNAntdata(file, criterium):
        data=open(file, 'r')
        PatNAntlist=list()
        for line in data:
            cont=False
            newline=line.split()
            #print(newline)
            newPatHLA=PatHLA(newline)
            #print(newline[2], criterium)
            #break
            if newline[2]==criterium:
                if "HLA-A" in newline[1]:
                    newPatHLA.NALA+=1
                if "HLA-B" in newline[1]:
                    newPatHLA.NALB+=1
                if "HLA-C" in newline[1]:
                    newPatHLA.NALC+=1
                for pat in PatNAntlist:
                    if pat.ID==newPatHLA.ID:
                        print(pat.ID, newPatHLA.ID)
                        pat.NALA+=newPatHLA.NALA
                        pat.NALB+=newPatHLA.NALB
                        pat.NALC+=newPatHLA.NALC
                        cont=True
                        
                        
                print(len(PatNAntlist))    
                if cont==False:
                    PatNAntlist.append(newPatHLA)
        newfile=open(Option+"loads"+criterium+".txt", 'w')
        newfile.write("Patient\tHLA-A\tHLA-B\tHLA-C")
        
        for pat in PatNAntlist:
            newfile.write("\n"+str(pat.ID)+"\t"+str(pat.NALA)+"\t"+str(pat.NALB)+"\t"+str(pat.NALC))
        newfile.close()
        return
    if option=="Neo":
        loadNAntdata("testneo.txt", "response")
        loadNAntdata("testneo.txt", "nonresponse")
        loadNAntdata("testneo.txt", "long-survival")
    if option=="Mut":
        loadNAntdata("testneo2.txt", "response")
        loadNAntdata("testneo2.txt", "nonresponse")
        loadNAntdata("testneo2.txt", "long-survival")
    else:
        print("Option unknown, please use the option Neo or Mut")

def Equal6mer():
    def loadNAdata(file, mincount=6):
        data=open(file, 'r')
        merlist=list()
        
        #neoantigencount=0
        for line in data:

            newline=line.split()
            #print(line, newline)
            anti=PatHLA(newline)
            
            mer=list()
            i=0
            for char in anti.ID:
                if i>1 and i<8:
                    mer.append(char)
                i+=1
            anti.Nean=mer 
            #print(anti.Nean)
            merlist.append(anti)
            
            #if neoantigencount==2:
            #    break  
            #neoantigencount+=1
        i=0
        j=0
        mnlist=list()
        molist=list()
        for i in range(0,len(merlist)-1):
            if merlist[i].Nean in mnlist:
                continue
            for j in range(i+1,len(merlist)-1):
                aacount=0
                #print(merlist[i].Nean, merlist[j].Nean)
                for k in range(0,6):
                    #print(k)
                    #print(merlist[i].Nean[k], merlist[j].Nean[k])
                    
                    if merlist[i].Nean[k]==merlist[j].Nean[k]:
                        aacount+=1
                if aacount==mincount:
                    merlist[i].matchcount+=1
                    #print("bier")
            if merlist[i].matchcount>0:
                mnlist.append(merlist[i].Nean)
                molist.append(merlist[i])
        newfile=open("6mermatch.txt", 'w')
        newfile.write("Neoantigen\tMatchcount")
        for mer in molist:
            newfile.write("\n"+mer.ID+"\t"+str(mer.matchcount))
        newfile.close()
    loadNAdata("test6mer.txt")

def MatchOcc():
    def loadCompareData(matchfile, neofile, runcount):
        def  mercreate(line, merlist):
            newline=line.split()
            #print(line, newline)
            anti=PatHLA(newline)
            mer=list()
            i=0
            for char in anti.ID:
                if i>1 and i<8:
                    mer.append(char)
                i+=1
            anti.Nean=mer
            merlist.append(anti)
            return merlist
        match=open(matchfile, 'r')
        neo=open(neofile, 'r')
        newfile=open("uniqueNeos"+str(runcount)+".text", 'w')
        matchlist=list()
        neolist=list()
        for line in match:
            matchlist=mercreate(line,matchlist)
            #break
        for line in neo:
            neolist=mercreate(line,neolist)
            #break
        #print(len(matchlist), len(neolist))
        #print(matchlist[0].ID, neolist[0].ID)
        i=0
        j=0
        mnlist=list()
        molist=list()
        for i in range(0,len(matchlist)-1):
            if matchlist[i].Nean in mnlist:
                print("whisky")
                continue
            for j in range(0,len(neolist)-1):
                aacount=0
                #print(merlist[i].Nean, merlist[j].Nean)
                for k in range(0,6):
                    #print(k)
                    #print(merlist[i].Nean[k], merlist[j].Nean[k])
                    
                    if matchlist[i].Nean[k]==neolist[j].Nean[k]:
                        aacount+=1
                #print("wijn")        
                if aacount>=5:
                    matchlist[i].matchcount+=1
                    #print("bier")
            if matchlist[i].matchcount==0:
                #print("bier")
                mnlist.append(matchlist[i].Nean)
                molist.append(matchlist[i])    
        
        newfile.write("Neoantigen\tMatchcount")
        for mer in molist:
            newfile.write("\n"+mer.ID+"\t"+str(mer.matchcount))
        newfile.close()
        print("Number of unique recurrent signatures discovered :", len(molist))
        return
    runcount=1
    loadCompareData("matchedneosRes.txt", "neosNonres.txt", runcount)
    runcount+=1
    loadCompareData("matchedneosRes.txt", "neosLS.txt", runcount)
    runcount+=1
    loadCompareData("matchedneosNonres.txt", "neosRes.txt", runcount)
    runcount+=1
    loadCompareData("matchedneosNonres.txt", "neosLS.txt", runcount)
    runcount+=1
    loadCompareData("matchedneosLS.txt", "neosRes.txt", runcount)
    runcount+=1
    loadCompareData("matchedneosLS.txt", "neosNonres.txt", runcount)
    return

def Equalpropmer():
    def loadpropdata(file, mincount=6):
        data=open(file, 'r')
        merlist=list()
        
        #neoantigencount=0
        for line in data:

            newline=line.split()
            #print(line, newline)
            anti=PatHLA(newline)
            
            mer=list()
            i=0
            for char in anti.ID:
                if i>1 and i<8:
                    mer.append(char)
                    
                i+=1
            anti.Nean=mer 
            anti.Prop=anti.CreatingProperties(anti.Nean)
            #print(anti.Nean)
            merlist.append(anti)
            
            #if neoantigencount==2:
            #    break  
            #neoantigencount+=1
        i=0
        j=0
        mnlist=list()
        molist=list()
        for i in range(0,len(merlist)-1):
            if merlist[i].Prop in mnlist:
                continue
            for j in range(i+1,len(merlist)-1):
                aacount=0
                #print(merlist[i].Nean, merlist[j].Nean)
                for k in range(0,6):
                    #print(k)
                    #print(merlist[i].Nean[k], merlist[j].Nean[k])
                    
                    if merlist[i].Nean[k]==merlist[j].Nean[k]:
                        aacount+=1
                if aacount==5:
                    propcount=0
                    for k in range(0,6):
                        #print(k)
                        #print(merlist[i].Nean[k], merlist[j].Nean[k])
                        #print(k)
                        if merlist[i].Prop[k]==merlist[j].Prop[k]:
                            propcount+=1
                    if propcount==6:
                        merlist[i].matchcount+=1
                    #print("bier")
            if merlist[i].matchcount>0:
                mnlist.append(merlist[i].Prop)
                molist.append(merlist[i])
        newfile=open("6propmermatch.txt", 'w')
        newfile.write("Neoantigen\tProperty\tMatchcount")
        print("Number of recurrent antigens discovered :",len(molist))
        for mer in molist:
            prop=str()
            for char in mer.Prop:
                prop+=str(char)
            newfile.write("\n"+mer.ID+"\t"+prop+"\t"+str(mer.matchcount))
        newfile.close()
    loadpropdata("testpropmer.txt")
    
def MatchPropOcc():
    def loadCompareData(matchfile, neofile, runcount):
        def  mercreate(line, merlist):
            newline=line.split()
            #print(line, newline)
            anti=PatHLA(newline)
            mer=list()
            i=0
            for char in anti.ID:
                if i>1 and i<8:
                    mer.append(char)
                i+=1
            anti.Nean=mer
            anti.Prop=anti.CreatingProperties(anti.Nean)
            merlist.append(anti)
            return merlist
        match=open(matchfile, 'r')
        neo=open(neofile, 'r')
        newfile=open("uniqueNeos"+str(runcount)+".text", 'w')
        matchlist=list()
        neolist=list()
        for line in match:
            matchlist=mercreate(line,matchlist)
            
            #break
        for line in neo:
            neolist=mercreate(line,neolist)
            #break
        #print(len(matchlist), len(neolist))
        #print(matchlist[0].ID, neolist[0].ID)
        i=0
        j=0
        mnlist=list()
        molist=list()
        for i in range(0,len(matchlist)-1):
            if matchlist[i].Prop in mnlist:
                print("whisky")
                continue
            for j in range(0,len(neolist)-1):
                aacount=0
                #print(merlist[i].Nean, merlist[j].Nean)
                for k in range(0,6):
                    #print(k)
                    #print(merlist[i].Nean[k], merlist[j].Nean[k])
                    
                    if matchlist[i].Prop[k]==neolist[j].Prop[k]:
                        aacount+=1
                #print("wijn")        
                if aacount>=5:
                    matchlist[i].matchcount+=1
                    #print("bier")
            if matchlist[i].matchcount==0:
                #print("bier")
                mnlist.append(matchlist[i].Prop)
                molist.append(matchlist[i])    
        
        newfile.write("Neoantigen\tProperty\tMatchcount")
        for mer in molist:
            prop=str()
            for char in mer.Prop:
                prop+=str(char)
            newfile.write("\n"+mer.ID+"\t"+ prop+"\t"+str(mer.matchcount))
        newfile.close()
        print("Number of unique recurrent signatures discovered :", len(molist))
        return
    runcount=1
    loadCompareData("matchedneosRes.txt", "neosNonres.txt", runcount)
    runcount+=1
    loadCompareData("matchedneosRes.txt", "neosLS.txt", runcount)
    runcount+=1
    loadCompareData("matchedneosNonres.txt", "neosRes.txt", runcount)
    runcount+=1
    loadCompareData("matchedneosNonres.txt", "neosLS.txt", runcount)
    runcount+=1
    loadCompareData("matchedneosLS.txt", "neosRes.txt", runcount)
    runcount+=1
    loadCompareData("matchedneosLS.txt", "neosNonres.txt", runcount)
    return    
    
#Generation()            
#Equal6mer()
#MatchOcc()
#Equalpropmer()
#MatchPropOcc()
#print("Starting HLA-A")
#Curvefit("HLA-A",10)
#print("Starting HLA-B")
#Curvefit("HLA-B10000")
#print("Starting HLA-C")
#Curvefit()
#print("Done")
#print("Total Computation Time :", time.perf_counter()-Tb) 


#input()