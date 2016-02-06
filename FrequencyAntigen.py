import sys
import time
#from Epitope import *
from GenHLA import *

Tb=time.perf_counter()

def loaddataFA(file, option,start=0, stop=8):
    datafile=file+".txt"
    data=open(datafile, 'r')
    patlist=list()
    newlist=list()
    for line in data:
        line=line.split()
        newpat=PatHLA(line)
        
        mer=str()
        i=0
        origin=line[1]
        if option=="Prop":
            line[1]=newpat.CreatingProperties(line[1])
        if line[0] in patlist:
            for pat in newlist:
                if pat.ID==line[0]:
                    for char in line[1]:
                        if i>=start and i<=stop:
                            mer+=str(char)
                            
                        i+=1
                    pat.Nean.append(mer)
                    newlist2=list()
                    newlist2.append(origin)
                    pat.Prop.append(newlist2)
                    pat.Matchlist.append(0)
                    break
        else:
            for char in line[1]:
                if i>1 and i<8:
                    mer+=str(char)
                i+=1
            newpat.Nean.append(mer)
            newlist2=list()
            newlist2.append(origin)
            newpat.Prop.append(newlist2)
            newpat.Matchlist.append(0)
            newlist.append(newpat)
            patlist.append(newpat.ID)
        #break    
    return newlist    
    
def loaddataPG(file):
    datafile=file+".txt"
    data=open(datafile, 'r')
    patlist=list()
    newlist=list()
    for line in data:
        line=line.split()
        newpat=PatHLA(line)
        newpat.Group=line[1]
        newlist.append(newpat)
    return newlist

def RecMatch(pa, minimumrecurrency=1):
    newlist=list()
    merlist=list()
    originlist=list()
    count=0
    for i in range(0,len(pa)-1):
        #if count>10:
        #    break
        for k in range(0,len(pa[i].Nean)-1):
            if pa[i].Nean[k] in merlist:
                #print("double 6mer recurrent!")
                continue
            for j in range(i+1, len(pa)-1):
                for l in range(0,len(pa[j].Nean)-1):
                    if pa[i].Nean[k]==pa[j].Nean[l]:
                        pa[i].Matchlist[k]+=1
                        count+=1
                        merlist.append(pa[i].Nean[k])
                        pa[i].Prop[k].append(pa[j].Prop[l])
                        #print(i, j, "'Match' found on:  ", pa[i].ID, pa[i].Nean[k], pa[j].ID, pa[j].Nean[l], pa[i].Prop)
    print(count, " Matches on 6 mer discovered over :  ", len(merlist),"  antigens of different patients")  
    for pat in pa:
        for i in range(0,len(pat.Matchlist)-1):
            if pat.Matchlist[i]>=minimumrecurrency:
                newlist.append(pat.Nean[i])
                originlist.append(pat.Prop[i])
                #print("usable match")
    print("Match discovered on : ", len(newlist), "  patients")       
    #return newlist, originlist
    return pa

def Freqcount(recurrent, pg, pa, ID, origin):
    freq=list()
    freq.append(list())
    freq[0].append("6mer")
    #freq.append(list())
    freq[0].append("Non response")
    #freq.append(list())
    freq[0].append("Long-Survival")
    #freq.append(list())
    freq[0].append("Response")
    newfile=open("freqcount"+str(ID)+".txt", 'w')
    
    for i in range(0,len(recurrent)):
        
        #print(rec)
        freq.append(list())
        freq[i+1].append(recurrent[i])
        freq[i+1].append(int(0))
        freq[i+1].append(int(0))
        freq[i+1].append(int(0))
        freq[i+1].append(origin[i])
        
    print(len(freq[0]), len(freq[1]))
    for i in range(0,len(recurrent)-1):
        #print(recurrent[i])
        for pat in pa:
            for antigen in pat.Nean:
                if recurrent[i]==antigen:
                    #print("matched recurrent antigen")
                    for group in pg:
                        if pat.ID==group.ID:
                            #print("found patient")
                            if group.Group=="nonresponse":
                                freq[i+1][1]=int(freq[i+1][1])+1
                            if group.Group=="long-survival":
                                freq[i+1][2]=int(freq[i+1][2])+1
                            if group.Group=="response":
                                freq[i+1][3]=int(freq[i+1][3])+1
    for line in freq:
        for arg in line:
            newfile.write(str(arg)+"\t")
        newfile.write("\n")    
    return freq
    

#pa=loaddataFA("DataMut1")
def MatchonAA(option, start=0, stop=8):
    pa=loaddataFA("DataMut1", option, start, stop)
    print(len(pa))
    #print(pg[0].Nean, pg[0].ID)
    pg=loaddataPG("DataGroups")
    print(len(pg))
    ID=1
    base=RecMatch(pa)
    while ID<12:
        Tr=time.perf_counter()
        #recurrent, origin=RecMatch(pa,ID)
        #print(recurrent[0],recurrent[1])
        #print(len(recurrent))
        #print(origin[0], origin[1])
        #print(len(origin))
        recurrent=list()
        origin=list()
        for pat in pa:
            for i in range(0,len(pat.Matchlist)-1):
                if pat.Matchlist[i]>=ID:
                    recurrent.append(pat.Nean[i])
                    origin.append(pat.Prop[i])
                    #print("usable match")
        print("Match discovered on : ", len(recurrent), "  antigens") 
        freq=Freqcount(recurrent, pg,pa, ID, origin)
        print("\t Next round starting :",ID, time.perf_counter()-Tr, " seconds")
        ID+=1
    print("Done")
    print("Total Computation Time :", time.perf_counter()-Tb)     
