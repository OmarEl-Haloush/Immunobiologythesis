import sys
import time
import math
from GenHLA import *
from GenHLAdis import *

Tb=time.perf_counter()

def MinAverage(file, HLAtable):
    data=open(file, 'r')
    av=list()
    var=list()
    datalist=[0.04286,0.29898,1.2702349990,3,28262,0.86711,1.29556,1.84724,0.49138,1.441746]
    chancelist=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
    #print("c")
    while len(av)<9:
        #print("a")
        av.append(float(0))
        var.append(float(0))
    #print(av)    
    linecount=0
    corrlist=list()
    for HLA in HLAtable:
        corrlist.append(HLA[0])
    for HLA in HLAtable:
        corrlist.append(HLA[1])
    for HLA in HLAtable:
        corrlist.append(HLA[2])
    #for i in range(0,len(corrlist)):
        #datalist[i]=datalist[i]-corrlist[i]-1
    print("Datalist :\n", datalist)    
    for line in data:
        linecount+=1
        line=line.split()
        #print(line)
        for i in range(0,len(line)):
            line[i]=float(corrlist[i])-float(line[i])-1
        #break    
        #line[0]=(17.0511-float(line[0]))
        #line[1]=(31.708-float(line[1]))
        #line[2]=(18.7228-float(line[2]))
        #line[3]=(8.7102-float(line[3]))
        #line[4]=(13.1346-float(line[4]))
        #line[5]=(10.6858-float(line[5]))
        #line[6]=(12.8232-float(line[6]))
        #line[7]=(22.5538-float(line[7]))
       # print(line[8])
        #line[8]=(15.4564-float(line[8]))
        #print("b")
        #print(line)


        #print(len(line))
        for i in range(0,len(line)):
            
            #print(av[i])
            av[i]+=line[i]
            var[i]+=line[i]**2
            if (line[i])>datalist[i]:
                #print(line[i], datalist[i])
                chancelist[i]+=1.0
            #print(av[i])
        #if linecount==10:
            #break    
        #print(i)    
        #print(av[8])    
        #break    
    newlist=list()
    varlist=list()
    chanlist=list()
    print("Corrected average  :\n")
    for fl in av:
       new=float(fl)/linecount
       print("\t", new)
       #new=math.sqrt(new)
       #ewlist.append(new)
    #print("Corrected average  :\n") 
    print("Variance of averages: \n") 
    for va in var:
        new=float(va)/linecount
        #arlist.append(new)
        print("\t", new)
    print("chances of occurance real data :\n")
    for chance in chancelist:
        new=float(chance)/linecount
        #hanlist.append(new)
        print("\t", new)
        
       
def Nantload():
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
                        print("bier")
                        
                print(len(PatNAntlist))    
                if cont==False:
                    PatNAntlist.append(newPatHLA)
        newfile=open("Mutationloads"+criterium+".txt", 'w')
        newfile.write("Patient\tHLA-A\tHLA-B\tHLA-C")
        
        for pat in PatNAntlist:
            newfile.write("\n"+str(pat.ID)+"\t"+str(pat.NALA)+"\t"+str(pat.NALB)+"\t"+str(pat.NALC))
        newfile.close()
        return
    loadNAntdata("testneo2.txt", "response")
    loadNAntdata("testneo2.txt", "nonresponse")
    loadNAntdata("testneo2.txt", "long-survival")

def CalcAveragecurvefit(HLA):
    avls=float()
    avrs=float()
    avnr=float()
    table=list()
    for runcount in HLA:
        avls+=float(runcount[9])
        avrs+=float(runcount[26])
        avnr+=float(runcount[72])
    #print(avnr, avls, avrs)    
    table.append(avnr/(len(HLA)))
    table.append(avls/(len(HLA)))
    table.append(avrs/(len(HLA)))
    return table
    
print("start")
#MinAverage("rundishla.txt")
#Nantload()
#HLAa=Curvefit("HLA-A", 10)
file=Generation(100000)
HLAtable=list()
HLAtable.append(CalcAveragecurvefit(Curvefit("HLA-A",1000000)))
HLAtable.append(CalcAveragecurvefit(Curvefit("HLA-B",1000000)))
HLAtable.append(CalcAveragecurvefit(Curvefit("HLA-C",1000000)))
#print(HLAtable)
MinAverage(file, HLAtable)



#print(len(HLAa))

print("Done")
print("Total Computation Time :", time.perf_counter()-Tb) 


input()