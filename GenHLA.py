import sys

class PatHLA:
    
    def __init__(self, PatHLAline=list("Empty")):
        self.ID=str(PatHLAline[0])
        self.Group=str()
        self.NALA=0
        self.NALB=0
        self.NALC=0
        self.matchcount=0
        self.Nean=list()
        self.Prop=list()
        self.Matchlist=list()
        try:
            self.HLAlist=list(PatHLAline[1:])
        except:
            self.HLAlist=list()
           
        
        return 
    
    def CreatingProperties(self, Mutation):
        #Propdb={'A':1, 'R':2, 'N':3, 'D':4, 'C' :5, 'E':6, 'Q':7, 'G':8, 'H':9, 'O':10, 'I':11, 'L':12, 'K':13, 'M':14, 'F':15, 'P':16, 'U':17, 'S':18, 'T':19, 'W':20, 'Y':21, 'V':22 }
        Propdb={'A':1, 'R':3, 'N':2, 'D':6, 'C' :2, 'E':6, 'Q':4, 'G':6, 'H':3, 'O':6, 'I':6, 'L':1, 'K':3, 'M':1, 'F':5, 'P':2, 'U':6, 'S':2, 'T':2, 'W':5, 'Y':5, 'V':1 }
        newPropMut=list()
        for char in Mutation:
            if char in Propdb:
                newPropMut.append(Propdb[char])
        return newPropMut 
    
    def __str__(self):
        print("Epitope __srt__")
        string=str(self.ID+'\n'+ '\t'+ self.Group+'\n')
        #mutcount=0
        for group in self.DNA:
            string+=str('\t'+'\t'+group+'\t'+str(self.DNA[group])+'\n')
        #string+=str(self.DNA)
        '''for Mutation in self.Mut:
            #for mutstring in Mutation:
                #string+=('\t\t'+mutstring+'\n')
            string+=('\t\t'+Mutation+'\n') 
            #mutcount+=1
            #print("Mutation number "+self.ID+":"+ str(mutcount))'''
        return string
        
     
        
    def Write(self, file) :
        
        for group in self.DNA:
            
            #print(type(group), group)
            
            for Mutation in self.DNA[group]:
                #print(type(Mutation), Mutation, len(Mutation))
                #print(Mutation, len(Mutation))
                file.write(self.ID) 
                file.write("\t")
                file.write(self.Group)
                file.write("\t")
                file.write(str(group))
                file.write(str('\t'+str(Mutation)))
                file.write("\n")
                #file.write(str(Mutation)) 
                
              
        return
    
    def WriteHLA(self, file):
        file.write(str(self.ID+'\t'+self.Group+'\t'+self.Pat+'\t'+ str(self.CountPat)+'\n'))
        #file.write()
        #file.write()
        return
        
def ClearEpitope(Epitope):
    newEpitope=str()
    #print(Epitope)
    for char in Epitope:
        if char not in ['[', "'" ,',', ']']:
            newEpitope+=char
            
    #print(newEpitope)        
    return newEpitope
    
           
        