#    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#    %                                                                                       %
#    %       SIMPLE : SImulation de Mutations doubles dans des Populations de LEvures        %
#    %                                                                                       %
#    %                           Nicolas Agier - V7.14 (Transient)                           %
#    %                                                                                       %
#    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# In this model :
# - only one of the two daughter cells can mutate.
# - the transient hypermutator population can be simulated over several consecutive generations
# - the generation at which the hypermutator population appear is chosen randomly during the colony development
# - use the setwd() command to define the folder in which you want to save the results 
# - the number of mutations at each generation is reported in the verbose.txt file
# - the total number of mutants and mutations produced at each realization is indicated in the res. txt file
# - this version allows parallel computing
# - this version also gives the generation and the number of events for simultaneous double mutants



##########################
# Simulation parameters  #
##########################
    # output: name used for the result files (output.verbose.txt and output.res.txt)
    # r1: rate of the mutation A
    # r2: rate of the mutation B
    # g: number of generations starting from 1 cell (note that more than 27 generations will generate huge files that will require more than 16Gb of RAM)
    # n: number of realizations (correspond to the number of colonies to simulate)
    # hyp: number of hypermutator cells in the population
    # mult: the fold change increase of mutation rates in the hypermutator subpopulation
    # ngen: number of consecutive hypermutator generations
    # t: numbre of threads to use (to be safe, use the following formula to calulate the number of threads you can use : Max amount of RAM/16Gb)


  Output="WT-DC_lowerSingle_1E5FC350_parallel"
  r1=2.03e-6
  r2=1.256e-7
  g=27
  n=2000
  hyp=100000
  mult=350
  ngen=1
  t=10



############################
# Simulation Functions     #
############################

SIMUeg1 <- function(x,y,z,c,h,m,tgen) 
# function if only one generation of hypermutator
# x=r1, y=r2, z=g, c = cycle, h = size transient pop, m = increase of rates, tgen = generation of transiant hypermutator

{
#-----------------------------------------
# Creating counters and starting objects |
#-----------------------------------------

# starting population (1 cell)
motherA=c(0)
motherB=c(0)
motherAB=c(0)
tgen=sample(1:z, size =  1)
while (h >= 2^(tgen-1))
{
    tgen=sample(1:z, size =  1)
}
print(paste("transient mutators appears at generation: ",tgen, " for the cycle:",as.character(c)))

#-------------------------------------------------------
# Loop of growth with appearance of mutation A and B   |
#-------------------------------------------------------

# Initiate counters of the numbers of mutations
nmutA=0
nmutB=0
nmutAB=0
ABsimu=0

# Simulation Loop
    for (i in 1:z)
    {
        if(i==tgen)
        {
        # calculate the number of mutants expected from division of previous mutants cells
        nmutAi = 2*sum(motherA)
        nmutBi = 2*sum(motherB)
        nmutABi = 2*length(motherAB[motherAB>1])

        # split both population at the mother level
        daughterAnm = sample(0:1,length(motherA)-h,prob = c(1-x,x), replace = T)
        daughterBnm = sample(0:1,length(motherB)-h,prob = c(1-y,y), replace = T)
        daughterAhm = sample(0:1,h,prob = c(1-x*m,x*m), replace = T)
        daughterBhm = sample(0:1,h,prob = c(1-y*m,y*m), replace = T)

        # merge population (normal and transient hypermutator)
        daughterA=c(daughterAnm,daughterAhm) 
        daughterB=c(daughterBnm,daughterBhm)

        # search of double simu in new mutants
        tempAB=daughterA+daughterB
        tempAB=replace(tempAB, tempAB>1,10)
        tempAB=tempAB+motherA+motherB
        SimuAB=tempAB[tempAB==10]

        if(length(SimuAB)>0)
            { print(paste(length(SimuAB)," simultaneous mutant at generation:",as.character(i), " for the cycle:",as.character(c)))
            ABsimu=ABsimu+length(SimuAB)
            }
         
        # Cleaning temp files
        rm(tempAB)
        rm(SimuAB)

        # transfert mother state to daughter state # mutants at previous generation remain mutant
        daughterA = motherA + daughterA
        daughterB = motherB + daughterB
   
        # junction of 2 "daughter" lists
        motherA = c(motherA, daughterA)
        motherB = c(motherB, daughterB) 

        # correct for mutation in mother already mutated
        motherA=replace(motherA, motherA>0,1)
        motherB=replace(motherB, motherB>0,1)

        # table for dble mutant
        motherAB=motherA+motherB

        # Cleaning temp files
        rm(daughterA)        
        rm(daughterB)

        # calculate the number of mutants at this generation # here we sum and remerge both populations
        nmutAf = sum(motherA)
        nmutBf = sum(motherB)
        nmutABf = length(motherAB[motherAB>1])
        nmutA=nmutA+nmutAf-nmutAi
        nmutB=nmutB+nmutBf-nmutBi
        nmutAB=nmutAB+nmutABf-nmutABi

        } else {
        daughterA = sample(0:1,length(motherA),prob = c(1-x,x), replace = T)
        daughterB = sample(0:1,length(motherB),prob = c(1-y,y), replace = T)

        # initiate the number of mutants at this generation
        nmutAi = 2*sum(motherA)
        nmutBi = 2*sum(motherB)
        nmutABi = 2*length(motherAB[motherAB>1])

        # search of double simu in new mutants
        tempAB=daughterA+daughterB
        tempAB=replace(tempAB, tempAB>1,10)
        tempAB=tempAB+motherA+motherB
        SimuAB=tempAB[tempAB==10]

        if(length(SimuAB)>0)
            { print(paste(length(SimuAB)," simultaneous mutant at generation:",as.character(i), " for the cycle:",as.character(c)))
            ABsimu=ABsimu+length(SimuAB)
            }
         
        # cleaning temp file
        rm(tempAB)

        # transfert mother state to daugther state
        daughterA = motherA + daughterA
        daughterB = motherB + daughterB
    
        # junction of 2 "daughter" lists
        motherA = c(motherA, daughterA)
        motherB = c(motherB, daughterB) 

        # correct for mutation in mother already mutated
        motherA=replace(motherA, motherA>0,1)
        motherB=replace(motherB, motherB>0,1)

        # table for dble mutant
        motherAB=motherA+motherB

        # cleaning temp file
        rm(daughterA)        
        rm(daughterB)

        # calculate the number of mutants at this generation
        nmutAf = sum(motherA)
        nmutBf = sum(motherB)
        nmutABf = length(motherAB[motherAB>1])
        nmutA=nmutA+nmutAf-nmutAi
        nmutB=nmutB+nmutBf-nmutBi
        nmutAB=nmutAB+nmutABf-nmutABi
        }
    }

#------------------------------------------------
# calculate the number of individuals colonies  |
#------------------------------------------------

resAB=motherAB[motherAB>1]
A=sum(motherA)
B=sum(motherB)
C=length(resAB)
reslist=c(A,B,C,nmutA,nmutB,nmutAB,ABsimu)

# cleaning temp file
rm(motherA)        
rm(motherB)
rm(motherAB)

#----------
# output  |
#----------
return(reslist)
}

SIMUsup1 <- function(x,y,z,c,h,m,tgen, ngen) 
# function if more than one generation of hypermutator
# x=r1, y=r2, z=g, c = cycle, h = size transient pop, m = increase of rates, tgen = 1st generation of transiant hypermutator, ngen=number of consecutive hypermutator generation

{
    #-----------------------
    # choosing random tgen |
    #-----------------------

    # Starting population (1 cell)
    motherA=c(0)
    motherB=c(0)
    motherAB=c(0)
    tgen=sample(1:(z-ngen), size =  1)
    while (h >= 2^(tgen-1))
    {
        tgen=sample(1:(z-ngen), size =  1)
    }
    print(paste("transient mutators appears at generation: ",tgen, " for the cycle:",as.character(c)))

    #-----------------------------------------------------
    # Loop of growth with appearance of mutation A and B |
    #-----------------------------------------------------
    # Initiate counter for the number of mutations
    nmutA=0
    nmutB=0
    nmutAB=0
    ABsimu=0

    # Simulation Loop
    for (i in 1:z)
    {
        # generation of entering in hypermutator state    
        if(i==tgen) # divide in two subpopulation (standard and hypermutator)  
        {
            # calculate the number of mutants expected from division of previous mutants cells
            nmutAi = 2*sum(motherA)
            nmutBi = 2*sum(motherB)
            nmutABi = 2*length(motherAB[motherAB>1])

            # split both population at the mother level
            popsizen=length(motherA)-h
            motherAnm = motherA[1:popsizen]
            motherBnm = motherB[1:popsizen]
            motherAhm = motherA[as.numeric(popsizen+1):length(motherA)]
            motherBhm = motherB[as.numeric(popsizen+1):length(motherB)]

            # split standard and hypermutator daughter cells and sample mutants
            daughterAnm = sample(0:1,length(motherA)-h,prob = c(1-x,x), replace = T) 
            daughterBnm = sample(0:1,length(motherB)-h,prob = c(1-y,y), replace = T) 
            daughterAhm = sample(0:1,h,prob = c(1-x*m,x*m), replace = T)
            daughterBhm = sample(0:1,h,prob = c(1-y*m,y*m), replace = T)

            # Find the simultaneous double mutants in the daughter cells (standard population)
            tempABnm=daughterAnm+daughterBnm
            tempABnm=replace(tempABnm, tempABnm>1,10) # set double (A=1 and B=1) to 10 to discriminate later with other case (mother=1 et daugther =1, flag later as stay mutant)
            tempABnm=tempABnm+motherAnm+motherBnm # add maternal character to separate simultaneous from sequentials (check that A and B mutations only appear at that generation)
            SimuABn=tempABnm[tempABnm==10]
            if(length(SimuABn)>0) {
                print(paste(length(SimuABn)," simultaneous mutant at generation:",as.character(i), " for the cycle:",as.character(c)))
                ABsimu=ABsimu+length(SimuABn)
                }

            # Find the simultaneous double mutants in the daughter cells (hypermutator population population)
            tempABhm=daughterAhm+daughterBhm
            tempABhm=replace(tempABhm, tempABhm>1,10)
            tempABhm=tempABhm+motherAhm+motherBhm 
            SimuABh=tempABhm[tempABhm==10]
            if(length(SimuABh)>0) {
                print(paste(length(SimuABh)," simultaneous mutant at generation:",as.character(i), " for the cycle:",as.character(c)))
                ABsimu=ABsimu+length(SimuABh)
                }

            # Cleaning temp files
            rm(tempABnm)
            rm(tempABhm)
            rm(SimuABn)
            rm(SimuABh)
                
            # transfert mother state to daughter state # mutants at previous generation remain mutant
            daughterAnm = motherAnm + daughterAnm
            daughterBnm = motherBnm + daughterBnm
            daughterAhm = motherAhm + daughterAhm
            daughterBhm = motherBhm + daughterBhm

            # junction of 2 "daughter" lists
            motherAnm = c(motherAnm, daughterAnm)
            motherBnm = c(motherBnm, daughterBnm) 
            motherAhm = c(motherAhm, daughterAhm)
            motherBhm = c(motherBhm, daughterBhm) 

            # correct for mutation in mother already mutated (which should have a value of 2)
            motherAnm=replace(motherAnm, motherAnm>0,1)
            motherBnm=replace(motherBnm, motherBnm>0,1)
            motherAhm=replace(motherAhm, motherAhm>0,1)
            motherBhm=replace(motherBhm, motherBhm>0,1)

            # table for dble mutant
            motherABnm=motherAnm+motherBnm
            motherABhm=motherAhm+motherBhm

            # cleaning temp files
            rm(daughterAnm)        
            rm(daughterBnm)
            rm(daughterAhm)        
            rm(daughterBhm)

            # calculate the number of mutants at this generation # here we sum and remerge both populations
            nmutAf = sum(motherAnm)+sum(motherAhm)
            nmutBf = sum(motherBnm)+sum(motherBhm)
            nmutABf = length(motherABnm[motherABnm>1])+length(motherABhm[motherABhm>1])
            nmutA=nmutA+nmutAf-nmutAi
            nmutB=nmutB+nmutBf-nmutBi
            nmutAB=nmutAB+nmutABf-nmutABi

            # Break if generation i = 27
            if(i==27){
                motherA = c(motherAnm,motherAhm)
                motherB = c(motherBnm,motherBhm)
                motherAB = c(motherABnm,motherABhm)
            }

        } else if (i>tgen && i<tgen+ngen) {
            # hypermutator state    
            # both populations (standard and hypermutator are treated separately for the number of chosen generations)
            # To initiate the number of mutants at this generation
            nmutAi = 2*sum(motherAnm)+2*sum(motherAhm)
            nmutBi = 2*sum(motherBnm)+2*sum(motherBhm)
            nmutABi = 2*length(motherABnm[motherABnm>1])+2*length(motherABhm[motherABhm>1])

            # split standard and hypermutator daughter cells and sample mutants
            daughterAnm = sample(0:1,length(motherAnm),prob = c(1-x,x), replace = T)
            daughterBnm = sample(0:1,length(motherBnm),prob = c(1-y,y), replace = T)
            daughterAhm = sample(0:1,length(motherAhm),prob = c(1-x*m,x*m), replace = T)
            daughterBhm = sample(0:1,length(motherBhm),prob = c(1-y*m,y*m), replace = T)

            # Find the double in the daughter cells (standard and hypermutators)
            tempABnm=daughterAnm+daughterBnm
            tempABnm=replace(tempABnm, tempABnm>1,10)
            tempABnm=tempABnm+motherAnm+motherBnm 
            SimuABn=tempABnm[tempABnm==10]
            if(length(SimuABn)>0) {
                print(paste(length(SimuABn)," simultaneous mutant at generation:",as.character(i), " for the cycle:",as.character(c)))
                ABsimu=ABsimu+length(SimuABn)
            }
            tempABhm=daughterAhm+daughterBhm
            tempABhm=replace(tempABhm, tempABhm>1,10)
            tempABhm=tempABhm+motherAhm+motherBhm 
            SimuABh=tempABhm[tempABhm==10]
            if(length(SimuABh)>0) {
                print(paste(length(SimuABh)," simultaneous mutant at generation:",as.character(i), " for the cycle:",as.character(c)))
                ABsimu=ABsimu+length(SimuABh)
            }

            # cleaning temp files
            rm(tempABnm)
            rm(tempABhm)
            rm(SimuABn)
            rm(SimuABh)


            # transfert mother state to daughter state
            daughterAnm = motherAnm + daughterAnm
            daughterBnm = motherBnm + daughterBnm
            daughterAhm = motherAhm + daughterAhm
            daughterBhm = motherBhm + daughterBhm

            # junction of 2 "daughter" lists
            motherAnm = c(motherAnm, daughterAnm)
            motherBnm = c(motherBnm, daughterBnm) 
            motherAhm = c(motherAhm, daughterAhm)
            motherBhm = c(motherBhm, daughterBhm) 

            # correct for mutation in mother already mutated
            motherAnm=replace(motherAnm, motherAnm>0,1)
            motherBnm=replace(motherBnm, motherBnm>0,1)
            motherAhm=replace(motherAhm, motherAhm>0,1)
            motherBhm=replace(motherBhm, motherBhm>0,1)

            # table for dble mutant
            motherABnm=motherAnm+motherBnm
            motherABhm=motherAhm+motherBhm

            # cleaning temp files
            rm(daughterAnm)        
            rm(daughterBnm)
            rm(daughterAhm)        
            rm(daughterBhm)

            # calculate the number of mutations at this generation # here we sum and remerge both populations
            nmutAf = sum(motherAnm)+sum(motherAhm)
            nmutBf = sum(motherBnm)+sum(motherBhm)
            nmutABf = length(motherABnm[motherABnm>1])+length(motherABhm[motherABhm>1])
            nmutA=nmutA+nmutAf-nmutAi
            nmutB=nmutB+nmutBf-nmutBi
            nmutAB=nmutAB+nmutABf-nmutABi

            # Break if generation i = 27
            if(i==27){
                motherA = c(motherAnm,motherAhm)
                motherB = c(motherBnm,motherBhm)
                motherAB = c(motherABnm,motherABhm)
            }

        } else if (i==tgen+ngen) {
            # Getting out of hypermutator state (we should remerge nm and hm population at that stage, to check)    
            # initiate the number of mutants at this generation
            nmutAi = 2*sum(motherAnm)+2*sum(motherAhm)
            nmutBi = 2*sum(motherBnm)+2*sum(motherBhm)
            nmutABi = 2*length(motherABnm[motherABnm>1])+2*length(motherABhm[motherABhm>1])

            # Do not split the populations but keep with the one created in the previous commands
            daughterAnm = sample(0:1,length(motherAnm),prob = c(1-x,x), replace = T)
            daughterBnm = sample(0:1,length(motherBnm),prob = c(1-y,y), replace = T)
            daughterAhm = sample(0:1,length(motherAhm),prob = c(1-x*m,x*m), replace = T)
            daughterBhm = sample(0:1,length(motherBhm),prob = c(1-y*m,y*m), replace = T)

            # merge population (normal and transient hypermutator)
            daughterA=c(daughterAnm,daughterAhm) 
            daughterB=c(daughterBnm,daughterBhm)
            motherA=c(motherAnm,motherAhm) # update of the list taking into account the hypermutators cycles
            motherB=c(motherBnm,motherBhm) # update of the list taking into account the hypermutators cycles

            # search of double simu in new mutants
            tempAB=daughterA+daughterB
            tempAB=replace(tempAB, tempAB>1,10) 
            tempAB=tempAB+motherA+motherB 
            SimuAB=tempAB[tempAB==10]

            if(length(SimuAB)>0) {
                print(paste(length(SimuAB)," simultaneous mutant at generation:",as.character(i), " for the cycle:",as.character(c)))
                ABsimu=ABsimu+length(SimuAB)
            }
         
            # cleaning temp files
            rm(tempAB)
            rm(SimuAB)

            # transfert mother state to daugther state
            daughterA = motherA + daughterA
            daughterB = motherB + daughterB
   
            # junction of 2 "daughter" lists
            motherA = c(motherA, daughterA)
            motherB = c(motherB, daughterB) 

            # correct for mutation in mother already mutated
            motherA=replace(motherA, motherA>0,1)
            motherB=replace(motherB, motherB>0,1)

            # table for dble mutant
            motherAB=motherA+motherB

            # cleaning temp file
            rm(daughterA)        
            rm(daughterB)

            # calculate the number of mutants at this generation
            nmutAf = sum(motherA)
            nmutBf = sum(motherB)
            nmutABf = length(motherAB[motherAB>1])
            nmutA=nmutA+nmutAf-nmutAi
            nmutB=nmutB+nmutBf-nmutBi
            nmutAB=nmutAB+nmutABf-nmutABi
        
        } else {
            # Standard mutation regime    
            daughterA = sample(0:1,length(motherA),prob = c(1-x,x), replace = T)
            daughterB = sample(0:1,length(motherB),prob = c(1-y,y), replace = T)
  
            # initiate the number of mutants at this generation
            nmutAi = 2*sum(motherA)
            nmutBi = 2*sum(motherB)
            nmutABi = 2*length(motherAB[motherAB>1])

            # search of double simu in new mutants
            tempAB=daughterA+daughterB
            tempAB=replace(tempAB, tempAB>1,10)
            tempAB=tempAB+motherA+motherB 
            SimuAB=tempAB[tempAB==10]
 
            if(length(SimuAB)>0) {
                print(paste(length(SimuAB)," simultaneous mutant at generation:",as.character(i), " for the cycle:",as.character(c)))
                ABsimu=ABsimu+length(SimuAB)
            }
         
            # cleaning temp file
            rm(tempAB)
            rm(SimuAB)

            # transfert mother state to daugther state
            daughterA = motherA + daughterA
            daughterB = motherB + daughterB
    
            # junction of 2 "daughter" lists
            motherA = c(motherA, daughterA)
            motherB = c(motherB, daughterB) 

            # correct for mutation in mother already mutated
            motherA=replace(motherA, motherA>0,1)
            motherB=replace(motherB, motherB>0,1)

            # table for dble mutant
            motherAB=motherA+motherB
  
            # cleaning temp files
            rm(daughterA)        
            rm(daughterB)

            # calculate the number of mutants at this generation
            nmutAf = sum(motherA)
            nmutBf = sum(motherB)
            nmutABf = length(motherAB[motherAB>1])
            nmutA=nmutA+nmutAf-nmutAi
            nmutB=nmutB+nmutBf-nmutBi
            nmutAB=nmutAB+nmutABf-nmutABi
        
        }
    }

#-----------------------------------------------
# calculate the number of individuals colonies  |
#------------------------------------------------

resAB=motherAB[motherAB>1]
A=sum(motherA)
B=sum(motherB)
C=length(resAB)
reslist=c(A,B,C,nmutA,nmutB,nmutAB,ABsimu)

# nettoyage des fichiers temp
rm(motherA)        
rm(motherB)
rm(motherAB)

#----------
# output  |
#----------
return(reslist)
}



################################################
# Run simulation function                      #
################################################

run_simulation <- function(i) {
    message("Running simulation :", i, " on ",n)
    if(ngen==1) {
        S=SIMUeg1(r1,r2,g,i,hyp,mult,tgen)
    } else {
        ngen=ngen-1 #remove the first generation which is not taking into account in the number of gen
        S=SIMUsup1(r1,r2,g,i,hyp,mult,tgen,ngen)
    }
}



#################################
# Main loop                     #
#################################

library(parallel)

#-----------------------------------------
# creating output files                   |
#-----------------------------------------

res=matrix(nrow = n, ncol = 7)

#--------------------------
# Simulation loop          |
#--------------------------

sink(paste(Output,".verbose.txt"), split = T)

# summary parameters
print(paste("The rate for mutation A is: ",r1))
print(paste("The rate for mutation B is: ",r2))
print(paste("The number of generation is: ",g))
print(paste(hyp," hypermutator cells (x",mult,") for ",ngen," generations")) 

# Run simulations in parallel
results <- mclapply(1:n, run_simulation, mc.cores = t)

# Redirect output to a file
for (i in 1:n) {
    S <- results[[i]]
    if(!is.null(S)){
       res[i,1]=S[1]
       res[i,2]=S[2]
       res[i,3]=S[3] 
       res[i,4]=S[4]
       res[i,5]=S[5]
       res[i,6]=S[6]
       res[i,7]=S[7]
    }
}

colnames(res)=c("A","B","AB","mutA","mutB","mutAB","mutABSimu")
write.table(res,paste(Output,".res.txt"),row.names = FALSE, sep = "\t", quote = FALSE)
sink()
