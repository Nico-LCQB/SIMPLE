#        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#        %                                                                                       %
#        %       SIMPLE : SImulation de Mutations doubles dans des Populations de LEvures        %
#        %                                                                                       %
#        %                            Nicolas Agier - V7.11 (Neutral)                            %
#        %                                                                                       %
#        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# In this model :
# - only one of the two daughter cells can mutate
# - use the setwd() command to define the folder in which you want to save the results 
# - the number of mutations at each generation is reported in the verbose.txt file
# - the total number of mutants and mutations produced at each realization is indicated in the res. txt file
# - this version allows parallel computing
# - this version also gives the generation and the number of events for simultaneous double mutants



###############################
# Simulation parameters       #
###############################
  # output: name used for the result files (output.verbose.txt and output.res.txt)
  # r1: rate of the mutation A
  # r2: rate of the mutation B
  # g: number of generations starting from 1 cell (note that more than 27 generations will generate huge files that will require more than 16Gb of RAM)
  # n: number of realizations (correspond to the number of colonies to simulate)
  # t: numbre of threads to use (to be safe, use the following formula to calulate the number of threads you can use : Max amount of RAM/16Gb)

  Output="Test_version_verbose_4"
  r1=1E-4
  r2=1E-4
  g=27
  n=100
  t=10



###########################
# Simulation Function     #
###########################

SIMUI <- function(x,y,z,c) 
# x=r1, y=r2, z=g, c = cycle

{
    #-----------------------------------------
    # Creating counters and starting objects |
    #-----------------------------------------

    # starting population (1 cell)
    motherA=c(0)
    motherB=c(0)
    motherAB=c(0)

    #-------------------------------------------------------
    # Loop of growth with appearance of mutation A and B   |
    #-------------------------------------------------------

    # Initiate counters of the numbers of mutations
    nmutA=0
    nmutB=0
    nmutAB=0
    ABsimu=0
    

    # Simulation Loop
    for (j in 1:z)
    {
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
            print(paste(length(SimuAB)," simultaneous mutant at generation:",as.character(j), " for the realization:",as.character(c)))
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



################################################
# Run simulation function                      #
################################################

run_simulation <- function(i) {
    message("Running simulation :", i, " on ",n)
    SIMUI(r1,r2,g,i)
    }



#################################
# Main loop                     #
#################################

library(parallel)

#-----------------------------------------
# creating output files                  |
#-----------------------------------------

res=matrix(nrow = n, ncol = 7)

#--------------------------
# Simulation loop         |
#--------------------------

sink(paste(Output,".parameters.txt"), split = T)
    print(paste("The rate for mutation A is: ",r1))
    print(paste("The rate for mutation B is: ",r2))
    print(paste("The number of generation is: ",g))

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
