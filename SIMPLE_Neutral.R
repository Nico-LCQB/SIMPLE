#cat("\n")
#cat("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
#cat("%                                                                                       %\n")
#cat("%       SIMPLE : SImulation de Mutations doubles dans des Populations de LEvures        %\n")
#cat("%                                                                                       %\n")
#cat("%                            Nicolas Agier - V7.10 (Neutral)                            %\n")
#cat("%                                                                                       %\n")
#cat("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
#cat("\n")
#cat("\n")

# In this model :
# - only one of the two daughter cells can mutate
# - use the setwd() command to define the folder in which you want to save the results 
# - the number of mutations at each generation is reported in the verbose.txt file
# - the total number of mutants and mutations produced at each realization is indicated in the res. txt file


###############################
# Simulation parameters #
###############################
    # output: name used for the result files (output.verbose.txt and output.res.txt)
    # r1: rate of the mutation A
    # r2: rate of the mutation B
    # g: number of generations starting from 1 cell (note that more than 27 generations will generate huge files that will require more than 16Gb of RAM)
    # n: number of realizations (correspond to the number of colonies to simulate)

   Output="testNeutral30"
   r1=3.27e-6
   r2=1.61e-7
   g=27
   n=30



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

    # Simulation Loop

        for (i in 1:z)
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
            candidate=tempAB[tempAB==10]

            if(length(candidate)>0)
              { print(paste(length(candidate)," simultaneous mutant at generation:",as.character(i), " for the cycle:",as.character(c)))
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

            # Printing results for each generation
            print(paste("generation: ", i, "res: ",sum(motherA)," ",sum(motherB)," ",length(motherAB[motherAB>1])," ",nmutA," ",nmutB," ",nmutAB))
          
        }

    #------------------------------------------------
    # calculate the number of individuals colonies  |
    #------------------------------------------------

    resAB=motherAB[motherAB>1]
    A=sum(motherA)
    B=sum(motherB)
    C=length(resAB)
    reslist=c(A,B,C,nmutA,nmutB,nmutAB)

    # cleaning temp file
    rm(motherA)        
    rm(motherB)
    rm(motherAB)

    #----------
    # output  |
    #----------
    return(reslist)
}



#################################
# Main loop                     #
#################################


#-----------------------------------------
# creating output files                  |
#-----------------------------------------

   res=matrix(nrow = n, ncol = 6)

   #--------------------------
   # Simulation loop         |
   #--------------------------

   sink(paste(Output,".verbose.txt"), split = T)
   # summary parameters
   print(paste("The rate for mutation A is: ",r1))
   print(paste("The rate for mutation B is: ",r2))
   print(paste("The number of generation is: ",g))

   start=Sys.time()
   for (i in 1:n)
   {
       # verbose
       print(paste("simulation cycle : ",i," on ",n))

       # output results
       S=SIMUI(r1,r2,g,i)
       res[i,1]=S[1]
       res[i,2]=S[2]
       res[i,3]=S[3] 
       res[i,4]=S[4]
       res[i,5]=S[5]
       res[i,6]=S[6]

       # estimating remaining time
       end=Sys.time()
       tdiff=difftime(end,start, units = "mins")                       
       print(paste("calculation will finish in:", round(tdiff/i*(n-i),2)," mins"))
   }

   colnames(res)=c("A","B","AB","mutA","mutB","mutAB")
   write.table(res,paste(Output,".res.txt"),row.names = FALSE, sep = "\t", quote = FALSE)
   sink()
