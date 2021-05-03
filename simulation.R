#################################################
#### Simulation setup for Circular Mediation ####
#################################################

#### Set working directory ####
#file.choose()
#setwd("/Users/skb/Dropbox/Uni/Student Assistant/Irene Klugkist")
getwd()

###########################
#### Preparing Samples ####
###########################

# Prepare true parameter
truea <- c(.1,.2,.4,0)
trueb <- c(.1,.2,.4,0)
truec <- c(.1,.2,.4,0)
truen <- c(30,100,200)
nsim=1000


saveDatasets <- function(truen,truea,trueb,truec,nsim, seed = 140689) {
  set.seed(seed)
  
  # prepare all possible designs
  Alldesigns <- expand.grid(a=truea,b=trueb,c=truec,n=truen,stringsAsFactors=FALSE)
  existingDesigns <- 0
  
  for (i in 1:nrow(Alldesigns)) {
    design <- Alldesigns[i,]
    
    curr_a <- design[,1]
    curr_b <- design[,2]
    curr_c <- design[,3]
    curr_n <- design[,4]
    
    # Prepare saving datasets
    DirName <- paste0(getwd(),
                      "/Data/Datasets_",
                      "n=", curr_n,
                      "a=", curr_a,
                      "b=", curr_b,
                      "c=", curr_c)
    dir.create(DirName, showWarnings = FALSE)
    
    # Simulate n datasets per design
    for (j in 1:nsim) {
      
      filename <- paste0(DirName,"/nr",j,".RDS")
      
      if (!file.exists(filename)) {
        dat <- simData(a=curr_a,b=curr_b,c=curr_c,n=curr_n)
        saveRDS(dat, filename)
      } else {
        existingDesigns <- existingDesigns + 1
      }
      
    }
    
  }
  if (existingDesigns > 0) {
    cat("\n[Data generation: ", existingDesigns, "/", nsim*nrow(Alldesigns),
        " datasets already existed.]\n")
  }
}

system.time(saveDatasets(truen,truea,trueb,truec,nsim))


