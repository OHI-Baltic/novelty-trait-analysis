
data_dir <- here::here("Untitled", "data")

library(vegan)
require(analogue)

library(EnvCpt)

###########################################################
#Degree of Novelty analysis
##########################################################
####this function calculates novelty and change
Novelty_function <- function(M, start_baseline, end_baseline, end_year, n, method){
  
  M <- as.data.frame(M)
  M <- cbind(M[,(1:2)], apply(M[,-(1:2)], 2,as.numeric))
  
  #M <- M[apply(M[,-c(1,2)], 1, function(x) !all(x==0)),]
  
  Ana_result = data.frame(cma_dist=NA, TS.ind=NA, FS.ind=NA, TS.area=NA, TS.year=NA, FS.area=NA, FS.year=NA)
  Ana_result <- Ana_result[-1,]
  
  
  FS = M[which(M[,1] %in% c(start_baseline:end_baseline)),]
  
  ## First Select only columns with pollen %
  
  
  FS.a <- FS[,-c(1:2)]
  names.fs <- rownames(FS.a)
  n.fs <- length(names.fs)
  
  ## Rbind two datasets into one
  
  
  for (t in (end_baseline+1):end_year) { 
      
    
    TS = M[which(M[,1]==t),]
    print (paste("Doing training set TS",t, "...patience..."))
    
   
    TS.a <- TS[,-c(1:2)]
   
      
    names.ts <- rownames(TS.a)
    n.ts <- length(names.ts)
    S.a <- rbind(TS.a, FS.a)
    
  
    if(method=="Hellinger"){ #### helinger distance dissimilarity based on the composition of the area=====
      
      S.a2 <- S.a
      for( i in 1:ncol(S.a)){
        for(j in 1:nrow(S.a)){
          S.a2[j,i] <- S.a[j,i]/rowSums(S.a[j,])
        }
      }
      
      S.a2[which(rowSums(TS.a)==FALSE),] <- 0
      
      S.d2 <- as.matrix(distance(S.a2, method="chord", dist=T))
      
      S.dd <- S.d2
      
      S.dd <- matrix( S.dd[ 1 : n.ts, (n.ts+1) : (n.ts+n.fs) ], nrow=n.ts, ncol=length((n.ts+1) : (n.ts+n.fs)))
      
      colnames(S.dd) <- rownames(FS.a)
      rownames(S.dd) <- rownames(TS.a)
      
    }
    
    else if (method=="EDlog"){  #### Euclidean distance dissimilarity based on the stock size of the area=====
      
      
        
      S.a1<- S.a
      #for( i in  1:ncol(S.a)){
      #  for(j in 1:nrow(S.a)){
      #    S.a1[j,i] <- (S.a[j,i])/max.col[i]#(S.a[j,i]- sd.mean[i])/sd.col[i]
      #  }
      #}
      for( i in  1:ncol(S.a)){  ## for to reduce to the log and divide by the numver of variables to avoid the influence of bigger dataset
        for(j in 1:nrow(S.a)){
          S.a1[j,i] <- log(S.a[j,i]+1)#/(ncol(M)-2)#(S.a[j,i]- sd.mean[i])/sd.col[i]
        }
      }
      
      
      S.d1 <- as.matrix(dist(S.a1, method="euclidean"))
      
      
      S.dd <- S.d1 
      
      
      S.dd <- matrix( S.dd[ 1 : n.ts, (n.ts+1) : (n.ts+n.fs) ], nrow=n.ts, ncol=length((n.ts+1) : (n.ts+n.fs)))
      
      colnames(S.dd) <- rownames(FS.a)
      rownames(S.dd) <- rownames(TS.a)
      
    }

    close <- vector("list", n.ts)
    for(i in seq_len(n.ts)) {
      if(is.matrix(S.dd)){
        close[[i]] <- sort(S.dd[i, ])
      }else{
        close[[i]] <- sort(S.dd[i])
      }
    }
    
    
    names(close) <- names.ts 
    
    Fuse.cma <- list(close)
    names(Fuse.cma)[1] <- "close"
    
    .call <- match.call()
    .call[[1]] <- as.name("cma")
    structure(list(close = close,
                   call = .call, 
                   quant = "none", probs = "none",
                   prob = "none",
                   method = "fuse",
                   n.analogs = "all"),
              class = "cma")
    class(Fuse.cma) <- "cma"
    
    Ana.cma <- Fuse.cma
    
    
    print (paste("...checking closest analogues in baseline ...patience..."))
    
    
    ## Extract the cma distance for each TS sample ###
    ## and writes them in data.frame Close.dat
    Close.dat = data.frame(cma_dist=NA, TS.ind=NA, FS.ind=NA,TS.area=NA, TS.year=NA)
    Close.dat <- as.data.frame(Close.dat)
    TS.ind = row.names(TS) # list to find the samples back
    TS.n = length(Ana.cma$close)
    for (i in 1:TS.n) {
      k = TS.ind[i]
      if(length(Ana.cma[["close"]] [[k]])>0){
        temp = as.data.frame(Ana.cma[["close"]] [[k]] [1])
        colnames(temp)[1] = "cma_dist"
        temp$TS.ind = as.numeric(TS.ind[i])
        if (nrow(FS)==1){
          temp$FS.ind = as.numeric(row.names(FS))
        }else{
          temp$FS.ind = as.numeric(row.names(temp))
        }
        
        temp$TS.area = TS[i,2]
        temp$TS.year = TS[i,1]
        Close.dat = rbind(Close.dat, temp)
      }
    }
    
    
    if(nrow(Close.dat)>1){
      Close.dat = Close.dat[-1, ]
      row.names(Close.dat) = seq(nrow(Close.dat)) # delete unused row.names 
      FS.ind = row.names(FS)            
      for (i in 1:nrow(Close.dat)) {
        k = as.numeric(Close.dat$FS.ind[i])
        Close.dat$FS.area [i] = as.character(FS[which(rownames(FS)==k),2])
        Close.dat$FS.year [i] = as.character(FS[which(rownames(FS)==k),1])
      }
      
      Ana_result <- rbind(Ana_result, Close.dat)
    }
    
  }
  print(Ana_result)
}

##########################################################

#####things to change!

path <-  data_dir

##choose method for composition "Hellinger" or stock size ""EDlog"
method= "Hellinger"

##choose number of years in target bin 
n=1

##########################################################


 ###phytoplankton   ######## 
load(paste0(path,'/phyto_class_mean_approx_5mean.Rdata')) ##phytoplankton mean value by class, with approximation to reduced NAs and a moving average of 5 years
phyto <- phyto_class_mean_approx_5mean
colnames(phyto)[2] <- "helcom_id"
colnames(phyto)[1] <- "year"
phyto <- cbind(phyto[,(1:2)], apply(phyto[,-(1:2)], 2,as.numeric))

#deleting almost empty columns 
phyto <- phyto[,-which(colnames(phyto) %in% c("V27","Mamiellophyceae", "Thecofilosea", "Oligotrichea", "Bicoecea"))]
for ( i in 3:ncol(phyto)){
  for ( j in 1:nrow(phyto)){
    if(is.na(phyto[j,i])){
      phyto[j,i] <- 0
    }
  }
}
View(phyto) ####check the data

##deleting areas with lots of missing data 
phyto <-  phyto[-which( phyto$helcom_id %in% c("SEA-001","SEA-002", "SEA-003","SEA-004","SEA-005","SEA-006",
                                               "SEA-008", "SEA-010", "SEA-011", "SEA-014", "SEA-016")),]

##Degree of Novelty analysis ====
## P.S. the stating and ending years can be changed
#the end_year here is the last year available in the dataset
phyto_novelty <- Novelty_function(M=phyto, start_baseline=1979, end_baseline=1985, end_year=2017, n, method) #n is the number of years by target bin
                                                                                                             # the baseline and the end year can be changed
## degree of change analysis ====
## P.S. the stating and ending years can be changed
phyto_change <- as.data.frame(phyto_novelty[1,])
phyto_change <- phyto_change[-1,]
helcom <- unique(phyto$helcom_id)
for (i in helcom){
  data <- phyto[which(phyto$helcom_id==i),]
  phyto_change <- rbind(phyto_change, 
                        Novelty_function(M=data, start_baseline=1979, end_baseline=1985, end_year=2017, n, method))
}


### zooplankton  ########

load(paste0(path,'Zooplankton_SYKE_BB_mean_Helcom_fix_approx_5mean.Rdata')) ##zoopkankton data from SYKE including Bornholm Basin (BB), mean by genera and Helcom Subdivision, approximation to reduce NAs and 5 years moving average
zoo <- Zooplankton_SYKE_BB_mean_Helcom_fix_approx_5mean

##columns used in the analysis
zoo <- zoo[, which(colnames(zoo) %in%  c("year", "helcom_id", "Synchaeta", "Keratella", "Asplanchna", "Euchlanis", 
                                         "Collotheca"  , "Polyarthra","Kellicottia" , "Filinia", "Notholca", "Bosmina",
                                         "Evadne",  "Pleopsis", "Podon", "Daphnia","Cercopagis", 
                                         "Ceriodaphnia",   "Leptodora", "Bythotrephes", "Eubosmina", "Polyphemus", 
                                         "Chydorus","Acartia","Centropages", "Eurytemora","Limnocalanus","Pseudocalanus", 
                                         "Temora", "Cyclops", "Oithona", "Calanoida", "Mysis", "Harmothoe", "Pygospio", 
                                         "Capitella","Marenzelleria", "Mytilus",  "Balanus", "Amphibalanus","Electra", 
                                         "Pontoporeia", "Monoporeia","Fritillaria"))]


zoo <- cbind(zoo[,(1:2)], apply(zoo[,-(1:2)], 2,as.numeric))

for(i in 1:nrow(zoo)){
  for(j in 3:ncol(zoo)){
    if(is.na(zoo[i,j])){
      zoo[i,j] <- 0
    }
  }
}
zoo <- as.data.frame(zoo)

#delete areas with lots of missing data
zoo <- zoo[-which(zoo$helcom_id %in% c("SEA-006", "SEA-007", "SEA-014","SEA-016")),]
zoo <- cbind(zoo[,(1:2)], apply(zoo[,-(1:2)], 2,as.numeric))

##Degree of Novelty analysis ====
## P.S. the stating and ending years can be changed
#the end_year here is the last year available in the dataset
zoo_novelty <- Novelty_function(M=zoo, start_baseline=1979, end_baseline=1985, end_year=2016, n, method) 

## degree of change analysis ====
## P.S. the stating and ending years can be changed
zoo_change <- as.data.frame(zoo_novelty[1,])
zoo_change <- zoo_change[-1,]
helcom <- unique(zoo$helcom_id)
for (i in helcom){
  data <- zoo[which(zoo$helcom_id==i),]
  zoo_change <- rbind(zoo_change, 
                        Novelty_function(M=data, start_baseline=1979, end_baseline=1985, end_year=2016, n, method))
}

#####fish  ########

load(paste0(path,"fish_Helcom_fix_approx_5mean.Rdata"))
fish <- fish_Helcom_fix_approx_5mean 

fish <- cbind(fish[,(1:2)], apply(fish[,-(1:2)], 2,as.numeric))

#homogenize the dataset
fish[,which(colnames(fish) %in% c("Chelidonichthys.lucerna"))] <- rowSums(fish[,which(colnames(fish) %in% c("Chelidonichthys.lucerna", "Chelidonichthys.lucernus"))], na.rm=TRUE)

# the columns used 
fish <- fish[,which(colnames(fish) %in% c("year", "helcom_id","Gadus.morhua", "Clupea.harengus", "Platichthys.flesus",
                                          "Sprattus.sprattus",  "Salmo.salar", "Limanda.limanda", "Pleuronectes.platessa",
                                          "Psetta.maxima", "Hippoglossoides.platessoides","Scophthalmus.rhombus",
                                          "Arnoglossus.laterna", "Solea.solea", "Scophthalmus.maximus","Microstomus.kitt",
                                          "Glyptocephalus.cynoglossus", "Lepidorhombus.whiffiagonis", "Zeugopterus.punctatus",
                                          "Ctenolabrus.rupestris", "Buglossidium.luteum", "Hippoglossus.hippoglossus",
                                          "Myoxocephalus.scorpius","Cyclopterus.lumpus", "Agonus.cataphractus", "Eutrigla.gurnardus",
                                          "Chelidonichthys.cuculus",  "Chelidonichthys.lucerna", "Liparis.liparis",
                                          "Triglopsis.quadricornis",  "Taurulus.bubalis", "Myoxocephalus.quadricornis" ,
                                          "Enchelyopus.cimbrius", "Pollachius.virens","Merlangius.merlangus", "Pollachius.pollachius",
                                          "Trisopterus.luscus", "Merluccius.merluccius", "Melanogrammus.aeglefinus",
                                          "Trisopterus.minutus", "Trisopterus.esmarkii", "Molva.molva",  "Pomatoschistus.minutus", 
                                          "Gobius.niger", "Neogobius.melanostomus", "Aphia.minuta","Gobiusculus.flavescens", 
                                          "Salmo.trutta",  "Coregonus.lavaretus", "Anguilla.anguilla",  "Gasterosteus.aculeatus",
                                          "Spinachia.spinachia", "Pungitius.pungitius", "Engraulis.encrasicolus",   "Alosa.fallax",
                                          "Sardina.pilchardus", "Alosa.agone", "Syngnathus.typhle", "Syngnathus.rostellatus",
                                           "Syngnathus.acus",  "Entelurus.aequoreus", "Nerophis.ophidion", "Lampetra.fluviatilis", 
                                          "Petromyzon.marinus", "Esox.lucius" ))]

fish <- as.data.frame(fish)
fish <- cbind(fish[,(1:2)], apply(fish[,-(1:2)], 2,as.numeric))
fish$year <- as.numeric(fish$year)
###replace NAs with 0
for(i in 1:nrow(fish)){
  for(j in 1:ncol(fish)){
    if(is.na(fish[i,j])){
      fish[i,j] <- 0
    }
  }
}

##delete areas with lots of missing data
fish <- fish[-which(fish$helcom_id %in% c("SEA-002", "SEA-003","SEA-004", "SEA-005","SEA-008" ,"SEA-011", "SEA-012" )),]

##Degree of Novelty analysis ====
## P.S. the stating and ending years can be changed
#the end_year here is the last year available in the dataset
fish_novelty <- Novelty_function(M=fish, start_baseline=1990, end_baseline=1995, end_year=2015, n, method)  #n is the number of years by target bin

## degree of change analysis ====
## P.S. the stating and ending years can be changed
fish_change <- as.data.frame(fish_novelty[1,])
fish_change <- fish_change[-1,]
helcom <- unique(fish$helcom_id)
for (i in helcom){
  data <- fish[which(fish$helcom_id==i),]
  fish_change <- rbind(fish_change, 
                      Novelty_function(M=data, start_baseline=1990, end_baseline=1995, end_year=2015, n, method))
}


###### combine phytoplankton and zooplankton  ########

load(paste0(path,'phyto_class_mean_approx_5mean.Rdata'))
phyto <- phyto_class_mean_approx_5mean
colnames(phyto)[2] <- "helcom_id"
colnames(phyto)[1] <- "year"
phyto <- cbind(phyto[,(1:2)], apply(phyto[,-(1:2)], 2,as.numeric))
#phyto <- phyto[,1:11]
phyto <- phyto[,-which(colnames(phyto) %in% c("V27","Mamiellophyceae", "Thecofilosea", "Oligotrichea", "Bicoecea"))]
for ( i in 3:ncol(phyto)){
  for ( j in 1:nrow(phyto)){
    if(is.na(phyto[j,i])){
      phyto[j,i] <- 0
    }
  }
}
phyto <- phyto[which(phyto$helcom_id %in% c("SEA-013", "SEA-017", "SEA-009", "SEA-012", "SEA-015")),]


load(paste0(path,'Zooplankton_SYKE_BB_mean_Helcom_fix_approx_5mean.Rdata'))
zoo <- Zooplankton_SYKE_BB_mean_Helcom_fix_approx_5mean

zoo <- zoo[, which(colnames(zoo) %in%  c("year", "helcom_id", "Synchaeta", "Keratella", "Asplanchna", "Euchlanis", 
                                         "Collotheca"  , "Polyarthra","Kellicottia" , "Filinia", "Notholca", "Bosmina",
                                         "Evadne",  "Pleopsis", "Podon", "Daphnia","Cercopagis", 
                                         "Ceriodaphnia",   "Leptodora", "Bythotrephes", "Eubosmina", "Polyphemus", 
                                         "Chydorus","Acartia","Centropages", "Eurytemora","Limnocalanus","Pseudocalanus", 
                                         "Temora", "Cyclops", "Oithona", "Calanoida", "Mysis", "Harmothoe", "Pygospio", 
                                         "Capitella","Marenzelleria", "Mytilus",  "Balanus", "Amphibalanus","Electra", 
                                         "Pontoporeia", "Monoporeia","Fritillaria"))]


zoo <- cbind(zoo[,(1:2)], apply(zoo[,-(1:2)], 2,as.numeric))

for(i in 1:nrow(zoo)){
  for(j in 3:ncol(zoo)){
    if(is.na(zoo[i,j])){
      zoo[i,j] <- 0
    }
  }
}
zoo <- as.data.frame(zoo)

zoo <- zoo[which(zoo$helcom_id %in% c("SEA-013", "SEA-017", "SEA-009", "SEA-012", "SEA-015")),]

phyto_zoo <- merge(phyto, zoo, by=c("year","helcom_id"))

phyto_zoo <- cbind(phyto_zoo[,(1:2)], apply(phyto_zoo[,-(1:2)], 2,as.numeric))


##Degree of Novelty analysis ====
## P.S. the stating and ending years can be changed
#the end_year here is the last year available in the dataset
phyto_zoo_novelty <- Novelty_function(M=phyto_zoo, start_baseline=1979, end_baseline=1985, end_year=2016, n, method)  #n is the number of years by target bin

## degree of change analysis ====
## P.S. the stating and ending years can be changed
phyto_zoo_change <- as.data.frame(phyto_zoo_novelty[1,])
phyto_zoo_change <- phyto_zoo_change[-1,]
helcom <- unique(phyto_zoo$helcom_id)
for (i in helcom){
  data <- phyto_zoo[which(phyto_zoo$helcom_id==i),]
  phyto_zoo_change <- rbind(phyto_zoo_change, 
                      Novelty_function(M=data, start_baseline=1979, end_baseline=1985, end_year=2016, n, method))
}

###### phytoplankton, zooplankton and fish ########

##SEA-009 is the only common area for the three groups
phyto_zoo <- phyto_zoo[which(phyto_zoo$helcom_id == "SEA-009"),]
phyto_zoo <- phyto_zoo[which(phyto_zoo$year %in% c(1990:2015)),]

load(paste0(path,"fish_Helcom_fix_approx_5mean.Rdata"))
fish <- fish_Helcom_fix_approx_5mean 

fish <- cbind(fish[,(1:2)], apply(fish[,-(1:2)], 2,as.numeric))

fish[,which(colnames(fish) %in% c("Chelidonichthys.lucerna"))] <- rowSums(fish[,which(colnames(fish) %in% c("Chelidonichthys.lucerna", "Chelidonichthys.lucernus"))], na.rm=TRUE)

fish <- fish[,which(colnames(fish) %in% c("year", "helcom_id","Gadus.morhua", "Clupea.harengus", "Platichthys.flesus",
                                          "Sprattus.sprattus",  "Salmo.salar", "Limanda.limanda", "Pleuronectes.platessa",
                                          "Psetta.maxima", "Hippoglossoides.platessoides","Scophthalmus.rhombus",
                                          "Arnoglossus.laterna", "Solea.solea", "Scophthalmus.maximus","Microstomus.kitt",
                                          "Glyptocephalus.cynoglossus", "Lepidorhombus.whiffiagonis", "Zeugopterus.punctatus",
                                          "Ctenolabrus.rupestris", "Buglossidium.luteum", "Hippoglossus.hippoglossus",
                                          "Myoxocephalus.scorpius","Cyclopterus.lumpus", "Agonus.cataphractus", "Eutrigla.gurnardus",
                                          "Chelidonichthys.cuculus",  "Chelidonichthys.lucerna", "Liparis.liparis",
                                          "Triglopsis.quadricornis",  "Taurulus.bubalis", "Myoxocephalus.quadricornis" ,
                                          "Enchelyopus.cimbrius", "Pollachius.virens","Merlangius.merlangus", "Pollachius.pollachius",
                                          "Trisopterus.luscus", "Merluccius.merluccius", "Melanogrammus.aeglefinus",
                                          "Trisopterus.minutus", "Trisopterus.esmarkii", "Molva.molva",  "Pomatoschistus.minutus", 
                                          "Gobius.niger", "Neogobius.melanostomus", "Aphia.minuta","Gobiusculus.flavescens", 
                                          "Salmo.trutta",  "Coregonus.lavaretus", "Anguilla.anguilla",  "Gasterosteus.aculeatus",
                                          "Spinachia.spinachia", "Pungitius.pungitius", "Engraulis.encrasicolus",   "Alosa.fallax",
                                          "Sardina.pilchardus", "Alosa.agone", "Syngnathus.typhle", "Syngnathus.rostellatus",
                                          "Syngnathus.acus",  "Entelurus.aequoreus", "Nerophis.ophidion", "Lampetra.fluviatilis", 
                                          "Petromyzon.marinus", "Esox.lucius" ))]

fish <- as.data.frame(fish)
fish <- cbind(fish[,(1:2)], apply(fish[,-(1:2)], 2,as.numeric))
fish$year <- as.numeric(fish$year)
###replace NAs with 0
for(i in 1:nrow(fish)){
  for(j in 1:ncol(fish)){
    if(is.na(fish[i,j])){
      fish[i,j] <- 0
    }
  }
}

fish <- fish[which(fish$helcom_id %in% c( "SEA-009" )),]

phyto_zoo_fish <- merge(phyto_zoo, fish, by=c("year","helcom_id"))

phyto_zoo_fish <- cbind(phyto_zoo_fish[,(1:2)], apply(phyto_zoo_fish[,-(1:2)], 2,as.numeric))
##Degree of Novelty analysis ====
## P.S. the stating and ending years can be changed
#the end_year is the last year available in the dataset
phyto_zoo_fish_novelty <- Novelty_function(M=phyto_zoo_fish, start_baseline=1990, end_baseline=1995, end_year=2015, n, method)  #n is the number of years by target bin

#it is one area so novelty and change are the same

#####change point analysis #############

##for novelty######
#choose between:  phyto_novelty, zoo_novelty, fish_novelty, phyto_zoo_novelty, phyto_zoo_fish_novelty
novelty_file <- zoo_novelty #change here!

##change point analysis for novelty ====
novelty_envCpt <- envcpt(novelty_file$cma_dist, minseglen=5)#
plot(novelty_envCpt)
plot(novelty_envCpt, "aic")
plot(novelty_envCpt, "bic")
which.min(BIC(novelty_envCpt))
novelty_envCpt$meancpt@cpts
novelty_envCpt$meanar1cpt@cpts
novelty_envCpt$meanar2cpt@cpts
novelty_envCpt$trendcpt@cpts


##for change#####
#choose between:  phyto_change, zoo_change, fish_change, phyto_zoo_change, phyto_zoo_fish_change
change_file <- phyto_change #change here!

unique(change_file$TS.area) ##choose one of the areas (e.g. SEA-009 , SEA-012 , SEA-013 , SEA-015 , SEA-017)
change_file <- change_file[which(change_file$TS.area == "SEA-015"),] #replace by the area chosen

##change point analysis for novelty ====
change_envCpt <- envcpt(change_file$cma_dist, minseglen=5)#
plot(change_envCpt)
plot(change_envCpt, "aic")
plot(change_envCpt, "bic")
which.min(BIC(change_envCpt))
change_envCpt$meancpt@cpts
change_envCpt$meanar1cpt@cpts
change_envCpt$meanar2cpt@cpts
change_envCpt$trendcpt@cpts



