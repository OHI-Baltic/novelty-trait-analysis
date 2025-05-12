## libraries ----
pkg <- c("here","vegan","analogue","EnvCpt","FD","dplyr","tidyr","ggplot2","zoo")
for(i in 1:length(pkg)){
  if(!(pkg[i] %in% installed.packages())){
    install.packages(pkg[i])
  }
}
for(i in 1:length(pkg)){
  library(pkg[i], character.only = TRUE)
}

## file.paths
data_dir <- here("data")
save_dir <- here("output")


## Yosr's function from COMFORT.R script ----

## M is a dataframe or matrix with colname, in wide format
## with numeric columns for Year, and each trait or species, 
## or other 'feature' to evaluate novelty across; 
## and categorical column for Area subdivisions or gridcells

## baseline start and end should be years used to compare
## end_year is the last year to be compared to the baseline...

## year params could instead refer to binned time periods 
## (but needs to be an unbroken seq of integers, and be same column for baseline and future)

## method (Hellinger or EDLog) is the metric to use to evaluate dissimilarity distance

Novelty_function <- function(M, start_baseline, end_baseline, end_year, method){
  
  ## make sure M is a dataframe
  ## if it's a tibble, some [] subsetting stuff breaks
  M <- as.data.frame(M)
  
  ## columns with features of the space we're 
  ## calculating novelty over
  novelty_feature_cols <- which(!stringr::str_detect(names(M), "Year|Area|Bin"))
  novelty_feature_names <- names(M)[novelty_feature_cols]
  
  ## initialize a dataframe to hold results
  ## cma_dist will be the distance, found using the 'close modern analogues' function
  ## TS.ind are rownames of TS, to keep track how results align with M entries, w/ Years etc
  ## FS.ind are rownames of FS (baseline in comparisons)
  ## TS.area is the area for which novelty is evaluated 
  ## TS.year is the year for which novelty is evaluated
  ## FS.area is the area which has minimum dissimilarity, 
  ## i.e. area or gridcell which is the closest analog
  ## FS.year is the year in the baseline corresponding to the closest analog
  Ana_result <- data.frame(matrix(
    NA, nrow = 0, ncol = 7,
    dimnames = list(NULL, c(
      "cma_dist", "TS.ind", "FS.ind",
      "TS.area", "TS.year", "FS.area", "FS.year"
    ))
  ))
  Ana_result <- Ana_result |> 
    mutate(across(everything(), as.numeric)) |> 
    mutate(across(matches("area"), as.character))
  
  ## data which correspond to the years in the baseline
  ## using FS rownames to keep track of ...
  FS <- M[which(M[,"Year"] %in% c(start_baseline:end_baseline)),]
  FS.a <- FS[,novelty_feature_cols]
  names.fs <- rownames(FS.a)
  n.fs <- length(names.fs)
  
  if(nrow(FS) > 0){
    for(t in (end_baseline+1):end_year){
      ## select data for the year to be compared
      ## (year to calculate novelty for) (???)
      TS <- M[which(M[,"Year"] == t),]
      TS.a <- TS[,novelty_feature_cols]
      names.ts <- rownames(TS.a)
      n.ts <- length(names.ts)
      
      S.a <- rbind(TS.a, FS.a)
      
      if(nrow(TS) > 0){
        if(method == "Hellinger"){
          
          ## CALCULATE HELINGER DISTANCE ----
          ## helinger distance dissimilarity 
          ## based on the composition of the area
          
          ## make a matrix S.a2 where entries
          ## are normalized by the total value in the row
          ## (sum of all features we're calculating novelty across)
          S.a2 <- S.a
          for(i in 1:ncol(S.a)){
            for(j in 1:nrow(S.a)){
              S.a2[j,i] <- S.a[j,i]/rowSums(S.a[j,])
              ## why do we normalize by row sum instead of global max sum??
            }
          }
          S.a2[which(rowSums(TS.a) == FALSE),] <- 0 ## when would rowSums be false??
          ## using distance function from 'analogue' package
          ## 'flexibly calculates distance or dissimilarity measures...'
          S.d2 <- as.matrix(distance(S.a2, method = "chord", dist = T))
          
          S.dd <- S.d2
          S.dd <- matrix(
            ## keep only rows corresponding to TS ie. year to evaluated per loop 
            ## and cols not with dissimilarity corresponding to the TS rows
            ## each row should correspond to one spatial area
            S.dd[ 1 : n.ts, (n.ts+1) : (n.ts+n.fs) ],
            nrow = n.ts,
            ncol = length((n.ts+1) : (n.ts+n.fs))
          )
          colnames(S.dd) <- rownames(FS.a)
          rownames(S.dd) <- rownames(TS.a)
          
        } else if(method == "EDlog"){
          ## Euclidean distance dissimilarity 
          ## based on the stock size of the area
          
          S.a1 <- S.a
          for(i in 1:ncol(S.a)){ 
            ## for to reduce to the log and divide by the number of variables
            ## to avoid the influence of bigger dataset
            for(j in 1:nrow(S.a)){
              S.a1[j,i] <- log(S.a[j,i]+1)
            }
          }
          S.d1 <- as.matrix(dist(S.a1, method = "euclidean"))
          
          S.dd <- S.d1 
          S.dd <- matrix(
            S.dd[ 1 : n.ts, (n.ts+1) : (n.ts+n.fs) ], 
            nrow = n.ts, 
            ncol = length((n.ts+1) : (n.ts+n.fs))
          )
          colnames(S.dd) <- rownames(FS.a)
          rownames(S.dd) <- rownames(TS.a)
        }
        
        ## make a structure 'close' that will contain various attributes, 
        ## used in evaluating closest analogue
        close <- vector("list", n.ts)
        for(i in seq_len(n.ts)){
          if(is.matrix(S.dd)){
            close[[i]] <- sort(S.dd[i,])
          } else {
            close[[i]] <- sort(S.dd[i])
          }
        }
        names(close) <- names.ts 
        
        Fuse.cma <- list(close)
        names(Fuse.cma)[1] <- "close"
        
        ## cma 'close modern analogues' function is from the 'analogue' R package
        .call <- match.call()
        .call[[1]] <- as.name("cma")
        structure(
          list(
            close = close,
            call = .call, 
            quant = "none", 
            probs = "none",
            method = "fuse",
            n.analogs = "all"
          ),
          class = "cma"
        )
        class(Fuse.cma) <- "cma"
        Ana.cma <- Fuse.cma
        
        ## Extract the cma distance for each TS sample
        ## and writes them in data.frame Close.dat
        Close.dat <- data.frame(matrix(
          NA, nrow = 0, ncol = 5,
          dimnames = list(NULL, c(
            "cma_dist", "TS.ind", "FS.ind",
            "TS.area", "TS.year"
          ))
        ))
        TS.ind <- row.names(TS) 
        TS.n <- length(Ana.cma$close)
        
        ## for each Area in TS, make a dataframe with the dissimilarity,
        ## as well as TS and FS indices and areas
        for(i in 1:TS.n){
          k <- TS.ind[i]
          if(length(Ana.cma[["close"]][[k]]) > 0){
            ## first entry of each Ana.cma[["close"]] entry is the smallest distance
            temp <- as.data.frame(Ana.cma[["close"]][[k]][1])
            colnames(temp)[1] <- "cma_dist"
            temp$TS.ind <- as.numeric(TS.ind[i])
            
            if(nrow(FS) == 1){
              temp$FS.ind <- as.numeric(row.names(FS))
            } else {
              temp$FS.ind <- as.numeric(row.names(temp))
            }
            temp$TS.area <- TS[i,"Area"]
            temp$TS.year <- TS[i,"Year"]
            Close.dat <- rbind(Close.dat, temp)
          }
        }
        row.names(Close.dat) <- seq(nrow(Close.dat))
        FS.ind <- row.names(FS)    
        
        for(i in 1:nrow(Close.dat)){
          k <- as.numeric(Close.dat$FS.ind[i])
          Close.dat$FS.area[i] <- as.character(FS[which(rownames(FS) == k), "Area"])
          Close.dat$FS.year[i] <- FS[which(rownames(FS) == k), "Year"]
        }
        
        Ana_result <- rbind(Ana_result, Close.dat)
      }
    }
  }
  
  return(Ana_result)
}

staircase_plot <- function(novelty_results){
  mnYr <- min(novelty_results$Year)
  mxYr <- max(novelty_results$Year)
  
  plot <- ggplot(novelty_results) + 
    geom_tile(aes(x = TS.year, y = Year, fill = cma_dist), color = "black") +
    scale_fill_viridis_c(option = "plasma", direction = -1) +
    labs(x = NULL, y = "Baseline", fill = "Novelty") +
    scale_x_continuous(
      breaks = seq(mnYr, mxYr+4, 2),
      expand = c(0,0),
      sec.axis = dup_axis()
    ) +
    scale_y_continuous(
      breaks = seq(mnYr, mxYr, 2),
      expand = c(0,0),
      sec.axis = dup_axis()
    ) +
    theme_bw() +
    theme(
      legend.key.height = unit(1.5,'cm'),
      legend.key.width = unit(0.3,'cm'),
      axis.text.x.top = element_blank(), 
      axis.text.y.right = element_blank(), 
      axis.title.x.top = element_blank(),
      axis.title.y.right = element_blank()
    ) +
    facet_grid(cols = vars(TS.area))
  
  return(plot)
}

## load and wrangle data ----
cwm_traits_data <- file.path(data_dir, "Multitrophic_CWM_traits_BalticSea.csv")

## for now we just use Gulf of Riga
df0_traits <- read.csv(cwm_traits_data) |> 
  filter(Area == "GoR")

View(df0_traits)

## which years apply for (all) basin(s) being analyzed
useYrs <- 1979:2016
startYrs <- useYrs[5:length(useYrs)-5]
useMinYr <- min(useYrs) - 1


## function to fill NAs using linear trend
fill_trend <- function(x, idx){
  res <- x[idx]
  if(is.na(res)){
    not_na <- which(!is.na(x))
    if(length(not_na) >= 2){
      x <- x[not_na]
      res <- predict(lm(x ~ not_na), data.frame(not_na = idx))[[1]]
    }
  }
  return(res)
}

## four options for fun argument:
## bins5yr will create discrete bins, so data looks like steps with 5yr width
## roll5yr will take rolling mean across all values, filling in NAs in 5yr window
## gfroll5yr replaces (not necessarily all) NAs with a value from rolling mean
## gftrend replaces only (not all) NAs with value from (5yr window) lm-trend

wrangle_data <- function(df0, x, taxa, fun){
  ## CWM option is for traits/ functional novelty analysis
  if(x == "CWM"){
    df <- df0 |> 
      filter(Taxa %in% taxa) |> 
      select(Area, Year, Factor = Trait, Value) |> 
      mutate(across(c("Year", "Value"), as.numeric))
  }
  ## Biomass dataset is surveys (change to pre-processed/modeled data...?)
  if(x == "Biomass"){
    df <- df0 |> 
      filter(Taxa == taxa) |>
      group_by(Area, Year, Month, ScientificName_accepted) |> 
      ## what to do about depths?
      ## at this point we just take per month averages over depths
      summarize(Value = mean(Biomass, na.rm = TRUE)) |> 
      rename(Factor = ScientificName_accepted) 
  }
  
  if(fun == "bins5yr"){
    ## 5-year bin averages
    df_yrs <- df |> 
      filter(Year %in% useYrs) |> 
      mutate(Bin = ceiling((Year-useMinYr)/5)) |> 
      group_by(Bin) |> 
      mutate(maxYear = max(Year)) |> 
      group_by(Bin, Area) |> 
      mutate(nYears = n_distinct(Year))
    
    if("Month" %in% names(df_yrs)){
      ## if multiple months, which to use...
      ## pick the month(s) most often occuring across year-bins??
      # df_yrs <- df_yrs |>
      #   group_by(Area, Factor, Month) |>
      #   mutate(nBins = n_distinct(Bin)) |>
      #   group_by(Area, Factor) |>
      #   mutate(mooMonth = nBins == max(nBins)) |>
      #   ungroup() |>
      #   filter(mooMonth) |>
      #   mutate(Factor = paste(Factor, Month))
      if(taxa %in% c("Phytoplankton", "Zooplankton")){
        df_yrs <- df_yrs |>
          mutate(Factor = case_when(
            Month %in% 3:5 ~ paste(Factor, "- Spring"),
            Month %in% 6:8 ~ paste(Factor, "- Summer"),
            Month %in% 9:11 ~ paste(Factor, "- Autumn")
          ))
      }
    }
    df_long <- df_yrs |> 
      select(Area, Factor, Bin, Year, maxYear, Value) |> 
      arrange(Area, Factor, Bin) |> 
      group_by(Area, Factor, Bin) |> 
      mutate(Value = mean(Value)) |> 
      ungroup()
    FactorLevels <- unique(df_long$Factor)
    df_long$Factor <- factor(df_long$Factor, levels = FactorLevels)
    
    ## expand so have missing years filled with same 5-year bin average
    ## in many cases 5-year bin averages are based on fewer than 5 values
    ## because of missing data
    df_short <- df_long |> 
      distinct(Area, Factor, maxYear, Value) |> 
      rowwise() |> 
      ## bc year list only works if nyears is mod5...
      mutate(maxYear = ifelse(
        maxYear == max(df$Year),
        maxYear + diff(range(df$Year)) %% 5 - 1,
        maxYear
      )) |> 
      mutate(Year = list((maxYear-4):maxYear)) |> 
      ungroup() |> 
      unnest(cols = c(Year)) |> 
      pivot_wider(
        names_from = "Factor", 
        values_from = "Value"
      )
  }
  
  if(fun %in% c("roll5yr", "gfroll5yr", "gftrend")){
    if("Month" %in% names(df)){
      ## if multiple months, which to use...
      ## pick the month(s) most often occuring across all years
      # df <- df |>
      #   group_by(Area, Factor, Month) |>
      #   mutate(nYrs = ifelse(is.na(Month), 0, n_distinct(Year))) |>
      #   group_by(Area, Factor) |>
      #   mutate(mooMonth = nYrs == max(nYrs)) |>
      #   ungroup() |>
      #   filter(mooMonth)
      if(taxa %in% c("Phytoplankton", "Zooplankton")){
        df <- df |>
          mutate(Factor = case_when(
            Month %in% 3:5 ~ paste(Factor, "- Spring"),
            Month %in% 6:8 ~ paste(Factor, "- Summer"),
            Month %in% 9:11 ~ paste(Factor, "- Autumn")
          ))
      }
    }
    FactorLevels <- unique(df$Factor)
    allYrs <- seq(min(df$Year), max(df$Year))
    
    df <- expand.grid(Area=unique(df$Area), Factor=FactorLevels, Year=allYrs) |> 
      left_join(df, by = join_by(Area, Factor, Year)) |>
      group_by(Area, Factor, Year) |> 
      summarize(Value = mean(Value, na.rm = TRUE)) |> 
      arrange(Area, Factor, Year) |> 
      group_by(Area, Factor)
    
    ## 5-year rolling means
    ## CWM data are already rolling means so should only fill with trend???
    if(fun == "roll5yr"){
      df_long <- df |> 
        mutate(Value = rollmean(
          Value, 5, fill = NA,
          align = "right",
          na.rm = TRUE
        )) 
    }
    if(fun == "gfroll5yr"){
      df_long <- df |> 
        mutate(meanFill = rollmean(
          Value, 5, fill = NA,
          align = "right",
          na.rm = TRUE
        )) |> 
        mutate(Value = ifelse(
          is.na(Value), 
          meanFill, Value
        ))
    }
    if(fun == "gftrend"){
      df_long <- df |> 
        mutate(trendFill = rollapply(
          Value, 5, fill = NA,
          FUN = function(x){fill_trend(x, 3)},
          align = "center"
        )) |> 
        ## if still have NAs, add right aligned trend...
        mutate(trendFill2 = rollapply(
          Value, 5, fill = NA,
          FUN = function(x){fill_trend(x, 5)},
          align = "right"
        )) |> 
        mutate(Value = ifelse(
          is.na(Value), 
          trendFill, Value
        )) |> 
        mutate(Value = ifelse(
          is.na(Value),
          trendFill2, Value
        ))
    }
    df_long <- df_long |> 
      filter(Year %in% useYrs) |> 
      select(Area, Factor, Year, Value) |> 
      ## set abundance value rolling means less than zero to zero
      ## (why are there means less than zero??)
      mutate(Value = ifelse(Value > 0, Value, 0)) |> 
      ungroup() |> 
      filter(Year %in% useYrs[5:length(useYrs)]) |> 
      mutate(Bin = ceiling((Year-useMinYr)/5)) |> 
      group_by(Bin) |> 
      mutate(maxYear = max(Year))
    
    df_long$Factor <- factor(df_long$Factor, levels = FactorLevels)
    
    df_short <- df_long |> 
      pivot_wider(names_from = "Factor", values_from = "Value") |> 
      ungroup()
  }
  # View(sort(colSums(is.na(df_short))))
  
  ## plots of the data, to check missing data etc
  p <- ggplot(df_long, aes(x = Year, y = Value, color = Area)) + 
    geom_point(size = 0.5, alpha = 0.5) +
    geom_vline(aes(xintercept = maxYear), linewidth = 0.1) +
    facet_wrap(~Factor, scales = "free_y", nrow = 4) +
    labs(x = "\nYear") +
    theme(legend.position = "bottom")
  
  if(fun == "bin5yr"){
    p <- p + 
      geom_point(shape = 95, size = 4) +
      labs(y = paste(x, "5-Year Bin Averages"))
  }
  if(fun == "gfroll5yr"){
    p <- p + 
      geom_line(linewidth = 0.2) +
      labs(y = paste(x, "Gapfilled 5-Year Rolling Averages"))
  }
  if(fun == "gftrend"){
    p <- p + 
      geom_line(linewidth = 0.2) + 
      labs(y = paste(x, "Gapfilled with 5-Year Window Trend"))
  }
  
  ## return data and plot
  ## long data only used to make plots,
  ## short data used to do the novelty calculation
  return(list(data = df_short, plot = p))
}


## run the function for CWM traits data ----

AllTaxa <- c("Phytoplankton", "Zooplankton", "Benthos", "Fish")

## four options to try for 'fun' argument 
## see 'wrangle-data' code above for details
## bins5yr, roll5yr, gfroll5yr, gftrend
df_traits <- df0_traits |> 
  filter(Variable == "biomass") |> 
  ## can test also with abundance...
  # filter(Variable == "abundance") |> 
  wrangle_data(x = "CWM", taxa = "Phytoplankton", fun = "gftrend")
# wrangle_data(x = "CWM", taxa = "Zooplankton", fun = "gftrend")
# wrangle_data(x = "CWM", taxa = "Benthos", fun = "gftrend")
# wrangle_data(x = "CWM", taxa = "Fish", fun = "gftrend")
# wrangle_data(x = "CWM", taxa = AllTaxa, fun = "gftrend")

View(df_traits$data)
df_traits$plot


## traits/functional novelty calculation ----
traits_novelty_list <- lapply(
  startYrs, FUN = function(x){
    Novelty_function(
      M = df_traits$data,
      start_baseline = x,
      end_baseline = (x+4),
      end_year = max(useYrs),
      method = "Hellinger"
    )
  }
)
names(traits_novelty_list) <- startYrs
traits_novelty <- traits_novelty_list |> 
  bind_rows(.id = "Year") |> 
  mutate(Year = as.numeric(Year))

# View(traits_novelty)
## cma_dist is 
## TS.ind and FS.ind are just rownames for indexing
## TS.area is the area for which novelty is evaluated 
## TS.year is the year for which novelty is evaluated
## FS.area is the area with minimum dissimilarity, i.e. the closest analog
## FS.year is the year in the baseline corresponding to the closest analog


## staircase plot of the novelty
staircase_plot(traits_novelty)

