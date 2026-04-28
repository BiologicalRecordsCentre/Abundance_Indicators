#---
#title: "tidy_functions"
#author: "fiona burns"
#opened: "16.02.2026"
#output: R script
#root.dir: "/data"
#---


#The functions help 
#-   make the datasets consistently formatted to allow combining
#-   process time-series values so they are suitable for inclusion in different indicator models




##### function to remove species with no data all NAs or a comb of NAs and zeros


no_data <- function(data) {
  
  data <- ddply(data, .(scientific_name), mutate,
              any_data = sum(!(is.na(index) | index == 0)))
              
 rm_no_data <- data[!(data$any_data == 0 | data$any_data == 1), ]
  
  rm_no_data <- select(rm_no_data, -any_data)
  
  return(rm_no_data)
}


##### function to remove years with no data 


no_data_yr <- function(data){
  
  if(length(grep("index_sm", names(data))) > 0) {
   data_out <- data %>% 
      group_by(year) %>% 
      mutate(any_data = sum(!(is.na(index)|index == 0)) + 
               sum(!(is.na(index_sm)))) %>% 
     filter(!any_data == 0) %>% 
      select(-any_data)
} else {
    data_out <- data %>% 
    group_by(year) %>% 
      mutate(any_data = sum(!(is.na(index)| index == 0))) %>% 
      filter(!any_data == 0) %>% 
      select(-any_data)
 }
 return(data_out) 
}


##### for datasets of absolute counts (rbbp/scarabbs) remove very rare species

#-   retain species with a median count of at least three (across non-zero values) and max of at least ten


scarce_data <- function(data) {
  
  data$scientific_name <- as.character(data$scientific_name)
  data$index <- as.numeric(data$index)
  data <- ddply(data, .(scientific_name), mutate,
               max = max(index, na.rm = TRUE),
               med = median(index[index > 0], na.rm = TRUE)
             )
  rm_sc_data <- data[!(is.na(data$med)) & data$med >= 3 & data$max >= 10, ]
  
 rm_sc_data <- select(rm_sc_data, -c(med, max))
  
  return(rm_sc_data)
}


##### function to remove species with a year range of less than ten years


rm_dpoor <- function(data) {
#add in yr range column 
  data$scientific_name <- as.character(data$scientific_name)
  data <- ddply(data, .(scientific_name), mutate,
              yr_range = max(year[!is.na(index)], na.rm = TRUE) - min(year[!is.na(index)], na.rm = TRUE))
              
# remove species that have a year range of less than ten 
 data$scientific_name <- as.character(data$scientific_name)
  data_atlst_10 <- data[!data$yr_range < 10, ]
 
  data_atlst_10 <- select(data_atlst_10, -yr_range)
  
  return(data_atlst_10)
}




### TIDY DATA

#-   make sure data are in the correct format to combine in a composite species index
#-   tidying required function produces a df giving the number of species where:
#-   there are NAs at the start/within/at the end of species' time-series
#-   there are zeros at the start/within/at the end of species' time-series
#-   are there missing rows in the time-series

#### TIDYING REQUIRED FUNCTION-----------------------------------------------------


tidying_rqd <- function(data) {
    
nspp <- length(unique(data$scientific_name))
uspp <- sort(unique(data$scientific_name))
nyr  <- max(data$year) - min(data$year) + 1

#are there missing rows in the species time series?
miss_rows <- function(data) {
    #are there missing species year combinations in the data?
miss_sppyrs <- with(data, table(factor((table(scientific_name) - 1) ==
    (max(year) - min(year)), levels = c("FALSE", "TRUE"))))
#are there missing years within a species time series?
miss_sppyrs_in <- 
    with(data, table(factor(table(scientific_name) - 1 == 
    tapply(year, scientific_name, function(x) max(x) - min(x)), 
    levels = c("FALSE", "TRUE"))))
#are there missing years at the end or start of a species time series?
miss_sppyrs_end <- 
   with(data, table(factor(max(year) == tapply(year, scientific_name, max),
                           levels = c("FALSE", "TRUE"))))
miss_sppyrs_str <- 
    with(data, table(factor(min(year) == tapply(year, scientific_name, min), 
                            levels = c("FALSE", "TRUE"))))

# test do all the outputs have two factor levels (TRUE and FALSE)?
if(length(miss_sppyrs) != 2 | length(miss_sppyrs_end) != 2 |
   length(miss_sppyrs_str) != 2 | length(miss_sppyrs_in) != 2) {
    stop("error 1+ tests don't have FALSE and TRUE values")
} else {
miss_rws <- rbind(miss_sppyrs, miss_sppyrs_in, miss_sppyrs_str, miss_sppyrs_end)
} 
return(miss_rws)
}

# are there zero values in the species time series?
zero_v <- function(data) {
# does the number of rows per species equal the number of rows excluding those
    # where index == zero?
zeros <- with(data, table(factor(table(scientific_name) == 
                           table(scientific_name[is.na(index)| index != 0]), 
                                levels = c("FALSE", "TRUE"))))

# zeros at start of time series
zeros_str <- with(data, table(factor( 
     tapply(year[!is.na(index) & index > 0 ], 
                 scientific_name[!is.na(index) & index > 0], min, na.rm = T) == 
    tapply(year[!is.na(index)], scientific_name[!is.na(index)], min, na.rm = T),
   levels = c("FALSE", "TRUE"))))
# zero values within or at the end of a time series
minyr_gzero <- with(data, tapply(year[!is.na(index) & index > 0 ], 
      scientific_name[!is.na(index) & index > 0], min, na.rm = T))
data$minyr_gzero <- minyr_gzero[match(data$scientific_name, names(minyr_gzero))]   

data$scientific_name <- as.factor(data$scientific_name)   
zeros_in <- with(data, table(factor(
    table(scientific_name[year > minyr_gzero & !is.na(index)]) ==
   table(scientific_name[year > minyr_gzero & index != 0 & !is.na(index)]),
     levels = c("FALSE", "TRUE"))))

# test do all the outputs have two factor levels (TRUE and FALSE)?
if(length(zeros) != 2 | length(zeros_str) != 2 |length(zeros_in) != 2) {
    print("error 1+ tests don't have FALSE and TRUE values")
} else {
    zero_values <- rbind(zeros, zeros_str, zeros_in)
} 
return(zero_values)
}

# are there NA values in the species time series?
NA_v <- function(data) {
NA_str <- with(data, table(factor(tapply(year, scientific_name, min) == 
             tapply(year[!is.na(index)], scientific_name[!is.na(index)], min),
             levels = c("FALSE", "TRUE"))))
NA_end <- with(data, table(factor(tapply(year, scientific_name, max) == 
             tapply(year[!is.na(index)], scientific_name[!is.na(index)], max),
                                  levels = c("FALSE", "TRUE"))))

max_sp_yr <- with(data, tapply(year[!is.na(index)],
                                 scientific_name[!is.na(index)], max, na.rm = T))
min_sp_yr <- with(data, tapply(year[!is.na(index)],
                                 scientific_name[!is.na(index)], min, na.rm = T))
data$max_sppyr <- max_sp_yr[match(data$scientific_name, names(max_sp_yr))]
data$min_sppyr <- min_sp_yr[match(data$scientific_name, names(min_sp_yr))]

most_yrs <- max(table(data$scientific_name))
most_yrs <- max(max_sp_yr) - min(min_sp_yr) + 1


NA_in <- with(data, table(factor(table(scientific_name[year >= min_sppyr & year <= max_sppyr ], 
                                       factor(is.na(index[year >= min_sppyr & year <= max_sppyr]), 
                                             levels = c("FALSE", "TRUE")))[, "TRUE"] == 0, 
                                   levels = c("FALSE", "TRUE"))))

# test do all the outputs have two factor levels (TRUE and FALSE)?
    if(length(NA_str) != 2 | length(NA_end) != 2 | length(NA_in) != 2) {
        stop("error 1+ tests don't have FALSE and TRUE values")
    } else {
        NA_values <- rbind(NA_str, NA_end, NA_in)
    } 
return(NA_values)
}    

# join the three reports together with dataset name
tidy_report <- as.data.frame(rbind(miss_rows(data), zero_v(data), NA_v(data)))
tidy_report$test <- row.names(tidy_report)
names(tidy_report) <- c("Fail", "Pass", "test")
tidy_report$survey <- unique(data$.id) 
 
return(tidy_report)
}



#### TIDYING FUNCTION

#-   allows user defined processing of datasets to;
#-   add rows for missing species years
#-   convert zeros at the start of time-series to NAs
#-   where zeros are present within time-series add 1% of the average value of a species' time-series to all elements of it.
#-   use linear interpolation to fill in NAs within species' time-series
#-   hold final time-series values constant where they end before the final year of other species


tidy_d <- function(data, tidy_rpt, exp.g = "N", interp = "N", zero_st = "N", zero_in = "N", exp_tail = "N") {

   nspp   <- length(unique(data$scientific_name))
   uspp   <- sort(unique(data$scientific_name))
   nyr    <- max(data$year) - min(data$year) + 1
   data_t <- data
   data_t$index_orig <- data_t$index
    
  if(tidy_rpt$Fail[tidy_rpt$test == "miss_sppyrs"] > 0 & exp.g == "Y") {
         #expand grid
        
       e_data <- expand.grid(min(data_t$year):max(data_t$year), uspp)
        names(e_data) <- c("year", "scientific_name")
        data_t <- merge(data_t, e_data, by = c("scientific_name", "year"), 
                      sort = T, all = T)
      }

  if(tidy_rpt$Fail[tidy_rpt$test == "zeros_str"] > 0 & zero_st == "Y"|
     tidy_rpt$Fail[tidy_rpt$test == "zeros_in"] > 0  & zero_st == "Y") {
        #change zeros at str to NAs
        
     minyr_gzero <- with(data_t, tapply(year[!is.na(index) & index > 0 ], 
                                   scientific_name[!is.na(index) & index > 0], 
                                     min, na.rm = T))
    data_t$minyr_gzero <- minyr_gzero[match(data_t$scientific_name, 
                                             names(minyr_gzero))]
     data_t$index<- with(data_t, ifelse(index == 0 & year < minyr_gzero,
                                        NA, index))
        
        #test if still any zeros at start of time series
      if (with(data_t, table(factor( 
              tapply(year[!is.na(index) & index > 0 ], 
                scientific_name[!is.na(index) & index > 0], min, na.rm = T) == 
                  tapply(year[!is.na(index)], 
                          scientific_name[!is.na(index)], min, na.rm = T),
              levels = c("FALSE", "TRUE"))))["FALSE"] > 0)   {
          stop("zeros at start of time series")
      } else {  
     
        #add 1% of mean time series to all values for species with zeros within time series
        #get 1% of mean time series per species with zeros
      if(zero_in == "Y") {
          mn_ind <- with(data_t, tapply(index, scientific_name, mean,
                                      na.rm = TRUE)/100)
      for (i in 1:nspp) {
        if(any(data_t$index[data_t$scientific_name== uspp[i]] == 0, na.rm = TRUE)) {    
             data_t$index[data_t$scientific_name == uspp[i]] <-
              data_t$index[data_t$scientific_name == uspp[i]] +
                  mn_ind[match(uspp[i], names(mn_ind))]
          }
        }
       #test if any zeros remain in index
        if(any(data_t$index == 0, na.rm = T) == TRUE) {
           stop("zeros still in index")
             }
           }
        
       data_t <- subset(data_t, select = -c(minyr_gzero)) #remove column
         }
        }
   
       if(tidy_rpt$Fail[tidy_rpt$test =="NA_in"] > 0 & interp == "Y"|
           tidy_rpt$Fail[tidy_rpt$test== "miss_sppyrs_in"] > 0 & interp == "Y") {
        #interpolate NAs within time series, not at start or end, so need to 
        #exclude rows less than minyr per species and greater than max year per species
        
       max_sppyr <- with(data_t, tapply(year[!is.na(index)],
                                scientific_name[!is.na(index)], max, na.rm = T))
       min_sppyr <- with(data_t, tapply(year[!is.na(index)],
                               scientific_name[!is.na(index)], min, na.rm = T))
      data_t$max_sppyr <- max_sppyr[match(data_t$scientific_name, names(max_sppyr))]
       data_t$min_sppyr <- min_sppyr[match(data_t$scientific_name, names(min_sppyr))]
        
    
       interp_ind <- expand.grid(min(data_t$year):max(data_t$year), uspp) 
       names(interp_ind) <- c("year", "scientific_name")
       interp_ind$index  <- 0

       for (i in 1:nspp) {
           interp_ind$index[interp_ind$scientific_name == uspp[i]] <- 
               with(data_t, approx(year[scientific_name == uspp[i]], 
               index[scientific_name == uspp[i]], xout = year[scientific_name == uspp[i]],
                    method = "linear"))$y
            }

       data_t$index <-       
            with(data_t, ifelse(year >= min_sppyr & year <= max_sppyr,
                    interp_ind$index[match(paste(data_t$scientific_name,data_t$year),
                               paste(interp_ind$scientific_name, interp_ind$year))],
                    data_t$index))

        #test to see if any NAs remain within index  
        if(with(data_t, table(factor(
            tapply(index[year<= max_sppyr & year >= min_sppyr],
                   scientific_name[year<= max_sppyr & year >= min_sppyr],
                  function(x) sum(is.na(x)) == 0),
           levels = c("FALSE", "TRUE"))))["FALSE"] != 0)  { 
           stop("still NAs within species time series")
           } 
         
      data_t <- subset(data_t, select = -c(max_sppyr, min_sppyr))
        
        }

        if(tidy_rpt$Fail[tidy_rpt$test == "NA_end"] > 0 & exp_tail == "Y"|
          tidy_rpt$Fail[tidy_rpt$test == "miss_sppyrs_end"] > 0 & exp_tail == "Y") {
        #hold final value const using tail function from BRC_Indicators
        
        fillTailNAs <- function(x){
        
            # Get trues and falses for locations of NAs
           na_true_false <- is.na(x)
            # Get the position of all falses
           na_position <- grep(FALSE, na_true_false)
            # If the max false is hte last year dont do anything...
           if(!max(na_position) == length(x)){
            # else give all the last years the value at the last false
               x[(max(na_position)+1):length(x)] <- x[max(na_position)]
           }
           return(x)    
         }
                
        #apply tail function
       tails <- with(data_t, tapply(index, scientific_name, fillTailNAs))
       data_t$index <- unlist(tails)
        }
    
 
  data_t$common_name <- with(data_t, ifelse(is.na(common_name), 
                                           as.character( data$common_name[match(data_t$scientific_name, data$scientific_name)]),
                                             as.character(data_t$common_name)))
    
return(data_t)
}




##### Function to rescale species time-series relative to first year

rescale <- function(data, ind_col = index, Index_root, ...) {
   data <- 
        data %>% 
        group_by(scientific_name, ...)  %>%  
        dplyr::mutate(min_yr = min(year[!is.na({{ind_col}})]))  %>%
        dplyr::mutate(Index_unscaled = {{ind_col}})  %>%
        dplyr::mutate(across({{ind_col}}, ~.x/.x[year == min_yr] * Index_root)) %>%
        dplyr::mutate(min_yr = NULL)
    return(data)
    }


##### function to calculate geometric mean


geometric <- function (x) exp(mean(log(x), na.rm = TRUE))


##### function to calculate indicator 
#-   Required: 
#-       dataframe of species indices where each row is a species year combination
#-       No missing values within or at the end of species time-series
#-       All species-years represented by a row
#-       No zeroes

#Can estimate the geometric mean using a user specified column with time-series data in (ind_col) and a user specified grouping variable (grp_col). The latter is often grouping by species, but could be by survey or subspecies to get specific values.



indFun1 <- function(ind_data, ind_col, grp_col) {
  
  idata <- ind_data %>% mutate_if(is.factor, as.character) %>% 
        group_by({{grp_col}}) %>% 
    dplyr::mutate(IND = {{ind_col}}) %>% 
    dplyr::mutate(min_yr = min(year[!is.na(IND)])) %>%
       filter(!year < min_yr) 
             
  str_yrs <- sort(as.vector(unique(idata$min_yr)))
  
   if (length(str_yrs) > 1) {
        idata$Iad <- ifelse(idata$min_yr == str_yrs[1], idata$IND, NA) 
        for (i in 2:length(str_yrs)) {
            a <- as.vector(ifelse(idata$min_yr == str_yrs[i],
                                geometric(idata$Iad[idata$year == str_yrs[i]]), NA))
            idata$Iad <- ifelse(idata$min_yr == str_yrs[i],idata$IND / (1 / a), idata$Iad)
        }} 
   if(length(str_yrs) == 1) idata$Iad <- idata$IND
  
  inds <- idata %>% group_by(year) %>% summarise(comp_index = geometric(Iad))
   
  return(inds)
}


