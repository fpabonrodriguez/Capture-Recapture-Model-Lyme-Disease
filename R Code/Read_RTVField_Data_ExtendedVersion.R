#########################################################################
### Felix Pabon-Rodriguez
### Bayesian Capture-Recapture Model for Mice/Tick RTV Field Data
#########################################################################

########################################################################
# Reading data
########################################################################

# RTV data
df2020 <- read.csv(file = "./Data/2020 RTV Mouse Data.csv", header = TRUE)
df2020$Date <- as.character(as.Date(df2020$Date,format = "%m/%d/%y"))
data2020 <- df2020 %>% arrange(as.Date(df2020$Date,format = "%m/%d/%y"))

df2021 <- read.csv(file = "./Data/2021 RTV Mouse Data.csv", header = TRUE)
df2021$Date <- as.character(as.Date(df2021$Date,format = "%m/%d/%Y"))
data2021 <- df2021 %>% arrange(Date)
## Since there is a new definition of recapture variable
## 0=new, 1=recaptured for 2nd-4th time, 2=recaptured for 1st time
## We modified this variable to be 0=new and 1=recapture
data2021$Recapture_Original <- data2021$Recapture
data2021$Recapture <- ifelse(data2021$Recapture==0,0,
                             ifelse(is.na(data2021$Recapture),NA,1))

df2022 <- read.csv(file = "./Data/2022 RTV Mouse Data.csv", header = TRUE)
df2022$Date <- as.character(as.Date(df2022$Date,format = "%m/%d/%Y"))
data2022 <- df2022 %>% arrange(as.Date(df2022$Date))
## Since there is a new definition of recapture variable
## 0=new, 1=recaptured for 2nd-4th time, 2=recaptured for 1st time
## We modified this variable to be 0=new and 1=recapture
data2022$Recapture_Original <- data2022$Recapture
data2022$Recapture <- ifelse(data2022$Recapture==0,0,
                             ifelse(is.na(data2022$Recapture),NA,1))

# PCR data (for infectivity)
pcr_df2020 <- read.csv(file = "./Data/2020 PCR Mouse Data2.csv", header = TRUE)
pcr_df2020$IDfull <- ifelse(pcr_df2020$IDfull==999,NA,pcr_df2020$IDfull)
pcr_df2020$Date <- as.character(as.Date(pcr_df2020$Date,format = "%m/%d/%y"))
pcr_data2020 <- pcr_df2020 %>% arrange(Date)

pcr_df2021 <- read.csv(file = "./Data/2021 PCR Mouse Data2.csv", header = TRUE)
pcr_df2021$UniqueID <- as.numeric(ifelse(pcr_df2021$UniqueID==999,NA,pcr_df2021$UniqueID))
pcr_df2021$Date <- as.Date(pcr_df2021$Date,format = "%m/%d/%Y")
pcr_data2021 <- pcr_df2021 %>% arrange(Date)

pcr_df2022 <- read.csv(file = "./Data/2022 PCR Mouse Data2.csv", header = TRUE)
pcr_df2022$UniqueID <- as.numeric(ifelse(pcr_df2022$UniqueID==999,NA,pcr_df2022$UniqueID))
pcr_df2022$Date <- as.Date(pcr_df2022$Date,format = "%m/%d/%Y")
pcr_data2022 <- pcr_df2022 %>% arrange(Date)


# ELISA data
elisa_df2020 <- read.csv(file = "./Data/2020 ELISA Mouse Data.csv", header = TRUE)
elisa_df2020$UniqueID <- ifelse(elisa_df2020$UniqueID==999,NA,elisa_df2020$UniqueID)
elisa_df2020$Date <- as.character(as.Date(elisa_df2020$Date,format = "%m/%d/%y"))
elisa_data2020 <- elisa_df2020 %>% arrange(Date)

elisa_df2021 <- read.csv(file = "./Data/2021 ELISA Mouse Data.csv", header = TRUE)
elisa_df2021$UniqueID <- ifelse(elisa_df2021$UniqueID==999,NA,elisa_df2021$UniqueID)
elisa_df2021$Date <- as.Date(elisa_df2021$Date,format = "%m/%d/%Y")
elisa_data2021 <- elisa_df2021 %>% arrange(Date)

elisa_df2022 <- read.csv(file = "./Data/2022 ELISA Mouse Data.csv", header = TRUE)
elisa_df2022$UniqueID <- ifelse(elisa_df2022$UniqueID==999,NA,elisa_df2022$UniqueID)
elisa_df2022$Date <- as.Date(elisa_df2022$Date,format = "%m/%d/%Y")
elisa_data2022 <- elisa_df2022 %>% arrange(Date)


# Ticks data
ticks_df2020 <- read.csv(file = "./Data/2020 Tick Data.csv", header = TRUE)
ticks_df2020$UniqueID <- ifelse(ticks_df2020$UniqueID==999,NA,ticks_df2020$UniqueID)
ticks_df2020$Date <- as.character(as.Date(ticks_df2020$Date,format = "%m/%d/%y"))
ticks_data2020 <- ticks_df2020 %>% arrange(Date)

ticks_df2021 <- read.csv(file = "./Data/2021 Tick Data.csv", header = TRUE)
ticks_df2021$UniqueID <- ifelse(ticks_df2021$UniqueID==999,NA,ticks_df2021$UniqueID)
ticks_df2021$Date <- as.Date(ticks_df2021$Date,format = "%m/%d/%Y")
ticks_data2021 <- ticks_df2021 %>% arrange(Date)

ticks_df2022 <- read.csv(file = "./Data/2022 Tick Data.csv", header = TRUE)
ticks_df2022$UniqueID <- ifelse(ticks_df2022$UniqueID==999,NA,ticks_df2022$UniqueID)
ticks_df2022$Date <- as.Date(ticks_df2022$Date,format = "%m/%d/%Y")
ticks_data2022 <- ticks_df2022 %>% arrange(Date)

########################################################################
# Removing Rows with Missing IDs
########################################################################

# 2020
datayear2020 <- data2020[!is.na(data2020$UniqueID) & 
                           !(data2020$UniqueID==2999)&
                           !(data2020$UniqueID==999),]
datayear2020 <- datayear2020 %>% arrange(datayear2020$Date)
datayear2020$Detect <- ifelse((0 %in% datayear2020$Recapture) | 
                                (1 %in% datayear2020$Recapture), 1, 0)

# 2021
datayear2021 <- data2021[!is.na(data2021$UniqueID) & 
                           !(data2021$UniqueID==2999)&
                           !(data2021$UniqueID==999),]
datayear2021 <- datayear2021 %>% arrange(datayear2021$Date)
datayear2021$Detect <- ifelse((0 %in% datayear2021$Recapture) | 
                                (1 %in% datayear2021$Recapture), 1, 0)


#2022
datayear2022 <- data2022[!is.na(data2022$UniqueID) & 
                        !(data2022$UniqueID==2999)&
                        !(data2022$UniqueID==999),]
datayear2022 <- datayear2022 %>% arrange(datayear2022$Date)
datayear2022$Detect <- ifelse((0 %in% datayear2022$Recapture) | 
                             (1 %in% datayear2022$Recapture), 1, 0)


########################################################################
# Correcting trap labels
########################################################################

# Correct label 
# Mouse ID#2219 has a trap label AB1, which was corrected to B1 since 
# all of its encounter occurred on B trap line
location_AB1 <- which(datayear2020$Trap=="AB1")
datayear2020$Trap[location_AB1]<-"B1"

# Insert 0 to string (trap label)
# Example A1 to A01
ch_insert <- function(x, pos, insert) {       
  gsub(paste0("^(.{", pos, "})(.*)$"),
       paste0("\\1", insert, "\\2"),
       x)
}

# Modify trap label
for(i in 1:nrow(datayear2020)){
  datayear2020$Trap[i] <- ifelse(str_length(datayear2020$Trap[i])==3,
                                 datayear2020$Trap[i],
                                 ch_insert(datayear2020$Trap[i],1,"0"))
}

for(i in 1:nrow(datayear2021)){
  datayear2021$Trap[i] <- ifelse(str_length(datayear2021$Trap[i])==3,
                                 datayear2021$Trap[i],
                                 ch_insert(datayear2021$Trap[i],1,"0"))
}

for(i in 1:nrow(datayear2022)){
  datayear2022$Trap[i] <- ifelse(str_length(datayear2022$Trap[i])==3,
                                 datayear2022$Trap[i],
                             ch_insert(datayear2022$Trap[i],1,"0"))
}
 
########################################################################
# 2020 and 2021 general info
########################################################################

unique_dates_2020 <- unique(datayear2020$Date)
unique_number_dates_2020 <- length(unique(datayear2020$Date))
unique_dates_2021 <- unique(datayear2021$Date)
unique_number_dates_2021 <- length(unique(datayear2021$Date))
unique_dates_2022 <- unique(datayear2022$Date)
unique_number_dates_2022 <- length(unique(datayear2022$Date))

########################################################################
# Array of encounters per year
#           occ_1                          occ_2  ...        occ_K
#   ind_1   trap where mouse was trapped  
#   ...
#   ind_M 
########################################################################

## adding rows for data augmentation process, to get total of M
bigM.2020 <- 2800
bigM.2021 <- bigM.2022 <- 2500
add_rows2020 <- bigM.2020 - length(unique(datayear2020$UniqueID))
add_rows2021 <- bigM.2021 - length(unique(datayear2021$UniqueID))
add_rows2022 <- bigM.2022 - length(unique(datayear2022$UniqueID))

## z indicators (if belongs to population)
z_ind_2020 <- c(rep(1,length(unique(datayear2020$UniqueID))), rep(NA,add_rows2020))
z_ind_2021 <- c(rep(1,length(unique(datayear2021$UniqueID))), rep(NA,add_rows2021))
z_ind_2022 <- c(rep(1,length(unique(datayear2022$UniqueID))), rep(NA,add_rows2022))

# Sorted traps (these will be used as categories of where the mouse was trapped)
sorted_traps2020 <- str_sort(unique(datayear2020$Trap))
sorted_traps2021 <- str_sort(unique(datayear2021$Trap))
sorted_traps2022 <- str_sort(unique(datayear2022$Trap))
selected_variables <- c("Date","UniqueID","Trap","Detect")

# Sorted sites
sorted_sites <- str_sort(unique(datayear2020$Site))
sorted_sites2022 <- c("HHB",str_sort(unique(datayear2022$Site))[-c(5,7)])

# Records of selected variables
record2020 <- datayear2020[,selected_variables]
record2021 <- datayear2021[,selected_variables]
record2022 <- datayear2022[,selected_variables]

# Adding trap labels as categories (2020)
for(i in 1:nrow(record2020)){
  record2020$TrapCateg[i] <- match(record2020$Trap[i],sorted_traps2020,
                                   nomatch = length(sorted_traps2020)+1)
}

# Adding trap labels as categories (2021)
for(i in 1:nrow(record2021)){
  record2021$TrapCateg[i] <- match(record2021$Trap[i],sorted_traps2021,
                                   nomatch = length(sorted_traps2021)+1)
}

# Adding trap labels as categories (2022)
for(i in 1:nrow(record2022)){
  record2022$TrapCateg[i] <- match(record2022$Trap[i],sorted_traps2022,
                                   nomatch = length(sorted_traps2022)+1)
}

# Getting Capture History (function)
get_CaptHist <- function(input.data){
  input.data %>%
    distinct() %>%
    spread(Date, TrapCateg, fill = 0) %>% 
    group_by(UniqueID)
}
 
# Getting the 3D array - year 2020
CapHist2020 <- get_CaptHist(record2020)
tempHist2020 <- CapHist2020[,!names(CapHist2020) %in% c("Trap", "Detect")]

EncounterHist2020 <- tempHist2020 %>%
  group_by(UniqueID) %>%
  summarise_all(sum)
EncounterHist2020[EncounterHist2020==0] <- length(sorted_traps2020)+1

temp_EncHist2020 <- as.matrix(EncounterHist2020[,-1])
added_individual2020 <- matrix(data = length(sorted_traps2020)+1, 
                               add_rows2020,ncol(temp_EncHist2020))
EncHist2020 <- rbind(temp_EncHist2020,added_individual2020)

# Getting the 3D array - year 2021
CapHist2021 <- get_CaptHist(record2021)
tempHist2021 <- CapHist2021[,!names(CapHist2021) %in% c("Trap", "Detect")]

EncounterHist2021 <- tempHist2021 %>%
  group_by(UniqueID) %>%
  summarise_all(sum) 
EncounterHist2021[EncounterHist2021==0] <- length(sorted_traps2021)+1

temp_EncHist2021 <- as.matrix(EncounterHist2021[,-1])
added_individual2021 <- matrix(data = length(sorted_traps2021)+1, 
                               add_rows2021,ncol(temp_EncHist2021))
EncHist2021 <- rbind(temp_EncHist2021,added_individual2021)


# Getting the 3D array - year 2022
CapHist2022 <- get_CaptHist(record2022)
tempHist2022 <- CapHist2022[,!names(CapHist2022) %in% c("Trap", "Detect")]

EncounterHist2022 <- tempHist2022 %>%
  group_by(UniqueID) %>%
  summarise_all(sum) 
EncounterHist2022[EncounterHist2022==0] <- length(sorted_traps2022)+1

temp_EncHist2022 <- as.matrix(EncounterHist2022[,-1])
added_individual2022 <- matrix(data = length(sorted_traps2022)+1, 
                               add_rows2022,ncol(temp_EncHist2022))
EncHist2022 <- rbind(temp_EncHist2022,added_individual2022)


########################################################################
# Number of Encounters by Traps
########################################################################

cls2020 <- paste(1:length(sorted_traps2020))
cls2021 <- paste(1:length(sorted_traps2021))
cls2022 <- paste(1:length(sorted_traps2022))

# 2020
ct_2020_1 <- record2020 %>%
  distinct() %>%
  group_by(Date,UniqueID) %>%
  count(TrapCateg, sort = TRUE) %>%
  arrange(UniqueID)

ct_2020_2 <- ct_2020_1 %>%
  group_by(UniqueID,TrapCateg) %>%
  summarise(across(n, 
                   function(x) ifelse(all(is.na(x)),
                                      NA,
                                      sum(x,na.rm = TRUE)))) %>%
  pivot_wider(names_from = TrapCateg, values_from = n, values_fill = 0) 

CapTrap_2020 <- rbind(unname(as.matrix(ct_2020_2[,cls2020])),
             matrix(data = 0,add_rows2020,length(cls2020)))


# 2021
ct_2021_1 <- record2021 %>%
  distinct() %>%
  group_by(Date,UniqueID) %>%
  count(TrapCateg, sort = TRUE) %>%
  arrange(UniqueID)

ct_2021_2 <- ct_2021_1 %>%
  group_by(UniqueID,TrapCateg) %>%
  summarise(across(n, 
                   function(x) ifelse(all(is.na(x)),
                                      NA,
                                      sum(x,na.rm = TRUE)))) %>%
  pivot_wider(names_from = TrapCateg, values_from = n, values_fill = 0) 

CapTrap_2021 <- rbind(unname(as.matrix(ct_2021_2[,cls2021])),
                      matrix(data = 0,add_rows2021,length(cls2021)))


# 2022
ct_2022_1 <- record2022 %>%
  distinct() %>%
  group_by(Date,UniqueID) %>%
  count(TrapCateg, sort = TRUE) %>%
  arrange(UniqueID)

ct_2022_2 <- ct_2022_1 %>%
  group_by(UniqueID,TrapCateg) %>%
  summarise(across(n, 
                   function(x) ifelse(all(is.na(x)),
                                      NA,
                                      sum(x,na.rm = TRUE)))) %>%
  pivot_wider(names_from = TrapCateg, values_from = n, values_fill = 0) 

CapTrap_2022 <- rbind(unname(as.matrix(ct_2022_2[,cls2022])),
                      matrix(data = 0,add_rows2022,length(cls2022)))


########################################################################
# Behavior - Definition 2
# Behavior[i,j] = 1 if subj i is trapped more than once by trap j
#               = 0 if was not captured or capture only once
# First definition is given later
########################################################################

# 2020
tb2020 <- unname(as.matrix(ct_2020_2[,cls2020]))
tb2020 <- ifelse(tb2020 > 1, 1, 0)
Trap_Behavior2020 <- rbind(tb2020,
                           matrix(data = NA,add_rows2020,length(cls2020)))


# 2021
tb2021 <- unname(as.matrix(ct_2021_2[,cls2021]))
tb2021 <- ifelse(tb2021 > 1, 1, 0)
Trap_Behavior2021 <- rbind(tb2021,
                           matrix(data = NA,add_rows2021,length(cls2021)))

# 2022
tb2022 <- unname(as.matrix(ct_2022_2[,cls2022]))
tb2022 <- ifelse(tb2022 > 1, 1, 0)
Trap_Behavior2022 <- rbind(tb2022,
                           matrix(data = NA,add_rows2022,length(cls2022)))

########################################################################
# Tick (Infection) Data
########################################################################

# 2020 data
# adding Infected_Ticks indicator
ticks_data2020$Infected_Ticks <- rep(NA,nrow(ticks_data2020))
for(i in 1:nrow(ticks_data2020)){
  ticks_data2020$Infected_Ticks[i] <- ifelse(ticks_data2020$flaB_Ct[i] <= 37, 1, 0)
}

# Modify trap label
for(i in 1:nrow(ticks_data2020)){
  ticks_data2020$Trap[i] <- ifelse(str_length(ticks_data2020$Trap[i])==3,
                                   ticks_data2020$Trap[i],
                                   ch_insert(ticks_data2020$Trap[i],1,"0"))
}

# get subset of data (ticks extracted from mouse)
ticks_data2020_mouse <- ticks_data2020[ticks_data2020$From_Mouse==1,]
ticks_data2020_drag <- ticks_data2020[ticks_data2020$From_Drag==1,]

# Trying to get some IDs back from RTV field data
id_for_merge_2020 <- datayear2020[,c("Date","Site","Trap","UniqueID","IDfull")]
tmp4 <- left_join(ticks_data2020_mouse, id_for_merge_2020, by = c("Date","Site"))
for(i in 1:nrow(tmp4)){
  if(!(is.na(tmp4$Trap.x[i])) & 
     !(is.na(tmp4$Trap.y[i])) & 
     (tmp4$Trap.y[i]==tmp4$Trap.x[i])){
    tmp4$UniqueID.x[i] <- tmp4$UniqueID.y[i]
  }
}

for(i in 1:nrow(ticks_data2020_mouse)){
  if(!is.na(ticks_data2020_mouse$UniqueID[i])) next
  someIDs4 <- tmp4[(tmp4$Lab.assigned.number==i),] 
  keep4 <- someIDs4[someIDs4$UniqueID.x==someIDs4$UniqueID.y,]
  ticks_data2020_mouse$UniqueID[i] <- ifelse(NA %in% keep4$UniqueID.x,NA,keep4$UniqueID.x)
}

ticks_data2020_mouse <- ticks_data2020_mouse %>% arrange(as.Date(ticks_data2020_mouse$Date))
ticks_data2020_mouse <- ticks_data2020_mouse[!is.na(ticks_data2020_mouse$UniqueID),]

tempInfectedTicks2020 <- ticks_data2020_mouse[,c("UniqueID","Date","Infected_Ticks")] %>%
  group_by(UniqueID) %>%
  mutate(row = row_number()) %>%
  pivot_wider(names_from = Date, values_from = Infected_Ticks)


tempInfectedTicks2020_2 <- (tempInfectedTicks2020 %>% 
                              group_by(UniqueID) %>% 
                              summarise(across(everything(), 
                                               function(x) ifelse(all(is.na(x)),
                                                                  NA,
                                                                  sum(x,na.rm = TRUE)))) %>%
                              as.matrix)[,-2]

# Get same rows as in previous data components
id2020_table <- (left_join(EncounterHist2020[,1],
                           datayear2020[,1:5],
                           by="UniqueID") %>%
                   group_by(UniqueID) %>%
                   filter(row_number() == 1))[,1]


add_dates2020 <- setdiff(unique_dates_2020,colnames(tempInfectedTicks2020_2)[-1])
tempInfectedTicks2020_3 <- as.matrix(left_join(id2020_table,
                                               tempInfectedTicks2020_2,
                                               by="UniqueID", copy = TRUE)[,-1]) 
tempInfectedTicks2020_4 <- cbind(tempInfectedTicks2020_3,
                                 matrix(NA,nrow(tempInfectedTicks2020_3),length(add_dates2020)))            
colnames(tempInfectedTicks2020_4)[ncol(tempInfectedTicks2020_4):(ncol(tempInfectedTicks2020_4)-length(add_dates2020)+1)] <- add_dates2020                                    
InfectedTicks_Mouse2020 <- tempInfectedTicks2020_4[,unique_dates_2020]


tdrag1 <- ticks_data2020_drag[!is.na(ticks_data2020_drag$Date),c("Date","Site","Infected_Ticks")] %>%
  arrange(Site) %>%
  group_by(Site)

tmp_InfectedTicks_Drag2020 <- as.data.frame((ticks_data2020_drag[!is.na(ticks_data2020_drag$Date),c("Date","Site","Infected_Ticks")] %>%
                           group_by(Site) %>%
                           mutate(row = row_number()) %>%
                           pivot_wider(names_from = Date, values_from = Infected_Ticks) %>%
                           summarise(across(everything(), 
                                            function(x) ifelse(all(is.na(x)),
                                                               NA,
                                                               sum(x,na.rm = TRUE)))))[,-c(1,2)])
InfectedTicks_Drag2020 <- as.matrix(tmp_InfectedTicks_Drag2020)


# 2021 data
# adding Infected_Ticks indicator
ticks_data2021$Infected_Ticks <- rep(NA,nrow(ticks_data2021))
for(i in 1:nrow(ticks_data2021)){
  ticks_data2021$Infected_Ticks[i] <- ifelse(ticks_data2021$flaB_Ct[i] <= 37, 1, 0)
}


# Modify trap label
for(i in 1:nrow(ticks_data2021)){
  ticks_data2021$Trap[i] <- ifelse(str_length(ticks_data2021$Trap[i])==3,
                                   ticks_data2021$Trap[i],
                                   ch_insert(ticks_data2021$Trap[i],1,"0"))
}

# get subset of data (ticks extracted from mouse)
ticks_data2021_mouse <- ticks_data2021[ticks_data2021$From_Mouse==1,]
ticks_data2021_drag <- ticks_data2021[ticks_data2021$From_Drag==1,]

# Trying to get some IDs back from RTV field data
id_for_merge_2021 <- datayear2021[,c("Date","Site","Trap","UniqueID")]
id_for_merge_2021$Date <- as.Date(id_for_merge_2021$Date, format = "%m/%d/%Y")
tmp5 <- left_join(ticks_data2021_mouse, id_for_merge_2021, by = c("Date","Site"))

for(i in 1:nrow(ticks_data2021_mouse)){
  if(!is.na(ticks_data2021_mouse$UniqueID[i])) next
  someIDs5 <- tmp5[(tmp4$Lab.assigned.number==i),] 
  keep5 <- someIDs5[someIDs5$UniqueID.x==someIDs5$UniqueID.y,]
  ticks_data2021_mouse$UniqueID[i] <- ifelse(NA %in% keep5$UniqueID.x,NA,keep5$UniqueID.x)
}

ticks_data2021_mouse <- ticks_data2021_mouse %>% arrange(as.Date(ticks_data2021_mouse$Date))
ticks_data2021_mouse <- ticks_data2021_mouse[!is.na(ticks_data2021_mouse$UniqueID),]

tempInfectedTicks2021 <- ticks_data2021_mouse[,c("UniqueID","Date","Infected_Ticks")] %>%
  group_by(UniqueID) %>%
  mutate(row = row_number()) %>%
  pivot_wider(names_from = Date, values_from = Infected_Ticks)

tempInfectedTicks2021_2 <- (tempInfectedTicks2021 %>% 
                              group_by(UniqueID) %>% 
                              summarise(across(everything(), 
                                               function(x) ifelse(all(is.na(x)),
                                                                  NA,
                                                                  sum(x,na.rm = TRUE)))) %>%
                              as.matrix)[,-2]


# Get same rows as in previous data components
id2021_table <- (left_join(EncounterHist2021[,1],
                           datayear2021[,1:6],
                           by="UniqueID") %>%
                   group_by(UniqueID) %>%
                   filter(row_number() == 1))[,1]

add_dates2021 <- setdiff(unique_dates_2021,colnames(tempInfectedTicks2021_2)[-1])
tempInfectedTicks2021_3 <- as.matrix(left_join(id2021_table,
                                               tempInfectedTicks2021_2,
                                               by="UniqueID", copy = TRUE)[,-1]) 
tempInfectedTicks2021_4 <- cbind(tempInfectedTicks2021_3,
                                 matrix(NA,nrow(tempInfectedTicks2021_3),length(add_dates2021)))            
colnames(tempInfectedTicks2021_4)[ncol(tempInfectedTicks2021_4):(ncol(tempInfectedTicks2021_4)-length(add_dates2021)+1)] <- add_dates2021                                    
InfectedTicks_Mouse2021 <- tempInfectedTicks2021_4[,unique_dates_2021]

tdrag4 <- ticks_data2021_drag[!is.na(ticks_data2021_drag$Date),c("Date","Site","Infected_Ticks")] %>%
  arrange(Site) %>%
  group_by(Site)

tmp_InfectedTicks_Drag2021 <- as.data.frame((ticks_data2021_drag[!is.na(ticks_data2021_drag$Date),c("Date","Site","Infected_Ticks")] %>%
                           group_by(Site) %>%
                           mutate(row = row_number()) %>%
                           pivot_wider(names_from = Date, values_from = Infected_Ticks) %>%
                           summarise(across(everything(), 
                                            function(x) ifelse(all(is.na(x)),
                                                               NA,
                                                               sum(x,na.rm = TRUE)))))[,-c(1,2)])
InfectedTicks_Drag2021 <- as.matrix(tmp_InfectedTicks_Drag2021)


# 2022 data
# adding Infected_Ticks indicator
ticks_data2022$Infected_Ticks <- rep(NA,nrow(ticks_data2022))
for(i in 1:nrow(ticks_data2022)){
  ticks_data2022$Infected_Ticks[i] <- ifelse(ticks_data2022$flaB_Ct[i] <= 37, 1, 0)
}


# Modify trap label
for(i in 1:nrow(ticks_data2022)){
  ticks_data2022$Trap[i] <- ifelse(str_length(ticks_data2022$Trap[i])==3,
                                   ticks_data2022$Trap[i],
                                   ch_insert(ticks_data2022$Trap[i],1,"0"))
}

# get subset of data (ticks extracted from mouse)
ticks_data2022_mouse <- ticks_data2022[ticks_data2022$From_Mouse==1,]
ticks_data2022_drag <- ticks_data2022[ticks_data2022$From_Drag==1,]

# Trying to get some IDs back from RTV field data
id_for_merge_2022 <- datayear2022[,c("Date","Site","Trap","UniqueID")]
id_for_merge_2022$Date <- as.Date(id_for_merge_2022$Date, format = "%m/%d/%Y")
tmp5 <- left_join(ticks_data2022_mouse, id_for_merge_2022, by = c("Date","Site"))

for(i in 1:nrow(ticks_data2022_mouse)){
  if(!is.na(ticks_data2022_mouse$UniqueID[i])) next
  someIDs5 <- tmp5[(tmp4$Lab.assigned.number==i),] 
  keep5 <- someIDs5[someIDs5$UniqueID.x==someIDs5$UniqueID.y,]
  ticks_data2022_mouse$UniqueID[i] <- ifelse(NA %in% keep5$UniqueID.x,NA,keep5$UniqueID.x)
}

ticks_data2022_mouse <- ticks_data2022_mouse %>% arrange(as.Date(ticks_data2022_mouse$Date))
ticks_data2022_mouse <- ticks_data2022_mouse[!is.na(ticks_data2022_mouse$UniqueID),]

tempInfectedTicks2022 <- ticks_data2022_mouse[,c("UniqueID","Date","Infected_Ticks")] %>%
  group_by(UniqueID) %>%
  mutate(row = row_number()) %>%
  pivot_wider(names_from = Date, values_from = Infected_Ticks)

tempInfectedTicks2022_2 <- (tempInfectedTicks2022 %>% 
                              group_by(UniqueID) %>% 
                              summarise(across(everything(), 
                                               function(x) ifelse(all(is.na(x)),
                                                                  NA,
                                                                  sum(x,na.rm = TRUE)))) %>%
                              as.matrix)[,-2]


# Get same rows as in previous data components
id2022_table <- (left_join(EncounterHist2022[,1],
                           datayear2022[,1:6],
                           by="UniqueID") %>%
                   group_by(UniqueID) %>%
                   filter(row_number() == 1))[,1]

add_dates2022 <- setdiff(unique_dates_2022,colnames(tempInfectedTicks2022_2)[-1])
tempInfectedTicks2022_3 <- as.matrix(left_join(id2022_table,
                                               tempInfectedTicks2022_2,
                                               by="UniqueID", copy = TRUE)[,-1]) 
tempInfectedTicks2022_4 <- cbind(tempInfectedTicks2022_3,
                                 matrix(NA,nrow(tempInfectedTicks2022_3),length(add_dates2022)))            
colnames(tempInfectedTicks2022_4)[ncol(tempInfectedTicks2022_4):(ncol(tempInfectedTicks2022_4)-length(add_dates2022)+1)] <- add_dates2022                                    
InfectedTicks_Mouse2022 <- tempInfectedTicks2022_4[,unique_dates_2022]

tdrag4 <- ticks_data2022_drag[!is.na(ticks_data2022_drag$Date),c("Date","Site","Infected_Ticks")] %>%
  arrange(Site) %>%
  group_by(Site)

tmp_InfectedTicks_Drag2022 <- as.data.frame((ticks_data2022_drag[!is.na(ticks_data2022_drag$Date),c("Date","Site","Infected_Ticks")] %>%
                                               group_by(Site) %>%
                                               mutate(row = row_number()) %>%
                                               pivot_wider(names_from = Date, values_from = Infected_Ticks) %>%
                                               summarise(across(everything(), 
                                                                function(x) ifelse(all(is.na(x)),
                                                                                   NA,
                                                                                   sum(x,na.rm = TRUE)))))[,-c(1,2)])
InfectedTicks_Drag2022 <- as.matrix(tmp_InfectedTicks_Drag2022)


########################################################################
# Number of Dragged Ticks per Site and Time
########################################################################

# 2020
tmp_DraggedTicks_2020 <- as.data.frame((ticks_data2020_drag[!is.na(ticks_data2020_drag$Date),c("Date","Site","From_Drag")] %>%
                 group_by(Site) %>%
                 mutate(row = row_number()) %>%
                 pivot_wider(names_from = Date, values_from = From_Drag) %>%
                 summarise(across(everything(), 
                                  function(x) ifelse(all(is.na(x)),
                                                     NA,
                                                     sum(x,na.rm = TRUE)))))[,-c(1,2)])
DraggedTicks_2020 <- as.matrix(tmp_DraggedTicks_2020)

# 2021
tmp_DraggedTicks_2021 <- as.data.frame((ticks_data2021_drag[!is.na(ticks_data2021_drag$Date),c("Date","Site","From_Drag")] %>%
                                      group_by(Site) %>%
                                      mutate(row = row_number()) %>%
                                      pivot_wider(names_from = Date, values_from = From_Drag) %>%
                                      summarise(across(everything(), 
                                                       function(x) ifelse(all(is.na(x)),
                                                                          NA,
                                                                          sum(x,na.rm = TRUE)))))[,-c(1,2)])
DraggedTicks_2021 <- as.matrix(tmp_DraggedTicks_2021)


# 2022
tmp_DraggedTicks_2022 <- as.data.frame((ticks_data2022_drag[!is.na(ticks_data2022_drag$Date),c("Date","Site","From_Drag")] %>%
                                          group_by(Site) %>%
                                          mutate(row = row_number()) %>%
                                          pivot_wider(names_from = Date, values_from = From_Drag) %>%
                                          summarise(across(everything(), 
                                                           function(x) ifelse(all(is.na(x)),
                                                                              NA,
                                                                              sum(x,na.rm = TRUE)))))[,-c(1,2)])
DraggedTicks_2022 <- as.matrix(tmp_DraggedTicks_2022)


########################################################################
# Vector of Times (Use union of times between processes)
########################################################################

#2020
all_times_2020 <- unique(sort(c(colnames(InfectedTicks_Mouse2020),
                               as.vector(na.omit(unique(ticks_data2020_drag$Date))))))

full_dates_2020 <- seq.Date(as.Date(all_times_2020[1]),
                            as.Date(all_times_2020[length(all_times_2020)]),
                            "day")

time_captures_2020 <- unique_dates_2020
time_drags_2020 <- as.vector(na.omit(unique(ticks_data2020_drag$Date)))


#2021
all_times_2021 <- unique(sort(c(colnames(InfectedTicks_Mouse2021),
                               unique(format(ticks_data2021_drag$Date, "%Y-%m-%d")))))

full_dates_2021 <- seq.Date(as.Date(all_times_2021[1]),
                            as.Date(all_times_2021[length(all_times_2021)]),
                            "day")

time_captures_2021 <- unique_dates_2021
time_drags_2021 <- unique(format(ticks_data2021_drag$Date, "%Y-%m-%d"))


#2022
all_times_2022 <- unique(sort(unique(format(ticks_data2022$Date, "%Y-%m-%d"))))

full_dates_2022 <- seq.Date(as.Date(all_times_2022[1]),
                            as.Date(all_times_2022[length(all_times_2022)]),
                            "day")

time_captures_2022 <- unique_dates_2022
time_drags_2022 <- unique(format(ticks_data2022$Date[ticks_data2022$From_Drag==1], "%Y-%m-%d"))

# all
full_dates_2020 <- format(full_dates_2020, "%Y-%m-%d")
full_dates_2021 <- format(full_dates_2021, "%Y-%m-%d")
full_dates_2022 <- format(full_dates_2022, "%Y-%m-%d")
 

########################################################################
# Sites membership
########################################################################

# 2020
sitetable2020 <- left_join(EncounterHist2020[,1],
                           datayear2020[,1:4],
                           by="UniqueID") %>%
  group_by(UniqueID) %>%
  filter(row_number() == 1)

# Adding site labels as categories (2020)
for(i in 1:nrow(sitetable2020)){
  sitetable2020$SiteCateg[i] <- match(sitetable2020$Site[i],sorted_sites)
}
Site2020 <- c(sitetable2020$SiteCateg,rep(NA,add_rows2020))


# 2021
sitetable2021 <- left_join(EncounterHist2021[,1],
                           datayear2021[,1:6],
                           by="UniqueID") %>%
  group_by(UniqueID) %>%
  filter(row_number() == 1)

# Adding site labels as categories (2021)
for(i in 1:nrow(sitetable2021)){
  sitetable2021$SiteCateg[i] <- match(sitetable2021$Site[i],sorted_sites)
}
Site2021 <- c(sitetable2021$SiteCateg,rep(NA,add_rows2021))


# 2022
sitetable2022 <- left_join(EncounterHist2022[,1],
                           datayear2022[,1:6],
                           by="UniqueID") %>%
  group_by(UniqueID) %>%
  filter(row_number() == 1)

# Adding site labels as categories (2022)
for(i in 1:nrow(sitetable2022)){
  sitetable2022$SiteCateg[i] <- match(sitetable2022$Site[i],sorted_sites2022)
}
Site2022 <- c(sitetable2022$SiteCateg,rep(NA,add_rows2022))


########################################################################
# Other information
########################################################################

M.2020 <- nrow(EncHist2020)
M.2021 <- nrow(EncHist2021)
M.2022 <- nrow(EncHist2022)
nsites2020 <- nsites2021 <- nsites2022 <- 6
ntraps2020 <- length(sorted_traps2020)  
ntraps2021 <- length(sorted_traps2021)
ntraps2022 <- length(sorted_traps2022)

########################################################################
# Construct final version of previous objects based on the union dates
# between drag and capture event dates
########################################################################

# 2020
# Encounter
Final_EncHist2020_tmp <- cbind(EncHist2020,
                          matrix(NA,nrow(EncHist2020),
                                 length(setdiff(full_dates_2020,time_captures_2020))))
colnames(Final_EncHist2020_tmp) <- c(time_captures_2020,setdiff(full_dates_2020,time_captures_2020))

Final_EncHist2020 <- Final_EncHist2020_tmp[,full_dates_2020]

# Infected Ticks (Mice)
Final_InfTicks_Mice2020_tmp <- cbind(InfectedTicks_Mouse2020,
                                     matrix(NA,nrow(InfectedTicks_Mouse2020),
                                            length(setdiff(full_dates_2020,time_captures_2020))))
colnames(Final_InfTicks_Mice2020_tmp) <- c(time_captures_2020,setdiff(full_dates_2020,time_captures_2020))

Final_InfTicks_Mice2020 <- Final_InfTicks_Mice2020_tmp[,full_dates_2020]

# Infected Ticks (Drag)
Final_InfTicks_Drag2020_tmp <- cbind(InfectedTicks_Drag2020,
                                     matrix(NA,nrow(InfectedTicks_Drag2020),
                                            length(setdiff(full_dates_2020,time_drags_2020))))
colnames(Final_InfTicks_Drag2020_tmp) <- c(time_drags_2020,setdiff(full_dates_2020,time_drags_2020))

Final_InfTicks_Drag2020 <- Final_InfTicks_Drag2020_tmp[,full_dates_2020]

# 2021
# Encounter
Final_EncHist2021_tmp <- cbind(EncHist2021,
                               matrix(NA,nrow(EncHist2021),
                                      length(setdiff(full_dates_2021,time_captures_2021))))
colnames(Final_EncHist2021_tmp) <- c(time_captures_2021,setdiff(full_dates_2021,time_captures_2021))

Final_EncHist2021 <- Final_EncHist2021_tmp[,full_dates_2021]

# Infected Ticks (Mice)
Final_InfTicks_Mice2021_tmp <- cbind(InfectedTicks_Mouse2021,
                                     matrix(NA,nrow(InfectedTicks_Mouse2021),
                                            length(setdiff(full_dates_2021,time_captures_2021))))
colnames(Final_InfTicks_Mice2021_tmp) <- c(time_captures_2021,setdiff(full_dates_2021,time_captures_2021))

Final_InfTicks_Mice2021 <- Final_InfTicks_Mice2021_tmp[,full_dates_2021]
Final_InfTicks_Mice2021[16,16] <- 1 # fixing discrepancy

# Infected Ticks (Drag)
Final_InfTicks_Drag2021_tmp <- cbind(InfectedTicks_Drag2021,
                                     matrix(NA,nrow(InfectedTicks_Drag2021),
                                            length(setdiff(full_dates_2021,time_drags_2021))))
colnames(Final_InfTicks_Drag2021_tmp) <- c(time_drags_2021,setdiff(full_dates_2021,time_drags_2021))

Final_InfTicks_Drag2021 <- Final_InfTicks_Drag2021_tmp[,full_dates_2021]


# 2022
# Encounter
Final_EncHist2022_tmp <- cbind(EncHist2022,
                               matrix(NA,nrow(EncHist2022),
                                      length(setdiff(full_dates_2022,time_captures_2022))))
colnames(Final_EncHist2022_tmp) <- c(time_captures_2022,
                                            setdiff(full_dates_2022,time_captures_2022))
Final_EncHist2022 <- Final_EncHist2022_tmp[,full_dates_2022]
Final_EncHist2022[117,71] <- which(sorted_traps2022=="B06")

# Infected Ticks (Mice)
Final_InfTicks_Mice2022_tmp <- cbind(InfectedTicks_Mouse2022,
                                     matrix(NA,nrow(InfectedTicks_Mouse2022),
                                            length(setdiff(full_dates_2022,time_captures_2022))))
colnames(Final_InfTicks_Mice2022_tmp) <- c(time_captures_2022,setdiff(full_dates_2022,time_captures_2022))

Final_InfTicks_Mice2022 <- Final_InfTicks_Mice2022_tmp[,full_dates_2022]

# Infected Ticks (Drag)
Final_InfTicks_Drag2022_tmp <- cbind(InfectedTicks_Drag2022,
                                     matrix(NA,nrow(InfectedTicks_Drag2022),
                                            length(setdiff(full_dates_2022,time_drags_2022))))
colnames(Final_InfTicks_Drag2022_tmp) <- c(time_drags_2022,setdiff(full_dates_2022,time_drags_2022))

Final_InfTicks_Drag2022 <- Final_InfTicks_Drag2022_tmp[,full_dates_2022]

########################################################################
# Binary Encounter
########################################################################

Binary_Encounter2020 <- ifelse(Final_EncHist2020==ntraps2020+1,0,1)
Binary_Encounter2021 <- ifelse(Final_EncHist2021==ntraps2021+1,0,1)
Binary_Encounter2022 <- ifelse(Final_EncHist2022==ntraps2022+1,0,1)

########################################################################
# Last capture occasion for each
########################################################################

lp2020 <- ifelse(Final_EncHist2020==(ntraps2020+1),0,1)
lp_temp2020 <- unlist(apply(lp2020,1,function(x) last(which(x==1))))
lastcapture2020 <- ifelse(is.na(lp_temp2020),ncol(Final_EncHist2020),lp_temp2020)

lp2021 <- ifelse(Final_EncHist2021==(ntraps2021+1),0,1)
lp_temp2021 <- unlist(apply(lp2021,1,function(x) last(which(x==1))))
lastcapture2021 <- ifelse(is.na(lp_temp2021),ncol(Final_EncHist2021),lp_temp2021)

lp2022 <- ifelse(Final_EncHist2022==(ntraps2022+1),0,1)
lp_temp2022 <- unlist(apply(lp2022,1,function(x) last(which(x==1))))
lastcapture2022 <- ifelse(is.na(lp_temp2022),ncol(Final_EncHist2022),lp_temp2022)

########################################################################
# First capture occasion for each
# need this for next section of code
########################################################################

ft2020 <- ifelse(Final_EncHist2020==(ntraps2020+1),0,1)
ft_temp2020 <- unlist(apply(ft2020,1,function(x) first(which(x==1))))
firstcapture2020 <- ifelse(is.na(ft_temp2020),ncol(Final_EncHist2020),ft_temp2020)

ft2021 <- ifelse(Final_EncHist2021==(ntraps2021+1),0,1)
ft_temp2021 <- unlist(apply(ft2021,1,function(x) first(which(x==1))))
firstcapture2021 <- ifelse(is.na(ft_temp2021),ncol(Final_EncHist2021),ft_temp2021)

ft2022 <- ifelse(Final_EncHist2022==(ntraps2022+1),0,1)
ft_temp2022 <- unlist(apply(ft2022,1,function(x) first(which(x==1))))
firstcapture2022 <- ifelse(is.na(ft_temp2022),ncol(Final_EncHist2022),ft_temp2022)

########################################################################
# Covariates
########################################################################

# AgeAdult (binary)
#2020
datayear2020$AgeAdult <- ifelse(datayear2020$Age == 1, 1, 
                                ifelse(is.na(datayear2020$Age),NA,0))

AgeAdult2020 <- pull((datayear2020[,c("UniqueID","Date","AgeAdult")] %>%
                        arrange(UniqueID) %>% 
                        group_by(UniqueID) %>%
                        filter(row_number(AgeAdult) == 1))[,"AgeAdult"])
AgeAdult2020 <- c(AgeAdult2020,rep(NA,add_rows2020))

#2021
datayear2021$AgeAdult <- ifelse(datayear2021$Age == 1, 1, 
                                ifelse(is.na(datayear2021$Age),NA,0))

AgeAdult2021 <- pull((datayear2021[,c("UniqueID","Date","AgeAdult")] %>%
                        arrange(UniqueID) %>% 
                        group_by(UniqueID) %>%
                        filter(row_number(AgeAdult) == 1))[,"AgeAdult"])
AgeAdult2021 <- c(AgeAdult2021,rep(NA,add_rows2021))

#2022
datayear2022$AgeAdult <- ifelse(datayear2022$Age == 1, 1, 
                                ifelse(is.na(datayear2022$Age),NA,0))

AgeAdult2022 <- pull((datayear2022[,c("UniqueID","Date","AgeAdult")] %>%
                        arrange(UniqueID) %>% 
                        group_by(UniqueID) %>%
                        filter(row_number(AgeAdult) == 1))[,"AgeAdult"])
AgeAdult2022 <- c(AgeAdult2022,rep(NA,add_rows2022))



# AgeSubAdult (binary)
datayear2020$AgeSubAdult <- ifelse(datayear2020$Age == 2, 1, 
                                   ifelse(is.na(datayear2020$Age),NA,0))
AgeSubAdult2020 <- pull((datayear2020[,c("UniqueID","Date","AgeSubAdult")] %>%
                           arrange(UniqueID) %>% 
                           group_by(UniqueID) %>%
                           filter(row_number(AgeSubAdult) == 1))[,"AgeSubAdult"])
AgeSubAdult2020 <- c(AgeSubAdult2020,rep(NA,add_rows2020))


datayear2021$AgeSubAdult <- ifelse(datayear2021$Age == 2, 1, 
                                   ifelse(is.na(datayear2021$Age),NA,0))
AgeSubAdult2021 <- pull((datayear2021[,c("UniqueID","Date","AgeSubAdult")] %>%
                           arrange(UniqueID) %>% 
                           group_by(UniqueID) %>%
                           filter(row_number(AgeSubAdult) == 1))[,"AgeSubAdult"])
AgeSubAdult2021 <- c(AgeSubAdult2021,rep(NA,add_rows2021))


datayear2022$AgeSubAdult <- ifelse(datayear2022$Age == 2, 1, 
                                   ifelse(is.na(datayear2022$Age),NA,0))
AgeSubAdult2022 <- pull((datayear2022[,c("UniqueID","Date","AgeSubAdult")] %>%
                           arrange(UniqueID) %>% 
                           group_by(UniqueID) %>%
                           filter(row_number(AgeSubAdult) == 1))[,"AgeSubAdult"])
AgeSubAdult2022 <- c(AgeSubAdult2022,rep(NA,add_rows2022))


# SexMale (binary)
datayear2020$SexMale <- ifelse(datayear2020$Sex == 1, 1, 
                               ifelse(is.na(datayear2020$Sex),NA,0))
SexMale2020 <- pull((datayear2020[,c("UniqueID","Date","SexMale")] %>%
                       arrange(UniqueID) %>% 
                       group_by(UniqueID) %>%
                       filter(row_number(SexMale) == 1))[,"SexMale"])
SexMale2020 <- c(SexMale2020,rep(NA,add_rows2020))


datayear2021$SexMale <- ifelse(datayear2021$Sex == 1, 1, 
                               ifelse(is.na(datayear2021$Sex),NA,0))
SexMale2021 <- pull((datayear2021[,c("UniqueID","Date","SexMale")] %>%
                       arrange(UniqueID) %>% 
                       group_by(UniqueID) %>%
                       filter(row_number(SexMale) == 1))[,"SexMale"])
SexMale2021 <- c(SexMale2021,rep(NA,add_rows2021))


datayear2022$SexMale <- ifelse(datayear2022$Sex == 1, 1, 
                               ifelse(is.na(datayear2022$Sex),NA,0))
SexMale2022 <- pull((datayear2022[,c("UniqueID","Date","SexMale")] %>%
                       arrange(UniqueID) %>% 
                       group_by(UniqueID) %>%
                       filter(row_number(SexMale) == 1))[,"SexMale"])
SexMale2022 <- c(SexMale2022,rep(NA,add_rows2022))


# Behavior
Behavior2020 <- matrix(data = 0, 
                       nrow(temp_EncHist2020),
                       ncol(Final_EncHist2020))

for(i in 1:nrow(temp_EncHist2020)){
  Behavior2020[i,1:firstcapture2020[i]] <- rep(1,length(1:firstcapture2020[i]))
}

Behavior2021 <- matrix(data = 0, 
                       nrow(temp_EncHist2021),
                       ncol(Final_EncHist2021))

for(i in 1:nrow(temp_EncHist2021)){
  Behavior2021[i,1:firstcapture2021[i]] <- rep(1,length(1:firstcapture2021[i]))
}


Behavior2022 <- matrix(data = 0, 
                       nrow(temp_EncHist2022),
                       ncol(Final_EncHist2022))

for(i in 1:nrow(temp_EncHist2022)){
  Behavior2022[i,1:firstcapture2022[i]] <- rep(1,length(1:firstcapture2022[i]))
}


Final_Behavior2020 <- rbind(Behavior2020, matrix(1,add_rows2020,ncol = ncol(Behavior2020)))
Final_Behavior2021 <- rbind(Behavior2021, matrix(1,add_rows2021,ncol = ncol(Behavior2021)))
Final_Behavior2022 <- rbind(Behavior2022, matrix(1,add_rows2022,ncol = ncol(Behavior2022)))

########################################################################
# Dead Indicator
########################################################################

# 2020
for(i in 1:nrow(datayear2020)){
  string <- c("dead", "die", "died", "deceased","DEAD", "DIE", "DIED", "DECEASED")
  str_in <- sapply(string, function(x){grepl(x, datayear2020$Comments[i])})
  datayear2020$dead[i] <- ifelse(sum(str_in) > 0,1,0)
}

temp_dead2020 <- as.matrix(datayear2020[,c("Date","UniqueID","dead")] %>%
                             distinct() %>%
                             spread(Date, dead, fill = 0) %>% 
                             group_by(UniqueID))[,-1]


Dead2020_tmp <- cbind(temp_dead2020,
                            matrix(0,nrow(temp_dead2020),
                                   length(setdiff(full_dates_2020,time_captures_2020))))
colnames(Dead2020_tmp) <- c(time_captures_2020,setdiff(full_dates_2020,time_captures_2020))
Dead2020 <- Dead2020_tmp[,full_dates_2020]

# have dead after first indicator is 1
for(i in 1:nrow(Dead2020)){
  look_at_nonzero <- first(which(as.vector(Dead2020[i,]) != 0))
  if(is.na(look_at_nonzero)){
    next
  } else{
    Dead2020[i,look_at_nonzero:ncol(Dead2020)] <- rep(1,length(look_at_nonzero:ncol(Dead2020)))
  }
}


# 2021
for(i in 1:nrow(datayear2021)){
  string <- c("dead", "die", "died", "deceased","DEAD", "DIE", "DIED", "DECEASED")
  str_in <- sapply(string, function(x){grepl(x, datayear2021$Notes[i])})
  datayear2021$dead[i] <- ifelse(sum(str_in) > 0,1,0)
}

temp_dead2021 <- as.matrix(datayear2021[,c("Date","UniqueID","dead")] %>%
                             distinct() %>%
                             spread(Date, dead, fill = 0) %>% 
                             group_by(UniqueID))[,-1]


Dead2021_tmp <- cbind(temp_dead2021,
                      matrix(0,nrow(temp_dead2021),
                             length(setdiff(full_dates_2021,time_captures_2021))))
colnames(Dead2021_tmp) <- c(time_captures_2021,setdiff(full_dates_2021,time_captures_2021))
Dead2021 <- Dead2021_tmp[,full_dates_2021]

# have dead after first indicator is 1
for(i in 1:nrow(Dead2021)){
  look_at_nonzero <- first(which(as.vector(Dead2021[i,]) != 0))
  if(is.na(look_at_nonzero)){
    next
  } else{
    Dead2021[i,look_at_nonzero:ncol(Dead2021)] <- rep(1,length(look_at_nonzero:ncol(Dead2021)))
  }
}


# 2022
for(i in 1:nrow(datayear2022)){
  string <- c("dead", "die", "died", "deceased","DEAD", "DIE", "DIED", "DECEASED")
  str_in <- sapply(string, function(x){grepl(x, datayear2022$Comments[i])})
  datayear2022$dead[i] <- ifelse(sum(str_in) > 0,1,0)
}

temp_dead2022 <- as.matrix(datayear2022[,c("Date","UniqueID","dead")] %>%
                             distinct() %>%
                             spread(Date, dead, fill = 0) %>% 
                             group_by(UniqueID))[,-1]


Dead2022_tmp <- cbind(temp_dead2022,
                      matrix(0,nrow(temp_dead2022),
                             length(setdiff(full_dates_2022,time_captures_2022))))
colnames(Dead2022_tmp) <- c(time_captures_2022,
                                   setdiff(full_dates_2022,time_captures_2022))
Dead2022 <- Dead2022_tmp[,full_dates_2022]

# have dead after first indicator is 1
for(i in 1:nrow(Dead2022)){
  look_at_nonzero <- first(which(as.vector(Dead2022[i,]) != 0))
  if(is.na(look_at_nonzero)){
    next
  } else{
    Dead2022[i,look_at_nonzero:ncol(Dead2022)] <- rep(1,length(look_at_nonzero:ncol(Dead2022)))
  }
}


########################################################################
# (Tested) Ticks Count
########################################################################

# 2020
ticks_data2020_mouse$Count <- rep(1,nrow(ticks_data2020_mouse))
temp_TicksCount2020_1 <- ticks_data2020_mouse[,c("Date","UniqueID","Count")] %>% 
    arrange(Date) %>% 
    group_by(UniqueID,Date) %>% 
    summarise(across(everything(), 
                     function(x) ifelse(all(is.na(x)),
                                        NA,
                                        sum(x,na.rm = TRUE)))) %>%
  pivot_wider(names_from = Date, values_from = Count)

temp_TicksCount2020_2 <- as.matrix(left_join(id2020_table,
                                             temp_TicksCount2020_1,
                                               by="UniqueID", copy = TRUE)[,-1])

Final_Ticks_Mice2020_tmp <- cbind(temp_TicksCount2020_2,
                                     matrix(NA,nrow(temp_TicksCount2020_2),
                                            length(setdiff(full_dates_2020,unique(ticks_data2020_mouse$Date)))))
colnames(Final_Ticks_Mice2020_tmp) <- c(unique(ticks_data2020_mouse$Date),
                                        setdiff(full_dates_2020,unique(ticks_data2020_mouse$Date)))
Final_Ticks_Mice2020 <- Final_Ticks_Mice2020_tmp[,full_dates_2020]


# 2021
ticks_data2021_mouse$Count <- rep(1,nrow(ticks_data2021_mouse))
temp_TicksCount2021_1 <- ticks_data2021_mouse[,c("Date","UniqueID","Count")] %>% 
  arrange(Date) %>% 
  group_by(UniqueID,Date) %>% 
  summarise(across(everything(), 
                   function(x) ifelse(all(is.na(x)),
                                      NA,
                                      sum(x,na.rm = TRUE)))) %>%
  pivot_wider(names_from = Date, values_from = Count)

temp_TicksCount2021_2 <- as.matrix(left_join(id2021_table,
                                             temp_TicksCount2021_1,
                                             by="UniqueID", copy = TRUE)[,-1])

Final_Ticks_Mice2021_tmp <- cbind(temp_TicksCount2021_2,
                                  matrix(NA,nrow(temp_TicksCount2021_2),
                                         length(setdiff(full_dates_2021,
                                                        unique(format(ticks_data2021_mouse$Date, "%Y-%m-%d"))))))

colnames(Final_Ticks_Mice2021_tmp) <- c(unique(format(ticks_data2021_mouse$Date, "%Y-%m-%d")),
                                        setdiff(full_dates_2021,
                                                unique(format(ticks_data2021_mouse$Date, "%Y-%m-%d"))))
Final_Ticks_Mice2021 <- Final_Ticks_Mice2021_tmp[,full_dates_2021]


# 2022
ticks_data2022_mouse$Count <- rep(1,nrow(ticks_data2022_mouse))
temp_TicksCount2022_1 <- ticks_data2022_mouse[,c("Date","UniqueID","Count")] %>% 
  arrange(Date) %>% 
  group_by(UniqueID,Date) %>% 
  summarise(across(everything(), 
                   function(x) ifelse(all(is.na(x)),
                                      NA,
                                      sum(x,na.rm = TRUE)))) %>%
  pivot_wider(names_from = Date, values_from = Count)

temp_TicksCount2022_2 <- as.matrix(left_join(id2022_table,
                                             temp_TicksCount2022_1,
                                             by="UniqueID", copy = TRUE)[,-1])

Final_Ticks_Mice2022_tmp <- cbind(temp_TicksCount2022_2,
                              matrix(NA,nrow(temp_TicksCount2022_2),
                                     length(setdiff(full_dates_2022,
                                                    unique(format(ticks_data2022_mouse$Date, "%Y-%m-%d"))))))

colnames(Final_Ticks_Mice2022_tmp) <- c(unique(format(ticks_data2022_mouse$Date, "%Y-%m-%d")),
                                    setdiff(full_dates_2022,
                                            unique(format(ticks_data2022_mouse$Date, "%Y-%m-%d"))))
Final_Ticks_Mice2022 <- Final_Ticks_Mice2022_tmp[,full_dates_2022]


########################################################################
# Dragged Ticks Counts with Common Dates 
########################################################################

# 2020
tmp_DraggedTicks2020 <- cbind(DraggedTicks_2020,
                              matrix(NA,nrow(DraggedTicks_2020),
                                     length(setdiff(full_dates_2020,time_drags_2020))))
colnames(tmp_DraggedTicks2020) <- c(time_drags_2020,setdiff(full_dates_2020,time_drags_2020))
Final_DraggedTicks2020 <- tmp_DraggedTicks2020[,full_dates_2020]

# 2021
tmp_DraggedTicks2021 <- cbind(DraggedTicks_2021,
                              matrix(NA,nrow(DraggedTicks_2021),
                                     length(setdiff(full_dates_2021,time_drags_2021))))
colnames(tmp_DraggedTicks2021) <- c(time_drags_2021,setdiff(full_dates_2021,time_drags_2021))
Final_DraggedTicks2021 <- tmp_DraggedTicks2021[,full_dates_2021]

# 2022
tmp_DraggedTicks2022 <- cbind(DraggedTicks_2022,
                              matrix(NA,nrow(DraggedTicks_2022),
                                     length(setdiff(full_dates_2022,time_drags_2022))))
colnames(tmp_DraggedTicks2022) <- c(time_drags_2022,setdiff(full_dates_2022,time_drags_2022))
Final_DraggedTicks2022 <- tmp_DraggedTicks2022[,full_dates_2022]


########################################################################
# Percent of Vaccine Eaten by Mice
########################################################################

# PercentAte (numeric)
# 2020
datayear2020$Count <- 1:nrow(datayear2020)
tempPercentAte2020 <- (datayear2020[,c("Date","UniqueID","Count","PercentAte")] %>%
                         group_by(UniqueID) %>%
                         pivot_wider(names_from = Date, values_from = PercentAte))[,-2]

PercentAte2020 <- (tempPercentAte2020 %>% 
                     group_by(UniqueID) %>% 
                     summarise(across(everything(), 
                                      function(x) ifelse(all(is.na(x)),
                                                         NA,
                                                         sum(x,na.rm = TRUE)))) %>%
                     as.matrix)[,-1]

# There are two integers in the dataset that need to be change
# data node PercentAte[18, 24] is 2
# data node PercentAte[54, 10] is 2
# since we are not sure if the correct measure is 0.20 or 0.02, or a 1 (typo)
# changed it to a NA.
PercentAte2020 <- ifelse(PercentAte2020 > 1, NA, PercentAte2020)


# 2021
datayear2021$Count <- 1:nrow(datayear2021)
tempPercentAte2021 <- (datayear2021[,c("Date","UniqueID","Count","PercentAte")] %>%
                         group_by(UniqueID) %>%
                         pivot_wider(names_from = Date, values_from = PercentAte))[,-2]

PercentAte2021 <- (tempPercentAte2021 %>% 
                     group_by(UniqueID) %>% 
                     summarise(across(everything(), 
                                      function(x) ifelse(all(is.na(x)),
                                                         NA,
                                                         sum(x,na.rm = TRUE)))) %>%
                     as.matrix)[,-1]

# There are several integers in the dataset that need to be changed
# sum(1*(PercentAte2021 > 1),na.rm = TRUE)  =>>  17
# integers like 2, 5, 50, and 75 are present; changed to NA
PercentAte2021 <- ifelse(PercentAte2021 > 1, NA, PercentAte2021)


# 2022
datayear2022$Count <- 1:nrow(datayear2022)
tempPercentAte2022 <- (datayear2022[,c("Date","UniqueID","Count","PercentAte")] %>%
                         group_by(UniqueID) %>%
                         pivot_wider(names_from = Date, values_from = PercentAte))[,-2]

PercentAte2022 <- (tempPercentAte2022 %>% 
                     group_by(UniqueID) %>% 
                     summarise(across(everything(), 
                                      function(x) ifelse(all(is.na(x)),
                                                         NA,
                                                         sum(x,na.rm = TRUE)))) %>%
                     as.matrix)[,-1]
PercentAte2022 <- ifelse(PercentAte2022 > 1, NA, PercentAte2022)


# Final Versions (with all dates)

# 2020
tmp_PercentAte2020 <- cbind(PercentAte2020,
                       matrix(NA,nrow(PercentAte2020),
                              length(setdiff(full_dates_2020,time_captures_2020))))
colnames(tmp_PercentAte2020) <- c(time_captures_2020,setdiff(full_dates_2020,time_captures_2020))
Final_PercentAte2020 <- tmp_PercentAte2020[,full_dates_2020]

# 2021
tmp_PercentAte2021 <- cbind(PercentAte2021,
                            matrix(NA,nrow(PercentAte2021),
                                   length(setdiff(full_dates_2021,time_captures_2021))))
colnames(tmp_PercentAte2021) <- c(time_captures_2021,setdiff(full_dates_2021,time_captures_2021))
Final_PercentAte2021 <- tmp_PercentAte2021[,full_dates_2021]

# 2022
tmp_PercentAte2022 <- cbind(PercentAte2022,
                            matrix(NA,nrow(PercentAte2022),
                                   length(setdiff(full_dates_2022,time_captures_2022))))
colnames(tmp_PercentAte2022) <- c(time_captures_2022,setdiff(full_dates_2022,time_captures_2022))
Final_PercentAte2022 <- tmp_PercentAte2022[,full_dates_2022]


########################################################################
# PCR data (MICE)
########################################################################

get_BinomialObject_perTime <- function(data){
  
  number_dates <-  ncol(data)
  number_ID <- nrow(data)
  success <- matrix(NA, nrow = number_ID, ncol = number_dates)
  success <- ifelse(data==1,1,success) 
  
  ft_temp <- c()
  for(i in 1:number_ID){
    ft_temp[i] <- ifelse(is.integer(first(which(data[i,]==1))) && 
                           length(first(which(data[i,]==1))) == 0L,NA,
                         first(which(data[i,]==1)))
  }
  
  for(i in 1:number_ID){
    if(is.na(ft_temp[i])){
      success[i,] <- success[i,]
    } else{
      success[i,(ft_temp[i]):number_dates] <- rep(1,length((ft_temp[i]):number_dates))
    }
  }
  
  for(i in 1:number_ID){
    for(j in 1:number_dates){
      if(is.na(data[i,j])){
        next
      } else if(data[i,j]==0){
        success[i,j] <- 0
      }
    }
  } 
  
  return(success)
}


# 2020 
# adding infection indicator

# old criteria
#pcr_data2020$Infected <- rep(NA,nrow(pcr_data2020))
#for(i in 1:nrow(pcr_data2020)){
#  if(i %in% 1:128){
#    pcr_data2020$Infected[i] <- ifelse(pcr_data2020$qPCR[i] >= 37.5, 1, 0)
#  } else{
#    pcr_data2020$Infected[i] <- ifelse(pcr_data2020$qPCR[i] <= 36.5, 1, 0)
#  }
#}

# pcr_data2020$Infected <- ifelse(pcr_data2020$qPCR >= 100, 1, 0)
pcr_data2020$Infected <- ifelse(pcr_data2020$qPCR >= 10, 1, 0)


# Modify trap label
for(i in 1:nrow(pcr_data2020)){
  pcr_data2020$Trap[i] <- ifelse(str_length(str_conv(datayear2020$Trap[i],"UTF-8"))==3,
                                 pcr_data2020$Trap[i],
                                 ch_insert(pcr_data2020$Trap[i],1,"0"))
}


# Trying to get some IDs back from RTV field data
id_for_merge_2020 <- datayear2020[,c("Date","Site","Trap","UniqueID","IDfull")]
tmp <- left_join(pcr_data2020, id_for_merge_2020, by = c("Date","Site"))
for(i in 1:nrow(tmp)){
  if(!(is.na(tmp$Trap.x[i])) & 
     !(is.na(tmp$Trap.y[i])) & 
     (tmp$Trap.y[i]==tmp$Trap.x[i])){
    tmp$IDfull.x[i] <- tmp$IDfull.y[i]
  }
}

for(i in 1:nrow(pcr_data2020)){
  if(!is.na(pcr_data2020$IDfull[i])) next
  someIDs <- tmp[(tmp$UT.ear.number==i),] 
  keep <- someIDs[someIDs$IDfull.x==someIDs$IDfull.y,]
  pcr_data2020$IDfull[i] <- ifelse(NA %in% keep$IDfull.x,NA,keep$IDfull.x)
}

pcr_data2020 <- pcr_data2020 %>% arrange(as.Date(pcr_data2020$Date))
pcr_data2020 <- pcr_data2020[!is.na(pcr_data2020$IDfull),]
pcr_data2020$UniqueID <- 2000 + pcr_data2020$IDfull

tempInfected2020 <- pcr_data2020[,c("UniqueID","Date","Infected")] %>%
  distinct() %>%
  spread(Date, Infected, fill = NA) %>% 
  group_by(UniqueID)

tempPCR2020 <- pcr_data2020[,c("UniqueID","Date","qPCR")] %>%
  distinct() %>%
  spread(Date, qPCR, fill = NA) %>% 
  group_by(UniqueID)

# Get same rows as in previous data components
id2020_table <- (left_join(EncounterHist2020[,1],
                           datayear2020[,1:5],
                           by="UniqueID") %>%
                   group_by(UniqueID) %>%
                   filter(row_number() == 1))[,1]

Infected2020 <- as.matrix(left_join(id2020_table,
                                              tempInfected2020,
                                              by="UniqueID")[,-1])

PCR2020 <- as.matrix(left_join(id2020_table,
                                     tempPCR2020,
                                          by="UniqueID")[,-1])


# 2021 
# adding infection indicator

# old criteria
#pcr_data2021$Infected <- rep(NA,nrow(pcr_data2021))
#for(i in 1:nrow(pcr_data2021)){
#    pcr_data2021$Infected[i] <- ifelse(pcr_data2021$qPCR[i] <= 37, 1, 0)
#}

# pcr_data2021$Infected <- ifelse(pcr_data2021$qPCR >= 100, 1, 0)
pcr_data2021$Infected <- ifelse(pcr_data2021$qPCR >= 10, 1, 0)

# Modify trap label
for(i in 1:nrow(pcr_data2021)){
  pcr_data2021$Trap[i] <- ifelse(str_length(str_conv(datayear2021$Trap[i],"UTF-8"))==3,
                                 pcr_data2021$Trap[i],
                                 ch_insert(pcr_data2021$Trap[i],1,"0"))
}

# Trying to get some IDs back from RTV field data
id_for_merge_2021 <- datayear2021[,c("Date","Site","Trap","UniqueID")]
id_for_merge_2021$Date <- as.Date(id_for_merge_2021$Date, format = "%m/%d/%Y")
tmp2 <- left_join(pcr_data2021, id_for_merge_2021, by = c("Date","Site"))
for(i in 1:nrow(tmp2)){
  if(!(is.na(tmp2$Trap.x[i])) & 
     !(is.na(tmp2$Trap.y[i])) & 
     (tmp2$Trap.y[i]==tmp2$Trap.x[i])){
    tmp2$UniqueID.x[i] <- tmp2$UniqueID.y[i]
  }
}

for(i in 1:nrow(pcr_data2021)){
  if(!is.na(pcr_data2021$UniqueID[i])) next
  someIDs2 <- tmp2[(tmp2$UT.ear.number==i),] 
  keep2 <- someIDs2[someIDs2$UniqueID.x==someIDs2$UniqueID.y,]
  pcr_data2021$UniqueID[i] <- ifelse(NA %in% keep2$UniqueID.x,NA,keep2$UniqueID.x)
}

pcr_data2021 <- pcr_data2021 %>% arrange(as.Date(pcr_data2021$Date))
pcr_data2021 <- pcr_data2021[!is.na(pcr_data2021$UniqueID),]

tempInfected2021 <- pcr_data2021[,c("UniqueID","Date","Infected")] %>%
  distinct() %>%
  spread(Date, Infected, fill = NA) %>% 
  group_by(UniqueID)

tempPCR2021 <- pcr_data2021[,c("UniqueID","Date","qPCR")] %>%
  distinct() %>%
  spread(Date, qPCR, fill = NA) %>% 
  group_by(UniqueID)

# Get same rows as in previous data components
id2021_table <- (left_join(EncounterHist2021[,1],
                           datayear2021[,1:6],
                           by="UniqueID") %>%
                   group_by(UniqueID) %>%
                   filter(row_number() == 1))[,1]

Infected2021 <- as.matrix(left_join(id2021_table,
                                          tempInfected2021,
                                          by="UniqueID")[,-1])[,unique_dates_2021]
PCR2021 <- as.matrix(left_join(id2021_table,
                                     tempPCR2021,
                                     by="UniqueID")[,-1])[,unique_dates_2021]


# 2022 
# adding infection indicator

# old criteria 
#pcr_data2022$Infected <- rep(NA,nrow(pcr_data2022))
#for(i in 1:nrow(pcr_data2022)){
#  pcr_data2022$Infected[i] <- ifelse(pcr_data2022$qPCR[i] <= 37, 1, 0)
#}

# pcr_data2022$Infected <- ifelse(pcr_data2022$qPCR >= 100, 1, 0)
pcr_data2022$Infected <- ifelse(pcr_data2022$qPCR >= 10, 1, 0)

# Modify trap label
for(i in 1:nrow(pcr_data2022)){
  pcr_data2022$Trap[i] <- ifelse(str_length(str_conv(datayear2022$Trap[i],"UTF-8"))==3,
                                 pcr_data2022$Trap[i],
                                 ch_insert(pcr_data2022$Trap[i],1,"0"))
}

# Trying to get some IDs back from RTV field data
id_for_merge_2022 <- datayear2022[,c("Date","Site","Trap","UniqueID")]
id_for_merge_2022$Date <- as.Date(id_for_merge_2022$Date, format = "%m/%d/%Y")
tmp2 <- left_join(pcr_data2022, id_for_merge_2022, by = c("Date","Site"))
for(i in 1:nrow(tmp2)){
  if(!(is.na(tmp2$Trap.x[i])) & 
     !(is.na(tmp2$Trap.y[i])) & 
     (tmp2$Trap.y[i]==tmp2$Trap.x[i])){
    tmp2$UniqueID.x[i] <- tmp2$UniqueID.y[i]
  }
}

for(i in 1:nrow(pcr_data2022)){
  if(!is.na(pcr_data2022$UniqueID[i])) next
  someIDs2 <- tmp2[(tmp2$UT.ear.number==i),] 
  keep2 <- someIDs2[someIDs2$UniqueID.x==someIDs2$UniqueID.y,]
  pcr_data2022$UniqueID[i] <- ifelse(NA %in% keep2$UniqueID.x,NA,keep2$UniqueID.x)
}

pcr_data2022 <- pcr_data2022 %>% arrange(as.Date(pcr_data2022$Date))
pcr_data2022 <- pcr_data2022[!is.na(pcr_data2022$UniqueID),]

tempInfected2022 <- pcr_data2022[,c("UniqueID","Date","Infected")] %>%
  distinct() %>%
  spread(Date, Infected, fill = NA) %>% 
  group_by(UniqueID)

tempPCR2022 <- pcr_data2022[,c("UniqueID","Date","qPCR")] %>%
  distinct() %>%
  spread(Date, qPCR, fill = NA) %>% 
  group_by(UniqueID)

# Get same rows as in previous data components
id2022_table <- (left_join(EncounterHist2022[,1],
                           datayear2022[,1:6],
                           by="UniqueID") %>%
                   group_by(UniqueID) %>%
                   filter(row_number() == 1))[,1]

Infected2022 <- as.matrix(left_join(id2022_table,
                                    tempInfected2022,
                                    by="UniqueID")[,-1])[,unique_dates_2022]
PCR2022 <- as.matrix(left_join(id2022_table,
                               tempPCR2022,
                               by="UniqueID")[,-1])[,unique_dates_2022]


# Final versions (with all dates)

# 2020
tmp_InfectedMice2020 <- cbind(Infected2020,
                              matrix(NA,nrow(Infected2020),
                                     length(setdiff(full_dates_2020,colnames(Infected2020)))))
colnames(tmp_InfectedMice2020) <- c(time_captures_2020,setdiff(full_dates_2020,time_captures_2020))
Final_InfectedMice2020 <- get_BinomialObject_perTime(tmp_InfectedMice2020[,full_dates_2020])

# 2021
tmp_InfectedMice2021 <- cbind(Infected2021,
                              matrix(NA,nrow(Infected2021),
                                     length(setdiff(full_dates_2021,time_captures_2021))))
colnames(tmp_InfectedMice2021) <- c(time_captures_2021,setdiff(full_dates_2021,time_captures_2021))
Final_InfectedMice2021 <- get_BinomialObject_perTime(tmp_InfectedMice2021[,full_dates_2021])

# 2022
tmp_InfectedMice2022 <- cbind(Infected2022,
                              matrix(NA,nrow(Infected2022),
                                     length(setdiff(full_dates_2022,time_captures_2022))))
colnames(tmp_InfectedMice2022) <- c(time_captures_2022,setdiff(full_dates_2022,time_captures_2022))
Final_InfectedMice2022 <- get_BinomialObject_perTime(tmp_InfectedMice2022[,full_dates_2022])



########################################################################
# Antibody (ELISA) data
########################################################################

# 2020 
# adding protective_OspA indicator
elisa_data2020$protective_OspA <- rep(NA,nrow(elisa_data2020))
for(i in 1:nrow(elisa_data2020)){
  elisa_data2020$protective_OspA[i] <- ifelse(elisa_data2020$OspA.OD[i] > 0.80, 1, 0)
}

# Modify trap label
for(i in 1:nrow(elisa_data2020)){
  elisa_data2020$Trap[i] <- ifelse(str_length(elisa_data2020$Trap[i])==3,
                                   elisa_data2020$Trap[i],
                                 ch_insert(elisa_data2020$Trap[i],1,"0"))
}

# Trying to get some IDs back from RTV field data
id_for_merge_2020 <- datayear2020[,c("Date","Site","Trap","UniqueID","IDfull")]
tmp3 <- left_join(elisa_data2020, id_for_merge_2020, by = c("Date","Site"))
for(i in 1:nrow(tmp3)){
  if(!(is.na(tmp3$Trap.x[i])) & 
     !(is.na(tmp3$Trap.y[i])) & 
     (tmp3$Trap.y[i]==tmp3$Trap.x[i])){
    tmp3$IDfull.x[i] <- tmp3$IDfull.y[i]
  }
}

for(i in 1:nrow(elisa_data2020)){
  if(!is.na(elisa_data2020$IDfull[i])) next
  someIDs3 <- tmp3[(tmp3$UT.ear.number==i),] 
  keep3 <- someIDs3[someIDs3$IDfull.x==someIDs3$IDfull.y,]
  elisa_data2020$IDfull[i] <- ifelse(NA %in% keep3$IDfull.x,NA,keep3$IDfull.x)
}

elisa_data2020 <- elisa_data2020 %>% arrange(as.Date(elisa_data2020$Date))
elisa_data2020 <- elisa_data2020[!is.na(elisa_data2020$IDfull),]
elisa_data2020$UniqueID <- 2000 + elisa_data2020$IDfull

tempProtectiveOspA2020 <- elisa_data2020[,c("UniqueID","Date","protective_OspA")] %>%
  distinct() %>%
  spread(Date, protective_OspA, fill = NA) %>% 
  group_by(UniqueID)

templogOspA2020 <- elisa_data2020[,c("UniqueID","Date","log_OspA_OD")] %>%
  distinct() %>%
  spread(Date, log_OspA_OD, fill = NA) %>% 
  group_by(UniqueID)

# Get same rows as in previous data components
id2020_table <- (left_join(EncounterHist2020[,1],
                           datayear2020[,1:5],
                           by="UniqueID") %>%
                   group_by(UniqueID) %>%
                   filter(row_number() == 1))[,1]

tempProtectiveOspA2020_2 <- as.matrix(left_join(id2020_table,
                                                tempProtectiveOspA2020,
                                          by="UniqueID")) 
                                      
add_dates2020 <- setdiff(unique_dates_2020,colnames(tempProtectiveOspA2020_2)[-1])
tempProtectiveOspA2020_3 <- as.matrix(left_join(id2020_table,
                                                tempProtectiveOspA2020_2,
                                               by="UniqueID", copy = TRUE)[,-1]) 

tempProtectiveOspA2020_4 <- cbind(tempProtectiveOspA2020_3,
                                 matrix(NA,nrow(tempProtectiveOspA2020_3),length(add_dates2020)))            
colnames(tempProtectiveOspA2020_4)[ncol(tempProtectiveOspA2020_4):(ncol(tempProtectiveOspA2020_4)-length(add_dates2020)+1)] <- add_dates2020                                    
ProtectiveOspA2020 <- tempProtectiveOspA2020_4[,unique_dates_2020]
                                      


# 2021
# adding protective_OspA indicator
elisa_data2021$protective_OspA <- rep(NA,nrow(elisa_data2021))
for(i in 1:nrow(elisa_data2021)){
  elisa_data2021$protective_OspA[i] <- ifelse(elisa_data2021$OspA_OD[i] > 0.80, 1, 0)
}

# Modify trap label
for(i in 1:nrow(elisa_data2021)){
  elisa_data2021$Trap[i] <- ifelse(str_length(elisa_data2021$Trap[i])==3,
                                   elisa_data2021$Trap[i],
                                   ch_insert(elisa_data2021$Trap[i],1,"0"))
}

elisa_data2021 <- elisa_data2021 %>% arrange(as.Date(elisa_data2021$Date))
elisa_data2021 <- elisa_data2021[!is.na(elisa_data2021$UniqueID),]


tempProtectiveOspA2021 <- elisa_data2021[,c("UniqueID","Date","protective_OspA")] %>%
  distinct() %>%
  spread(Date, protective_OspA, fill = NA) %>% 
  group_by(UniqueID)

templogOspA2021 <- elisa_data2021[,c("UniqueID","Date","log_OspA_OD")] %>%
  distinct() %>%
  spread(Date, log_OspA_OD, fill = NA) %>% 
  group_by(UniqueID)

# Get same rows as in previous data components
id2021_table <- (left_join(EncounterHist2021[,1],
                           datayear2021[,1:6],
                           by="UniqueID") %>%
                   group_by(UniqueID) %>%
                   filter(row_number() == 1))[,1]

tempProtectiveOspA2021_2 <- as.matrix(left_join(id2021_table,
                                                tempProtectiveOspA2021,
                                                by="UniqueID")) 

add_dates2021 <- setdiff(unique_dates_2021,colnames(tempProtectiveOspA2021_2)[-1])
tempProtectiveOspA2021_3 <- as.matrix(left_join(id2021_table,
                                                tempProtectiveOspA2021_2,
                                                by="UniqueID", copy = TRUE)[,-1]) 

tempProtectiveOspA2021_4 <- cbind(tempProtectiveOspA2021_3,
                                  matrix(NA,nrow(tempProtectiveOspA2021_3),length(add_dates2021)))            
#colnames(tempProtectiveOspA2021_4)[ncol(tempProtectiveOspA2021_4):(ncol(tempProtectiveOspA2021_4)-length(add_dates2021)+1)] <- add_dates2021                                    
ProtectiveOspA2021 <- tempProtectiveOspA2021_4[,unique_dates_2021]



# 2022
# adding protective_OspA indicator
elisa_data2022$protective_OspA <- rep(NA,nrow(elisa_data2022))
for(i in 1:nrow(elisa_data2022)){
  elisa_data2022$protective_OspA[i] <- ifelse(elisa_data2022$OspA_OD[i] > 0.80, 1, 0)
}

# Modify trap label
for(i in 1:nrow(elisa_data2022)){
  elisa_data2022$Trap[i] <- ifelse(str_length(elisa_data2022$Trap[i])==3,
                                   elisa_data2022$Trap[i],
                                   ch_insert(elisa_data2022$Trap[i],1,"0"))
}
 
elisa_data2022 <- elisa_data2022 %>% arrange(as.Date(elisa_data2022$Date))
elisa_data2022 <- elisa_data2022[!is.na(elisa_data2022$UniqueID),]


tempProtectiveOspA2022 <- elisa_data2022[,c("UniqueID","Date","protective_OspA")] %>%
  distinct() %>%
  spread(Date, protective_OspA, fill = NA) %>% 
  group_by(UniqueID)

templogOspA2022 <- elisa_data2022[,c("UniqueID","Date","log_OspA_OD")] %>%
  distinct() %>%
  spread(Date, log_OspA_OD, fill = NA) %>% 
  group_by(UniqueID)

# Get same rows as in previous data components
id2022_table <- (left_join(EncounterHist2022[,1],
                           datayear2022[,1:6],
                           by="UniqueID") %>%
                   group_by(UniqueID) %>%
                   filter(row_number() == 1))[,1]

tempProtectiveOspA2022_2 <- as.matrix(left_join(id2022_table,
                                                tempProtectiveOspA2022,
                                                by="UniqueID")) 

add_dates2022 <- setdiff(unique_dates_2022,colnames(tempProtectiveOspA2022_2)[-1])
tempProtectiveOspA2022_3 <- as.matrix(left_join(id2022_table,
                                                tempProtectiveOspA2022_2,
                                                by="UniqueID", copy = TRUE)[,-1]) 

tempProtectiveOspA2022_4 <- cbind(tempProtectiveOspA2022_3,
                                  matrix(NA,nrow(tempProtectiveOspA2022_3),length(add_dates2022)))            
#colnames(tempProtectiveOspA2022_4)[ncol(tempProtectiveOspA2022_4):(ncol(tempProtectiveOspA2022_4)-length(add_dates2022)+1)] <- add_dates2022                                    
ProtectiveOspA2022 <- tempProtectiveOspA2022_4[,unique_dates_2022]


# Final versions (with all dates)

# 2020
tmp_ProtectiveOspA2020 <- cbind(ProtectiveOspA2020,
                                matrix(NA,nrow(ProtectiveOspA2020),
                                       length(setdiff(full_dates_2020,time_captures_2020))))
colnames(tmp_ProtectiveOspA2020) <- c(time_captures_2020,setdiff(full_dates_2020,time_captures_2020))
Final_ProtectiveOspA2020 <- get_BinomialObject_perTime(tmp_ProtectiveOspA2020[,full_dates_2020])


# 2021
tmp_ProtectiveOspA2021 <- cbind(ProtectiveOspA2021,
                                matrix(NA,nrow(ProtectiveOspA2021),
                                       length(setdiff(full_dates_2021,time_captures_2021))))
colnames(tmp_ProtectiveOspA2021) <- c(time_captures_2021,setdiff(full_dates_2021,time_captures_2021))
Final_ProtectiveOspA2021 <- get_BinomialObject_perTime(tmp_ProtectiveOspA2021[,full_dates_2021])


# 2022
tmp_ProtectiveOspA2022 <- cbind(ProtectiveOspA2022,
                                matrix(NA,nrow(ProtectiveOspA2022),
                                       length(setdiff(full_dates_2022,time_captures_2022))))
colnames(tmp_ProtectiveOspA2022) <- c(time_captures_2022,setdiff(full_dates_2022,time_captures_2022))
Final_ProtectiveOspA2022 <- get_BinomialObject_perTime(tmp_ProtectiveOspA2022[,full_dates_2022])


########################################################################
# Overall Counts
########################################################################

# MICE - Incidence
ll2020 <- length(unique(datayear2020$UniqueID))
Final_Incidence_Mice_2020 <- rep(NA, nsites2020)
SampleSize2020 <- c()
for(s in 1:nsites2020){
  v1 <- Site2020[1:ll2020]==s
  v2 <- c()
  for(i in 1:ll2020){
    v2[i] <- Final_InfectedMice2020[i,lastcapture2020[i]]
  }
  Final_Incidence_Mice_2020[s] <- sum(v1*v2,na.rm = TRUE)
  SampleSize2020[s] <- sum(v1)
}

ll2021 <- length(unique(datayear2021$UniqueID))
Final_Incidence_Mice_2021 <- rep(NA, nsites2021)
SampleSize2021 <- c()
for(s in 1:nsites2021){
  v1 <- Site2021[1:ll2021]==s
  v2 <- c()
  for(i in 1:ll2021){
    v2[i] <- Final_InfectedMice2021[i,lastcapture2021[i]]
  }
  Final_Incidence_Mice_2021[s] <- sum(v1*v2,na.rm = TRUE)
  SampleSize2021[s] <- sum(v1)
}

ll2022 <- length(unique(datayear2022$UniqueID))
Final_Incidence_Mice_2022 <- rep(NA, nsites2022)
SampleSize2022 <- c()
for(s in 1:nsites2022){
  v1 <- Site2022[1:ll2022]==s
  v2 <- c()
  for(i in 1:ll2022){
    v2[i] <- Final_InfectedMice2022[i,lastcapture2022[i]]
  }
  Final_Incidence_Mice_2022[s] <- sum(v1*v2,na.rm = TRUE)
  SampleSize2022[s] <- sum(v1)
}

# MICE - Protective OspA
Final_Protection_Mice_2020 <- rep(NA, nsites2020)
for(s in 1:nsites2020){
  v1 <- Site2020[1:ll2020]==s
  v2 <- c()
  for(i in 1:ll2020){
    v2[i] <- Final_ProtectiveOspA2020[i,lastcapture2020[i]]
  }
  Final_Protection_Mice_2020[s] <- sum(v1*v2,na.rm = TRUE)
}

Final_Protection_Mice_2021 <- rep(NA, nsites2021)
for(s in 1:nsites2021){
  v1 <- Site2021[1:ll2021]==s
  v2 <- c()
  for(i in 1:ll2021){
    v2[i] <- Final_ProtectiveOspA2021[i,lastcapture2021[i]]
  }
  Final_Protection_Mice_2021[s] <- sum(v1*v2,na.rm = TRUE)
}

Final_Protection_Mice_2022 <- rep(NA, nsites2022)
for(s in 1:nsites2022){
  v1 <- Site2022[1:ll2022]==s
  v2 <- c()
  for(i in 1:ll2022){
    v2[i] <- Final_ProtectiveOspA2022[i,lastcapture2022[i]]
  }
  Final_Protection_Mice_2022[s] <- sum(v1*v2,na.rm = TRUE)
}


# MICE - Nymphal Ticks
Final_Nymphals_Mice_2020 <- rep(NA, nsites2020)
SampleNymphals2020 <- c()
for(s in 1:nsites2020){
  v1 <- ticks_data2020_mouse[which(ticks_data2020_mouse$Site == sorted_sites[s]),]
  Final_Nymphals_Mice_2020[s] <- sum(v1$Infected_Ticks,na.rm = TRUE)
  SampleNymphals2020[s] <- length(which(ticks_data2020_mouse$Site == sorted_sites[s]))
}

Final_Nymphals_Mice_2021 <- rep(NA, nsites2021)
SampleNymphals2021 <- c()
for(s in 1:nsites2021){
  v1 <- ticks_data2021_mouse[which(ticks_data2021_mouse$Site == sorted_sites[s]),]
  Final_Nymphals_Mice_2021[s] <- sum(v1$Infected_Ticks,na.rm = TRUE)
  SampleNymphals2021[s] <- length(which(ticks_data2021_mouse$Site == sorted_sites[s]))
}

Final_Nymphals_Mice_2022 <- rep(NA, nsites2022)
SampleNymphals2022 <- c()
for(s in 1:nsites2022){
  v1 <- ticks_data2022_mouse[which(ticks_data2022_mouse$Site == sorted_sites2022[s]),]
  Final_Nymphals_Mice_2022[s] <- sum(v1$Infected_Ticks,na.rm = TRUE)
  SampleNymphals2022[s] <- length(which(ticks_data2022_mouse$Site == sorted_sites2022[s]))
}


# MICE - Dragging Ticks
Final_Nymphals_Drag_2020 <- rep(NA, nsites2020)
SampleDrag2020 <- c()
for(s in 1:nsites2020){
  v1 <- ticks_data2020_drag[which(ticks_data2020_drag$Site == sorted_sites[s]),]
  Final_Nymphals_Drag_2020[s] <- sum(v1$Infected_Ticks,na.rm = TRUE)
  SampleDrag2020[s] <- length(which(ticks_data2020_drag$Site == sorted_sites[s]))
}

Final_Nymphals_Drag_2021 <- rep(NA, nsites2021)
SampleDrag2021 <- c()
for(s in 1:nsites2021){
  v1 <- ticks_data2021_drag[which(ticks_data2021_drag$Site == sorted_sites[s]),]
  Final_Nymphals_Drag_2021[s] <- sum(v1$Infected_Ticks,na.rm = TRUE)
  SampleDrag2021[s] <- length(which(ticks_data2021_drag$Site == sorted_sites[s]))
}

Final_Nymphals_Drag_2022 <- rep(NA, nsites2022)
SampleDrag2022 <- c()
for(s in 1:nsites2022){
  v1 <- ticks_data2022_drag[which(ticks_data2022_drag$Site == sorted_sites2022[s]),]
  Final_Nymphals_Drag_2022[s] <- sum(v1$Infected_Ticks,na.rm = TRUE)
  SampleDrag2022[s] <- length(which(ticks_data2022_drag$Site == sorted_sites2022[s]))
}

########################################################################
# Ticks Count
########################################################################

# Ticks (counts-integer)
# 2020
datayear2020$Count <- 1:nrow(datayear2020)
tempTicks2020 <- (datayear2020[,c("Date","UniqueID","Count","Ticks")] %>%
                    group_by(UniqueID) %>%
                    pivot_wider(names_from = Date, values_from = Ticks))[,-2]

Ticks2020 <- (tempTicks2020 %>% 
                group_by(UniqueID) %>% 
                summarise(across(everything(), 
                                 function(x) ifelse(all(is.na(x)),
                                                    NA,
                                                    sum(x,na.rm = TRUE)))) %>%
                as.matrix)[,-1]

# Look at R script "Discrepancies_Number_of_Ticks_2020_2021.R" for explanation
# Run code on R script before running following line
Ticks2020[161,19] <- NA

tmp_Ticks2020 <- cbind(Ticks2020,
                       matrix(NA,nrow(Ticks2020),
                              length(setdiff(full_dates_2020,time_captures_2020))))
colnames(tmp_Ticks2020) <- c(time_captures_2020,setdiff(full_dates_2020,time_captures_2020))
Final_Ticks2020 <- tmp_Ticks2020[,full_dates_2020]


# 2021
datayear2021$Count <- 1:nrow(datayear2021)
tempTicks2021 <- (datayear2021[,c("Date","UniqueID","Count","Ticks")] %>%
                    group_by(UniqueID) %>%
                    pivot_wider(names_from = Date, values_from = Ticks))[,-2]

Ticks2021 <- (tempTicks2021 %>% 
                group_by(UniqueID) %>% 
                summarise(across(everything(), 
                                 function(x) ifelse(all(is.na(x)),
                                                    NA,
                                                    sum(x,na.rm = TRUE)))) %>%
                as.matrix)[,-1]

# Look at R script "Discrepancies_Number_of_Ticks_2020_2021.R" for explanation
# Run code on R script before running following lines
Ticks2021[15,6] <- NA
Ticks2021[64,3] <- NA
Ticks2021[286,6] <- NA
Ticks2021[290,5] <- NA
Ticks2021[303,14] <- NA

tmp_Ticks2021 <- cbind(Ticks2021,
                       matrix(NA,nrow(Ticks2021),
                              length(setdiff(full_dates_2021,time_captures_2021))))
colnames(tmp_Ticks2021) <- c(time_captures_2021,setdiff(full_dates_2021,time_captures_2021))
Final_Ticks2021 <- tmp_Ticks2021[,full_dates_2021]


# 2022
datayear2022$Count <- 1:nrow(datayear2022)
tempTicks2022 <- (datayear2022[,c("Date","UniqueID","Count","Ticks")] %>%
                    group_by(UniqueID) %>%
                    pivot_wider(names_from = Date, values_from = Ticks))[,-2]
tempTicks2022 <- as.data.frame(tempTicks2022)
tempTicks2022[tempTicks2022==999]<-NA

Ticks2022 <- (tempTicks2022 %>% 
                group_by(UniqueID) %>% 
                summarise(across(everything(), 
                                 function(x) ifelse(all(is.na(x)),
                                                    NA,
                                                    sum(x,na.rm = TRUE)))) %>%
                as.matrix)[,-1]


tmp_Ticks2022 <- cbind(Ticks2022,
                       matrix(NA,nrow(Ticks2022),
                              length(setdiff(full_dates_2022,time_captures_2022))))
colnames(tmp_Ticks2022) <- c(time_captures_2022,setdiff(full_dates_2022,time_captures_2022))
Final_Ticks2022 <- tmp_Ticks2022[,full_dates_2022]

########################################################################
# Producing Timeline Plots (Mouse ticks vs Dragged ticks)
########################################################################
   
need_plots <- FALSE
if(need_plots == TRUE){
  time_ticks_2020 <- data.frame(date = c(colnames(InfectedTicks_Mouse2020),
                                         as.vector(na.omit(unique(ticks_data2020_drag$Date)))),
                                type = c(rep("Mouse_Skin",unique_number_dates_2020),
                                         rep("Dragged_Field",length(as.vector(na.omit(unique(ticks_data2020_drag$Date)))))))
  
  g1 <- ggplot(data = time_ticks_2020, aes(x = date, y = type, shape = factor(type))) +
    geom_point(aes(colour = factor(type)), size = 3) +
    labs(x = "Date",
         y = "",
         subtitle = "Year 2020",
         colour = "Tick Collection Type") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    theme(legend.position="none")
  
  
  time_ticks_2021 <- data.frame(date = c(colnames(InfectedTicks_Mouse2021),
                                         unique(format(ticks_data2021_drag$Date, "%Y-%m-%d"))),
                                type = c(rep("Mouse_Skin",unique_number_dates_2021),
                                         rep("Dragged_Field",length(unique(format(ticks_data2021_drag$Date, "%Y-%m-%d"))))))
  
  g2 <- ggplot(data = time_ticks_2021, aes(x = date, y = type, shape = factor(type))) +
    geom_point(aes(colour = factor(type)), size = 3) +
    labs(x = "Date",
         y = "",
         subtitle = "Year 2021",
         colour = "Tick Collection Type") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    theme(legend.position="none")
  
  
  
  time_ticks_2022 <- data.frame(date = c(colnames(InfectedTicks_Mouse2022),
                                         unique(format(ticks_data2022_drag$Date, "%Y-%m-%d"))),
                                type = c(rep("Mouse_Skin",unique_number_dates_2022),
                                         rep("Dragged_Field",length(unique(format(ticks_data2022_drag$Date, "%Y-%m-%d"))))))
  
  g3 <- ggplot(data = time_ticks_2022, aes(x = date, y = type, shape = factor(type))) +
    geom_point(aes(colour = factor(type)), size = 3) +
    labs(x = "Date",
         y = "",
         subtitle = "Year 2022",
         colour = "Tick Collection Type") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    theme(legend.position="none")
  
  
  plot1 <-  annotate_figure(ggarrange(g1,g2,g3, labels = c("A","B", "C"), 
                                      common.legend = TRUE, legend = "none",
                                      ncol = 1), 
                            top = text_grob("RTV Field Study - Sampling/Dragging Dates by Tick Collection Type\n", 
                                            color = "black", face = "bold", size = 14))
  plot1
  
}


