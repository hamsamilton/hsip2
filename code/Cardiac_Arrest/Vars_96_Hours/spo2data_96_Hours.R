library(TSEntropies)
library(dplyr)
library(lubridate)
library(ggplot2)
#library(plyr)

#Calculates entropy using sliding window approach
CalcEnSlide <- function(entmeasure = "FastSampEn",df,grouplist,r_thrsh = .2,size = 100, step = 1) {
  
  #calculate  sd for all ts 
  allsd <- sd(df$value,na.rm = T)
  print(allsd)
  #Calculates entropy for a TS across a set of sliding windows
  calc_En_slidingwindow <- function(entmeasure,varofint,size,step,r_thresh,allsd) {
    slidingentlist <- c() #init list 
    print(length(varofint))
    for(i in seq(from = 1, to = (length(varofint) - size), by = step)){
      
      window <- varofint[i:(i+size)]
      EN <- do.call(entmeasure,list(window, r = r_thresh*allsd))
      slidingentlist <- append(slidingentlist,EN)
    }
    return(slidingentlist)
  }
  
  #Uses above to calc mean sliding window entropy function for each var in varslist. 
  #
  k <- ddply(df,
             grouplist,
             here(summarise),
             entropy = median(calc_En_slidingwindow(entmeasure = entmeasure,
                                                    varofint = value,
                                                    size = size,
                                                    step = step,
                                                    r_thresh = r_thrsh,
                                                    allsd = allsd),
                              na.rm = T
             )
  )
  return(k)
}


# Data
spo2 <- read.csv("/Users/User/Box Sync/Projects/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/Chartevents_646_spo2.csv")

## Cleaning
# Transform Value to numeric
spo2$value <- as.numeric(as.character(spo2$value))

# remove null values
spo2 <- spo2 %>% filter(!value == "NULL")

# Patient IDs
ca_ids <- read.csv("/Users/User/Box Sync/Projects/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/ca_ids.csv")

# Sample Size
nrow(unique(data.frame(spo2$subject_id)))
nrow(unique(data.frame(spo2$hadm_id))) # checking that there are unique admissions per patient - looks good!

## Data Pre-processing
#Format time
spo2$charttime_conv <- strptime(as.character(spo2$charttime), format = "%Y-%m-%d %H:%M:%S", tz = "EST")
spo2$charttime_conv <- as.POSIXct(spo2$charttime_conv, tz="EST")

# Combine physiologic data with subject ids
spo2 = spo2 %>% inner_join(ca_ids, by = c("subject_id","hadm_id"))

# clean-up dataframe
spo2 <- spo2 %>% select(subject_id, hadm_id, hospital_expire_flag, charttime_conv, value)

# remove any duplicate measurements
spo2 <- spo2 %>% distinct

# Indentify first admission date
# Import Admissions
admissions <- read.csv("C:/Users/User/Box Sync/Projects/Mimic_HSIP/Mimic_Data/admissions.csv")
admissions$admittime_conv <- strptime(as.character(admissions$admittime), format = "%Y-%m-%d %H:%M:%S", tz = "EST")
admissions$admittime_conv <- as.POSIXct(admissions$admittime_conv, tz="EST")

# add admission time to the temperature dataset
admissions <- admissions %>% select(subject_id, hadm_id, admittime_conv) %>% distinct
spo2 <- merge(spo2, admissions, by = c("subject_id", "hadm_id"))

# Select only first 96 hours of admission
spo2 <- spo2 %>% group_by(subject_id) %>% arrange(subject_id, charttime_conv) %>% 
  mutate(time_diff = difftime(charttime_conv, admittime_conv, unit = "hours")) %>% 
  filter(time_diff <= 96) %>% 
  mutate(count = n()) 

# Determine median count - how many measurements do most patients have 
# within the first 96 hours
median(spo2$count)
mean(spo2$count)

spo2 <- spo2 %>% filter(count > 100) 

# Make sure these patients have about 4 days worth of measurements
spo2 <- spo2 %>% group_by(subject_id) %>% arrange(subject_id, charttime_conv) %>% 
  mutate(last_time_diff = last(time_diff)) %>% filter(last_time_diff >= 92)

nrow(unique(data.frame(spo2$subject_id)))

## Create time plots
# x = subject_id
plot_spo2 <- function(x) { 
  
  # arrange by Date
  df <- spo2 %>% filter(subject_id == x) %>% arrange(charttime_conv)
  
  #Title: create a vector that will contain the ID number and can be used to title the    graph by patient ID 
  title <- x
  
  #Create limits to identify the start and end dates of the time series
  first <- first(df$charttime_conv)
  last <- last(df$charttime_conv)
  
  spo2 <- df %>% filter(!is.na(value)) %>% select(subject_id, charttime_conv, value) %>%
    ggplot(aes(x = charttime_conv, y = value)) + 
    geom_point(color = "slateblue") + 
    geom_line(color = "slateblue") +
    #geom_hline(yintercept=3) + #could demark dangerous (low or high) spo2s
    scale_y_continuous(name = "spo2") + #define the y-axis limits
    #scale_x_continuous(name = "Time", date_break = "1 hour", date_labels = "%H", limits = c(first, last)) + #, limits = c(-6, 0), breaks = c(-6:0)) +
    scale_x_datetime(name = "Time", date_break = "6 hour", date_labels = "%H", limits = c(first, last)) +
    ggtitle(paste0('ID:', title))  #add a title to identify patients 
  
  #remove legend and x-axis
  spo2 <- spo2  + theme(legend.position = "none")
  
  #remove gray background and add a border around the plot
  spo2 <- spo2 + theme(panel.border = element_rect(fill=NA,color="black", size=0.5,
                                                   linetype="solid"))
  
  
  spo2 <- spo2 + theme_bw()
  
  return(spo2)
}

# arrange by time
spo2 <- spo2 %>% group_by(subject_id) %>% arrange(charttime_conv)

library(plyr)
library(dplyr)
m <- CalcEnSlide(df = spo2, 
                 grouplist = c("hospital_expire_flag","subject_id","hadm_id"), 
                 size = 100,
                 step = 1,
                 r_thrsh = .2)

lm1 <- glm(hospital_expire_flag ~ entropy, family = "binomial", data = m)
summary(lm1)

m$hospital_expire_flag <- as.factor(m$hospital_expire_flag)
ggplot(m,aes(x = hospital_expire_flag, y = entropy, fill = hospital_expire_flag)) +
  geom_boxplot() +
  geom_point()

#plot_spo2() 

# Select the first 80 measurements for each patients
spo2_96 <- spo2 %>% group_by(subject_id) %>% 
  arrange(charttime_conv) %>% 
  slice(1:100) 

# check counts 
View(spo2_96 %>% group_by(subject_id) %>% dplyr::summarize(count = n()))
nrow(unique(data.frame(spo2_96$subject_id)))

# Determine the max, min, mean, sd, and variance of the first 80 measurements
library(dplyr) #reload, interferring with plyr
spo2_96 <- spo2_96 %>% group_by(subject_id) %>% 
  dplyr::mutate(min_spo2 = min(value, na.rm = TRUE)) %>% 
  dplyr::mutate(max_spo2 = max(value, na.rm = TRUE)) %>% 
  dplyr::mutate(mean_spo2 = mean(value, na.rm = TRUE)) %>%
  dplyr::mutate(sd_spo2 = sd(value, na.rm = TRUE)) %>% 
  dplyr::mutate(var_spo2 = var(value, na.rm = TRUE)) %>% 
  select(subject_id, min_spo2, max_spo2, mean_spo2, sd_spo2, var_spo2) %>%
  slice(1L)

# format entropy data
m_ent <- m %>% select(subject_id, entropy)
colnames(m_ent)[2] <- "entropy_spo2"

# add entropy + demographics
spo2_96 <- merge(spo2_96, m_ent, by = "subject_id")

# save as csv
write.csv(spo2_96, "/Users/User/Box Sync/Projects/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/spo2_96.csv", row.names = FALSE)
