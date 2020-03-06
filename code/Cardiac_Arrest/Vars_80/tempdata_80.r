library(TSEntropies)
library(dplyr)
library(lubridate)
library(ggplot2)
#library(plyr)

#Calculates entropy using sliding window approach
CalcEnSlide <- function(entmeasure = "FastSampEn",df,grouplist,r_thrsh = .2,size = 80, step = 1) {
  
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
temp <- read.csv("/Users/User/Box Sync/Projects/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/Chartevents_678_679_temp_F.csv")

## Cleaning
# Transform Value to numeric
temp$value <- as.numeric(as.character(temp$value))

# remove null values
temp <- temp %>% filter(!value == "NULL")

# Clean up temp data by removing values < 80
temp <- temp %>% filter(!value < 80)

# Patient IDs
ca_ids <- read.csv("/Users/User/Box Sync/Projects/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/ca_ids.csv")

# Sample Size
nrow(unique(data.frame(temp$subject_id)))
nrow(unique(data.frame(temp$hadm_id))) # checking that there are unique admissions per patient - looks good!

## Data Pre-processing
#Format time
temp$charttime_conv <- strptime(as.character(temp$charttime), format = "%Y-%m-%d %H:%M:%S", tz = "EST")
temp$charttime_conv <- as.POSIXct(temp$charttime_conv, tz="EST")

# Combine physiologic data with subject ids
temp = temp %>% inner_join(ca_ids, by = c("subject_id","hadm_id"))

# clean-up dataframe
temp <- temp %>% select(subject_id, hadm_id, hospital_expire_flag, charttime_conv, value)

# remove any duplicate measurements
temp <- temp %>% distinct

# Select only first 48 hours of admission
temp <- temp %>% group_by(subject_id) %>% arrange(subject_id, charttime_conv) %>% 
  mutate(first_date = first(charttime_conv)) %>% 
  mutate(time_diff = difftime(charttime_conv, first_date, unit = "hours")) %>% 
  filter(time_diff <= 72) %>% 
  mutate(count = n()) %>% 
  filter(count > 80)

## Create time plots
# x = subject_id
plot_temp <- function(x) { 
  
  # arrange by Date
  df <- temp %>% filter(subject_id == x) %>% arrange(charttime_conv)
  
  #Title: create a vector that will contain the ID number and can be used to title the    graph by patient ID 
  title <- x
  
  #Create limits to identify the start and end dates of the time series
  first <- first(df$charttime_conv)
  last <- last(df$charttime_conv)
  
  temp <- df %>% filter(!is.na(value)) %>% select(subject_id, charttime_conv, value) %>%
    ggplot(aes(x = charttime_conv, y = value)) + 
    geom_point(color = "slateblue") + 
    geom_line(color = "slateblue") +
    #geom_hline(yintercept=3) + #could demark dangerous (low or high) temps
    scale_y_continuous(name = "Temp (F)") + #define the y-axis limits
    #scale_x_continuous(name = "Time", date_break = "1 hour", date_labels = "%H", limits = c(first, last)) + #, limits = c(-6, 0), breaks = c(-6:0)) +
    scale_x_datetime(name = "Time", date_break = "6 hour", date_labels = "%H", limits = c(first, last)) +
    ggtitle(paste0('ID:', title))  #add a title to identify patients 
  
  #remove legend and x-axis
  temp <- temp  + theme(legend.position = "none")
  
  #remove gray background and add a border around the plot
  temp <- temp + theme(panel.border = element_rect(fill=NA,color="black", size=0.5,
                                               linetype="solid"))
  
  
  temp <- temp + theme_bw()
  
  return(temp)
}

# arrange by time
temp <- temp %>% group_by(subject_id) %>% arrange(charttime_conv)

library(plyr)
library(dplyr)
m <- CalcEnSlide(df = temp, 
                 grouplist = c("hospital_expire_flag","subject_id","hadm_id"), 
                 size = 80,
                 step = 1,
                 r_thrsh = .2)

lm1 <- glm(hospital_expire_flag ~ entropy, family = "binomial", data = m)
summary(lm1)

m$hospital_expire_flag <- as.factor(m$hospital_expire_flag)
ggplot(m,aes(x = hospital_expire_flag, y = entropy, fill = hospital_expire_flag)) +
  geom_boxplot() +
  geom_point()

# both deceased
plot_temp(10993) #low ent
plot_temp(27850) # high entropy

# deceased
plot_temp(11438) # high

# Select the first 80 measurements for each patients
temp_80 <- temp %>% group_by(subject_id) %>% 
  arrange(charttime_conv) %>% 
  slice(1:80) 

# check counts 
View(temp_80 %>% group_by(subject_id) %>% dplyr::summarize(count = n()))
nrow(unique(data.frame(temp_80$subject_id)))

# Determine the max, min, mean, sd, and variance of the first 80 measurements
library(dplyr) #reload, interferring with plyr
temp_80 <- temp_80 %>% group_by(subject_id) %>% 
  dplyr::mutate(min_temp = min(value, na.rm = TRUE)) %>% 
  dplyr::mutate(max_temp = max(value, na.rm = TRUE)) %>% 
  dplyr::mutate(mean_temp = mean(value, na.rm = TRUE)) %>%
  dplyr::mutate(sd_temp = sd(value, na.rm = TRUE)) %>% 
  dplyr::mutate(var_temp = var(value, na.rm = TRUE)) %>% 
  select(subject_id, min_temp, max_temp, mean_temp, sd_temp, var_temp) %>%
  slice(1L)

# format entropy data
m_ent <- m %>% select(subject_id, entropy)
colnames(m_ent)[2] <- "entropy_temp"

# add entropy + demographics
temp_80 <- merge(temp_80, m_ent, by = "subject_id")

# save as csv
write.csv(temp_80, "/Users/User/Box Sync/Projects/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/Temp_80.csv", row.names = FALSE)
