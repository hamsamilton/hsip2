library(TSEntropies)
library(dplyr)
library(lubridate)
library(ggplot2)
library(plyr)

#Calculates entropy using sliding window approach
CalcEnSlide <- function(entmeasure = "FastSampEn",df,grouplist,r_thrsh = .2,size = 200, step = 1) {
  
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
spo2 <- read.csv("C:/Users/User/Box Sync/Projects/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/Chartevents_646_spo2.csv")

# Patient IDs
ca_ids <- read.csv("/Users/User/Box Sync/Projects/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/ca_ids.csv")

# Sample Size
nrow(unique(data.frame(spo2$subject_id)))
nrow(unique(data.frame(spo2$hadm_id))) # checking that there are unique admissions per patient - looks good!

# Select patients with at least 200 observations 
spo2_counts <- spo2 %>% group_by(subject_id) %>% dplyr::count(num_obs = n())
spo2_counts_200 <- spo2_counts %>% filter(num_obs > 200) #(specify num_obs > 200 to avoid error with the entropy function!)
spo2 <- spo2 %>% filter(subject_id %in% spo2_counts_200$subject_id)

## Data Pre-processing

#Format time
spo2$charttime_conv <- strptime(as.character(spo2$charttime), format = "%Y-%m-%d %H:%M:%S", tz = "EST")
spo2$charttime_conv <- as.POSIXct(spo2$charttime_conv, tz="EST")

# Transform Value
spo2$value <- as.numeric(as.character(spo2$value))

# Combine physiologic data with subject ids
spo2 = spo2 %>% inner_join(ca_ids, by = c("subject_id","hadm_id"))

# remove null values
spo2 <- spo2 %>% filter(!value == "NULL")

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
    #geom_hline(yintercept=3) + #could demark dangerous (low or high) spo2
    scale_y_continuous(name = "SPo2") + #define the y-axis limits
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

m <- CalcEnSlide(df = spo2, 
                 grouplist = c("hospital_expire_flag","subject_id","hadm_id"), 
                 size = 200,
                 step = 1,
                 r_thrsh = .2)

lm1 <- glm(hospital_expire_flag ~ entropy, family = "binomial", data = m)
summary(lm1)

m$hospital_expire_flag <- as.factor(m$hospital_expire_flag)
ggplot(m,aes(x = hospital_expire_flag, y = entropy, fill = hospital_expire_flag)) +
  geom_boxplot() +
  geom_point()

plot_spo2(4082) #high ent - mortality
plot_spo2(23120) # low entropy - mortality

# Select the firest 200 measurements for each patients
spo2_200 <- spo2 %>% group_by(subject_id) %>% 
  arrange(charttime_conv) %>% 
  slice(1:200) 

# check counts 
View(spo2_200 %>% group_by(subject_id) %>% dplyr::summarize(count = n()))
nrow(unique(data.frame(spo2_200$subject_id)))

# Determine the max, min, mean, sd, and variance of the first 200 measurements
library(dplyr) #reload, interferring with plyr
spo2_200 <- spo2_200 %>% group_by(subject_id) %>% 
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
spo2_200 <- merge(spo2_200, m_ent, by = "subject_id")

# save as csv
write.csv(spo2_200, "/Users/User/Box Sync/Projects/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/spo2_200.csv", row.names = FALSE)
