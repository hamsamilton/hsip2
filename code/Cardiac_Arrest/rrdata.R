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
rr <- read.csv("/Users/User/Box Sync/Projects/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/Chartevents_618_rr.csv")

# Patient IDs
ca_ids <- read.csv("/Users/User/Box Sync/Projects/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/ca_ids.csv")

# Sample Size
nrow(unique(data.frame(rr$subject_id)))
nrow(unique(data.frame(rr$hadm_id))) # checking that there are unique admissions per patient - looks good!

# Select patients with at least 200 observations 
rr_counts <- rr %>% group_by(subject_id) %>% dplyr::count(num_obs = n())
rr_counts_200 <- rr_counts %>% filter(num_obs > 200) #(specify num_obs > 200 to avoid error with the entropy function!)
rr <- rr %>% filter(subject_id %in% rr_counts_200$subject_id)

## Data Pre-processing

#Format time
rr$charttime_conv <- strptime(as.character(rr$charttime), format = "%Y-%m-%d %H:%M:%S", tz = "EST")
rr$charttime_conv <- as.POSIXct(rr$charttime_conv, tz="EST")

# Transform Value
rr$value <- as.numeric(as.character(rr$value))

# Combine physiologic data with subject ids
rr = rr %>% inner_join(ca_ids, by = c("subject_id","hadm_id"))

# remove null values
rr <- rr %>% filter(!value == "NULL")

## Create time plots
# x = subject_id
plot_rr <- function(x) { 
  
  # arrange by Date
  df <- rr %>% filter(subject_id == x) %>% arrange(charttime_conv)
  
  #Title: create a vector that will contain the ID number and can be used to title the    graph by patient ID 
  title <- x
  
  #Create limits to identify the start and end dates of the time series
  first <- first(df$charttime_conv)
  last <- last(df$charttime_conv)
  
  rr <- df %>% filter(!is.na(value)) %>% select(subject_id, charttime_conv, value) %>%
    ggplot(aes(x = charttime_conv, y = value)) + 
    geom_point(color = "slateblue") + 
    geom_line(color = "slateblue") +
    #geom_hline(yintercept=3) + #could demark dangerous (low or high) rrs
    scale_y_continuous(name = "Resp Rate") + #define the y-axis limits
    #scale_x_continuous(name = "Time", date_break = "1 hour", date_labels = "%H", limits = c(first, last)) + #, limits = c(-6, 0), breaks = c(-6:0)) +
    scale_x_datetime(name = "Time", date_break = "6 hour", date_labels = "%H", limits = c(first, last)) +
    ggtitle(paste0('ID:', title))  #add a title to identify patients 
  
  #remove legend and x-axis
  rr <- rr  + theme(legend.position = "none")
  
  #remove gray background and add a border around the plot
  rr <- rr + theme(panel.border = element_rect(fill=NA,color="black", size=0.5,
                                               linetype="solid"))
  
  
  rr <- rr + theme_bw()
  
  return(rr)
}

# arrange by time
rr <- rr %>% group_by(subject_id) %>% arrange(charttime_conv)

m <- CalcEnSlide(df = rr, 
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

plot_rr(20866) #high ent
plot_rr(6469) # low entropy

# Select the firest 200 measurements for each patients
rr_200 <- rr %>% group_by(subject_id) %>% 
  arrange(charttime_conv) %>% 
  slice(1:200) 

# check counts 
View(rr_200 %>% group_by(subject_id) %>% dplyr::summarize(count = n()))
nrow(unique(data.frame(rr_200$subject_id)))

# Determine the max, min, mean, sd, and variance of the first 200 measurements
library(dplyr) #reload, interferring with plyr
rr_200 <- rr_200 %>% group_by(subject_id) %>% 
  dplyr::mutate(min_rr = min(value, na.rm = TRUE)) %>% 
  dplyr::mutate(max_rr = max(value, na.rm = TRUE)) %>% 
  dplyr::mutate(mean_rr = mean(value, na.rm = TRUE)) %>%
  dplyr::mutate(sd_rr = sd(value, na.rm = TRUE)) %>% 
  dplyr::mutate(var_rr = var(value, na.rm = TRUE)) %>% 
  select(subject_id, min_rr, max_rr, mean_rr, sd_rr, var_rr) %>%
  slice(1L)

# format entropy data
m_ent <- m %>% select(subject_id, entropy)
colnames(m_ent)[2] <- "entropy_rr"

# add entropy + demographics
rr_200 <- merge(rr_200, m_ent, by = "subject_id")

# save as csv
write.csv(rr_200, "/Users/User/Box Sync/Projects/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/rr_200.csv")
