
library(dplyr)
library(lubridate)
library(ggplot2)

#Calculates entropy using sliding window approach
CalcEnSlide <- function(entmeasure = "FastSampEn",df,grouplist,r_thrsh = .2,size = 200, step = 1) {
  
  #calculate  sd for all ts 
  allsd <- sd(df$value,na.rm = T)
  #Calculates entropy for a TS across a set of sliding windows
  calc_En_slidingwindow <- function(entmeasure,size,ofint,step,r_thresh,allsd) {
    slidingentlist <- c() #init list 
    for(i in seq(from = 1, to = (length(ofint) - size), by = step)){
      
      window <- ofint[i:(i+size)]
      EN <- do.call(entmeasure,list(window, r = r_thresh*allsd))
      slidingentlist <- append(slidingentlist,EN)
    }
    return(slidingentlist)
  }
  
  #Uses above to calc mean sliding window entropy function for each var in varslist. 
  #
  k <- ddply(df,
             grouplist,
             summarise,
             entropy = median(calc_En_slidingwindow(entmeasure = entmeasure,
                                                    ofint = value,
                                                    size = size,
                                                    step = step,
                                                    r_thresh = r_thrsh,
                                                    allsd = allsd),
                              na.rm = T
             )
  )
  return(k)
}

# Heart Rate Data
hr <- read.csv("/Users/samuelhamilton/Downloads/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/Chartevents_211_hr.csv")

# Patient IDs
ca_ids <- read.csv("/Users/samuelhamilton/Downloads/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/ca_ids.csv")

hr <- merge(hr,ca_ids,by = "subject_id")


 nrow(unique(data.frame(hr$subject_id)))
nrow(unique(data.frame(hr$hadm_id))) # checking that there are unique admissions per patient - looks good!


hr_counts <- hr %>% group_by(subject_id) %>% count(num_obs = n())

hr_counts_200 <- hr_counts %>% filter(num_obs >= 220)

hr <- hr %>% 
  filter(subject_id %in% hr_counts_200$subject_id)
## Data Pre-processing


# CHF Patients
hr$charttime_conv <- strptime(as.character(hr$charttime), format = "%Y-%m-%d %H:%M:%S", tz = "EST")
hr$charttime_conv <- as.POSIXct(hr$charttime_conv, tz="EST")

hr$value <- as.numeric(as.character(hr$value))



hr <- hr %>% filter(!value == "NULL") %>% 
  filter(value > 0)



#Take the min, max, mean, and sd when measurements occur within the same time interval (if we want this step!)

#hr <- hr %>% 
# group_by(subject_id, charttime_conv) %>% 
#   mutate(max_hr = max(value, na.rm = TRUE)) %>%
#   mutate(min_hr = min(value, na.rm = TRUE)) %>%
#   mutate(mean_hr = mean(value, na.rm = TRUE)) %>% 
#   mutate(sd_hr = sd(value, na.rm = TRUE)) 


## Create time plots

# x = subject_id
plot_hr <- function(x) { 
  
  # arrange by Date
  df <- hr %>% filter(subject_id == x) %>% arrange(charttime_conv)
  
  #Title: create a vector that will contain the ID number and can be used to title the    graph by patient ID 
  title <- x
  
  #Create limits to identify the start and end dates of the time series
  first <- first(df$charttime_conv)
  last <- last(df$charttime_conv)
  
  hr <- df %>% filter(!is.na(value)) %>% select(subject_id, charttime_conv, value) %>%
    ggplot(aes(x = charttime_conv, y = value)) + 
    geom_point(color = "slateblue") + 
    geom_line(color = "slateblue") +
    #geom_hline(yintercept=3) + #could demark dangerous (low or high) HRs
    scale_y_continuous(name = "BMP") + #define the y-axis limits
    #scale_x_continuous(name = "Time", date_break = "1 hour", date_labels = "%H", limits = c(first, last)) + #, limits = c(-6, 0), breaks = c(-6:0)) +
    scale_x_datetime(name = "Time", date_break = "6 hour", date_labels = "%H", limits = c(first, last)) +
    ggtitle(paste0('ID:', title))  #add a title to identify patients 
  
  #remove legend and x-axis
  hr <- hr  + theme(legend.position = "none")
  
  #remove gray background and add a border around the plot
  hr <- hr + theme(panel.border = element_rect(fill=NA,color="black", size=0.5,
                                               linetype="solid"))
  
  
  hr <- hr + theme_bw()
  
  return(hr)
}



m <- CalcEnSlide(df = hr,grouplist = c("hospital_expire_flag","subject_id"), size = 200,step = 1,r_thrsh = .2)

plot_hr(31171) # low ent 
plot_hr(14244) #high ent

m$hospital_expire_flag <- as.factor(m$hospital_expire_flag)
ggplot(m,aes(x = hospital_expire_flag, y = entropy, fill = hospital_expire_flag)) +
  geom_boxplot()

lm1 <- glm(hospital_expire_flag ~ entropy, family = "binomial", data = m)
summary(lm1)
