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

# Heart Rate Data
temp <- read.csv("/Users/samuelhamilton/Downloads/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/Chartevents_678_temp_F.csv")

# Patient IDs
ca_ids <- read.csv("/Users/samuelhamilton/Downloads/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/ca_ids.csv")

nrow(unique(data.frame(temp$subject_id)))
nrow(unique(data.frame(temp$hadm_id))) # checking that there are unique admissions per patient - looks good!


temp_counts <- temp %>% group_by(hadm_id) %>% dplyr::count(num_obs = n())

temp_counts_200 <- temp_counts %>% filter(num_obs >= 200)
temp <- temp %>% 
  filter(hadm_id %in% temp_counts_200$hadm_id)
#We will need to check the time interval/perhaps standardize by time in some fashion (ie: if multiple measurements within the same time frame, do we want to take the mean, max, min? min or max he
## Data Pre-processing

# CHF Patients
temp$charttime_conv <- strptime(as.character(temp$charttime), format = "%Y-%m-%d %H:%M:%S", tz = "EST")
temp$charttime_conv <- as.POSIXct(temp$charttime_conv, tz="EST")

temp$value <- as.numeric(as.character(temp$value))

temp <- temp %>% inner_join(ca_ids, by = c("subject_id","hadm_id"))

temp <- temp %>% filter(!value == "NULL")



#Take the min, max, mean, and sd when measurements occur within the same time interval (if we want this step!)
#temp <- temp %>% group_by(subject_id, charttime_conv) %>% mutate(max_temp = max(value, na.rm = TRUE)) %>% mutate(min_temp = min(value, na.rm = TRUE)) %>% mutate(mean_temp = mean(value, na.rm = TRUE)) %>% mutate(sd_temp = sd(value, na.rm = TRUE)) 

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
    scale_y_continuous(name = "Degrees F") + #define the y-axis limits
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

m <- CalcEnSlide(df = temp, 
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
plot_temp(27850)
plot_temp(849) # high entropy
plot_temp(7223) # low entropy
plot_temp(849)# low entropy
