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
bp <- read.csv("C:/Users/User/Box Sync/Projects/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/Chartevents_51_arterial_sbp.csv")

# Patient IDs
ca_ids <- read.csv("/Users/User/Box Sync/Projects/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/ca_ids.csv")

# Sample Size
nrow(unique(data.frame(bp$subject_id)))
nrow(unique(data.frame(bp$hadm_id))) # checking that there are unique admissions per patient - looks good!

# Select patients with at least 200 observations 
bp_counts <- bp %>% group_by(subject_id) %>% dplyr::count(num_obs = n())
bp_counts_200 <- bp_counts %>% filter(num_obs > 201) #(specify num_obs > 201 to avoid error with the entropy function!)
bp <- bp %>% filter(subject_id %in% bp_counts_200$subject_id)

## Data Pre-processing

#Format time
bp$charttime_conv <- strptime(as.character(bp$charttime), format = "%Y-%m-%d %H:%M:%S", tz = "EST")
bp$charttime_conv <- as.POSIXct(bp$charttime_conv, tz="EST")

# Transform Value
bp$value <- as.numeric(as.character(bp$value))

# Combine physiologic data with subject ids
bp = bp %>% inner_join(ca_ids, by = c("subject_id","hadm_id"))

# remove null values
bp <- bp %>% filter(!value == "NULL")

## Create time plots
# x = subject_id
plot_bp <- function(x) { 
  
  # arrange by Date
  df <- bp %>% filter(subject_id == x) %>% arrange(charttime_conv)
  
  #Title: create a vector that will contain the ID number and can be used to title the    graph by patient ID 
  title <- x
  
  #Create limits to identify the start and end dates of the time series
  first <- first(df$charttime_conv)
  last <- last(df$charttime_conv)
  
  bp <- df %>% filter(!is.na(value)) %>% select(subject_id, charttime_conv, value) %>%
    ggplot(aes(x = charttime_conv, y = value)) + 
    geom_point(color = "slateblue") + 
    geom_line(color = "slateblue") +
    #geom_hline(yintercept=3) + #could demark dangerous (low or high) bp
    scale_y_continuous(name = "HHmg") + #define the y-axis limits
    #scale_x_continuous(name = "Time", date_break = "1 hour", date_labels = "%H", limits = c(first, last)) + #, limits = c(-6, 0), breaks = c(-6:0)) +
    scale_x_datetime(name = "Time", date_break = "6 hour", date_labels = "%H", limits = c(first, last)) +
    ggtitle(paste0('ID:', title))  #add a title to identify patients 
  
  #remove legend and x-axis
  bp <- bp  + theme(legend.position = "none")
  
  #remove gray background and add a border around the plot
  bp <- bp + theme(panel.border = element_rect(fill=NA,color="black", size=0.5,
                                                   linetype="solid"))
  
  
  bp <- bp + theme_bw()
  
  return(bp)
}

# arrange by time
bp <- bp %>% group_by(subject_id) %>% arrange(charttime_conv)

m <- CalcEnSlide(df = bp, 
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

plot_bp(14244) #high ent - mortality
plot_bp(11952) # low entropy - mortality - not a great looking plot....

# Select the firest 200 measurements for each patients
bp_200 <- bp %>% group_by(subject_id) %>% 
  arrange(charttime_conv) %>% 
  slice(1:200) 

# check counts 
View(bp_200 %>% group_by(subject_id) %>% dplyr::summarize(count = n()))
nrow(unique(data.frame(bp_200$subject_id)))

# Determine the max, min, mean, sd, and variance of the first 200 measurements
library(dplyr) #reload, interferring with plyr
bp_200 <- bp_200 %>% group_by(subject_id) %>% 
  dplyr::mutate(min_bp = min(value, na.rm = TRUE)) %>% 
  dplyr::mutate(max_bp = max(value, na.rm = TRUE)) %>% 
  dplyr::mutate(mean_bp = mean(value, na.rm = TRUE)) %>%
  dplyr::mutate(sd_bp = sd(value, na.rm = TRUE)) %>% 
  dplyr::mutate(var_bp = var(value, na.rm = TRUE)) %>% 
  select(subject_id, min_bp, max_bp, mean_bp, sd_bp, var_bp) %>%
  slice(1L)

# format entropy data
m_ent <- m %>% select(subject_id, entropy)
colnames(m_ent)[2] <- "entropy_bp"

# add entropy + demographics
bp_200 <- merge(bp_200, m_ent, by = "subject_id")

# save as csv
write.csv(bp_200, "/Users/User/Box Sync/Projects/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/bp_200.csv", row.names = FALSE)
