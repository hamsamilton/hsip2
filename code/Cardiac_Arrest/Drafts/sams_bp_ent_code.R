library(TSEntropies)
library(ggplot2)
library(dplyr)
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


# x = subject_id
plot_sbp <- function(x) { 
  
  # arrange by Date
  df <- sbp %>% filter(subject_id == x) %>% arrange(charttime_conv)
  
  #Title: create a vector that will contain the ID number and can be used to title the    graph by patient ID 
  title <- x
  
  #Create limits to identify the start and end dates of the time series
  first <- first(df$charttime_conv)
  last <- last(df$charttime_conv)
  
  sbp <- df %>% filter(!is.na(value)) %>% select(subject_id, charttime_conv, value) %>%
    ggplot(aes(x = charttime_conv, y = value)) + 
    geom_point(color = "slateblue") + 
    geom_line(color = "slateblue") +
    #geom_hline(yintercept=3) + #could demark dangerous (low or high) sbps
    scale_y_continuous(name = "mmHg") + #define the y-axis limits
    #scale_x_continuous(name = "Time", date_break = "1 hour", date_labels = "%H", limits = c(first, last)) + #, limits = c(-6, 0), breaks = c(-6:0)) +
    scale_x_datetime(name = "Time", date_break = "6 hour", date_labels = "%H", limits = c(first, last)) +
    ggtitle(paste0('ID:', title)) + #add a title to identify patients  
    coord_cartesian(ylim = c(0,300))
  
  #remove legend and x-axis
  sbp <- sbp  + theme(legend.position = "none")
  
  #remove gray background and add a border around the plot
  sbp <- sbp + theme(panel.border = element_rect(fill=NA,color="black", size=0.5,
                                                 linetype="solid"))
  
  
  sbp <- sbp + theme_bw()
  
  return(sbp)
}

# Heart Rate Data
sbp <- read.csv("/Users/samuelhamilton/Downloads/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/Chartevents_51_arterial_sbp.csv")

# Patient IDs
ca_ids <- read.csv("/Users/samuelhamilton/Downloads/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/ca_ids.csv")

sbp_counts <- sbp %>% group_by(hadm_id) %>% dplyr::count(num_obs = n())

sbp_counts_200 <- sbp_counts %>% filter(num_obs >= 210)
sbp <- sbp %>% 
  filter(hadm_id %in% sbp_counts_200$hadm_id)

# CHF Patients
sbp$charttime_conv <- strptime(as.character(sbp$charttime), format = "%Y-%m-%d %H:%M:%S", tz = "EST")
sbp$charttime_conv <- as.POSIXct(sbp$charttime_conv, tz="EST")

sbp$value <- as.numeric(as.character(sbp$value)) 

sbp = sbp %>% inner_join(ca_ids, by = c("subject_id","hadm_id"))

sbp <- sbp %>% filter(!value == "NULL") %>% filter(!is.na(value)) %>% 
  filter(value > 0)


m <- CalcEnSlide(df = sbp, 
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

