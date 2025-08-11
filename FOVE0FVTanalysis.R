rm(list = ls())

library(saccadr)
library(conflicted)
library(tidyverse)
library(saccades)
library(dplyr)
library(useful)
library(stats)
library("ggplot2")

#read file with the name
df0 <- read.csv("building1.csv")

# extract left eye data
df1 <- df0[,c(2,26,27,28,51,53)]

l = nrow(df0)
trial <- rep(1,l)
trial_df <- data.frame(trial)

#add frame number
fr_number <- rep(1:nrow(df0),)
fr_number_df <-data.frame(fr_number)

#make clean dataframe
df1 <- mutate(df1, trial_df, fr_number_df)
pre_blink =df1[,c(5,6)]
frame=df1[8]
time =df1[1]
eyeX =df1[2]
eyeY =df1[3]
eyeZ =df1[4]
trial =df1[7]

#convert FOVE eye ray unit vector to x and y coordinates
x <- eyeX/eyeZ
colnames(x) <- "x"
y <- eyeY/eyeZ
colnames(y) <- "y"

##calculate the difference in angle between current and the next frame
OcuRad <- acos(((eyeX)*(lead(eyeX))+(eyeY)*(lead(eyeY))+(eyeZ)*(lead(eyeZ)))/(sqrt((eyeX)^2+(eyeY)^2+(eyeZ)^2)*sqrt((lead(eyeX))^2+(lead(eyeY))^2+(lead(eyeZ))^2)))
OcuDeg <- (OcuRad*180)/pi
colnames(OcuDeg) <- "OcuDeg"
OcuDeg[is.na(OcuDeg)] <- 0
df1 <-mutate(df1,OcuDeg)

colnames(time) <- "time"
df3 <-cbind(frame,time, x, y, trial)

#Extract blink data from FOVE headset
#Change 'NotDetected' to 'Closed'
pre_blink[,1 == "NotDetected"] <- "Closed"
pre_blink[,2 == "NotDetected"] <- "Closed"
#If either left or right eye is detected as "closed", consider blink as "closed"
blink_bool <- rep(FALSE, l)
for(i in 1:l){
  if(pre_blink[i,1]=="Closed"|pre_blink[i,2]=="Closed"){
    blink_bool[i] <- TRUE
  }
}
post_blink <- ifelse(blink_bool, "Closed", "Opened")
blink <- mutate(time, post_blink)

#Consider 7 frames before and 7 frames after "Closed" as "Closed" (the eye data is noisy before and after FOVE detection of blinks)
beforeaft_blink <- rep(FALSE, l)
for(i in 8:(l-8)){
  if(blink$post_blink[i]=="Closed"){
    beforeaft_blink[i-(0:7)] <- TRUE
    beforeaft_blink[i+(0:7)] <- TRUE
  }
}
for(i in 1:7){
  if(blink$post_blink[i]=="Closed"){
    beforeaft_blink[1:i] <- TRUE
    beforeaft_blink[i+(0:7)] <- TRUE
  }
}
for(i in (l-7):l){
  if(blink$post_blink[i]=="Closed"){
    beforeaft_blink[i-(0:7)] <- TRUE
    beforeaft_blink[i:l] <- TRUE
  }
}
blink <- mutate(blink, beforeaft_blink)

#omit the first group if closed first - this is required because the latter script that detects change in eye state only works if the eyes are open from the beginning
blink$closedfirst <- rep(FALSE, l)
if(blink$beforeaft_blink[1]==TRUE){
  for(i in 1:l){
    if(blink$beforeaft_blink[i]==TRUE){
      blink$closedfirst[i] <- TRUE
    }else{
      break
  }
  }
}
#blink_closedfirst shows data when eyes were closed from the beginning of data acquisition. The rest of the data goes into "blink".
blink_closedfirst <- dplyr::filter(blink, blink[,4]==TRUE)
blink <- dplyr::filter(blink, blink[,4]==FALSE)

# create another column with eye state shifted one row downward
blink <- shift.column(data=blink, columns="beforeaft_blink", len=1, up=FALSE)
# TRUE if the shifted eye state is different from the original eye state column (i.e. blink_true=where open/closed boundary is)
blink_true <- blink$beforeaft_blink  != blink$beforeaft_blink.Shifted
blink_true_df <- mutate(blink, blink_true) 

# create a cumulative summation of the TRUE values
blink_true_cumsum <-  mutate(blink_true_df, cumsum(blink_true_df$blink_true))

# group the cumulative summation to show where open/closed changes
blink_group <-blink_true_cumsum %>% 
  group_by(`cumsum(blink_true_df$blink_true)`) %>% 
  summarise_at(vars(beforeaft_blink),
               list(min = min))

# select rows with TRUE with timestamp
blink_onlytrue <- dplyr::filter(blink_true_df, blink_true==TRUE)
blink_row <- nrow(blink_onlytrue)

#create new columns for blink_start and blink_end times. If there is only one blink and stays closed at the last frame, create the last blink
if(blink_row==1 & tail(blink$beforeaft_blink, n=1)==TRUE){
  blink_start <- blink_onlytrue[,1]
  blink_end <- max(blink[,1])
  blink_time <- data.frame(blink_start, blink_end)
  blink_num <- nrow(blink_time)
}else{
  if(blink_row==0){
    blink_num <- 0
  }else{
    b_frame_1 <- data.frame()
    b_frame_2 <- data.frame()
  for(i in 1:blink_row){
    blink_start <- ifelse(blink_onlytrue$beforeaft_blink[i]==TRUE, blink_onlytrue$time[i], "")
    blink_end <- ifelse(blink_onlytrue$beforeaft_blink[i]==FALSE, blink_onlytrue$time[i], "")
    b_frame_1 <- data.frame(blink_start, blink_end)
    b_frame_2 <- rbind(b_frame_2, b_frame_1,  factor.exclude = NA)}
  blink_onlytrue %>% cbind(blink_onlytrue, b_frame_2)
  #allocate blink numbers and blink time with original timestamp
  blink_end_shifted <- shift.column(data=b_frame_2, columns="blink_end", len=1, up=TRUE) %>% dplyr::select(-c("blink_end"))
  #get rid of blank spaces in between rows created after shifting
  is_blank <- function(x) {is.na(x)|x==""}
  unnecessary_row <- apply(blink_end_shifted, 1, function(x){
    all(is_blank(x))
  })
  blink_time <- blink_end_shifted[!unnecessary_row,]
  rownames(blink_time) <- c(1:nrow(blink_time))
  # calculate the number of blinks
  blink_num <- nrow(blink_time)
}
}
if(nrow(blink_closedfirst)>0){
  if(blink_num==0){
    blink_time <- data.frame()
  }
  closedfirst_time1 <- blink_closedfirst$time[1]
  closedfirst_time2 <- max(blink_closedfirst$time)
  closedfirst_row <- data.frame(closedfirst_time1, closedfirst_time2)
  colnames(closedfirst_row) <- c("blink_start", "blink_end.Shifted")
  blink_time <- rbind(closedfirst_row, blink_time)
}

#identify rows with y values of over 0.56 or under -0.56 and x values of over 1.13 or under -1.13 (as it is very unlikely to look outside the whole of VR visual area
outside_visual <- dplyr::filter(df3, df3[,4]>=0.56|df3[,4]<=-0.56|df3[,3]>=1.13|df3[,3]<=-1.13)

##detect fixation using R package
fixation <- detect.fixations(df3)
fixation_only <- subset(detect.fixations(df3), event=="fixation")

#Identify fixations that clash with blinks part1
if(nrow(blink_onlytrue)>0){
fixation_only$filter<-rep(FALSE,nrow(fixation_only))
for(i in 1:nrow(fixation_only)){
  for(j in 1:blink_num){
    if(fixation_only[i,2]>=blink_time[j,1]&fixation_only[i,2]<=blink_time[j,2]){
      fixation_only$filter[i] <- TRUE
    }else{
      if(fixation_only[i,3]>=blink_time[j,1]&fixation_only[i,3]<=blink_time[j,2]){
        fixation_only$filter[i] <- TRUE
      }else{
        if(fixation_only[i,2]<=blink_time[j,1]&fixation_only[i,3]>=blink_time[j,2]
           ){
          fixation_only$filter[i] <- TRUE
        }
      }
    }
  }
}
fixation_filtered <- dplyr::filter(fixation_only, fixation_only$filter == FALSE)
}else{
  fixation_filtered <- fixation_only
}

#Identify fixations that clash with blinks part2
if(nrow(blink_closedfirst>0)){
  fixation_filtered <- dplyr::filter(fixation_filtered, fixation_filtered[,2]>max(blink_closedfirst[,1]))
}

#identify fixations that are outside range (see line 165)
if(nrow(outside_visual)>0){
  fixation_filtered$filter_bool <- rep(FALSE, nrow(fixation_filtered))
  for(i in 1:nrow(fixation_filtered)){
    for(j in 1:nrow(outside_visual)){
      if(fixation_filtered[i,2]>=outside_visual[j,2]&fixation_filtered[i,3]<=outside_visual[j,2]){
        fixation_filtered$filter_bool[i] <- TRUE
      }else{
        if(fixation_filtered[i,2]==outside_visual[j,2]|fixation_filtered[i,3]==outside_visual[j,2]){
          fixation_filtered$filter_bool[i] <- TRUE
        }
      }
    }
  }
  fixation_filtered <- dplyr::filter(fixation_filtered, fixation_filtered$filter_bool == FALSE)
}

#calculate number of fixations
fixation_number_row<- rep(1:nrow(fixation_filtered),)
fixation_filtered <- mutate(fixation_filtered,fixation_number_row)
fixation_number <- nrow(fixation_filtered)

# calculate duration of fixation
fix_frame <- (fixation_filtered[3])-(fixation_filtered[2])
fixation_duration <- apply(fix_frame,2,mean)

#extract Ocudeg from fixation
frame1 <-data.frame()
frame2 <-data.frame()
for(i in 1:nrow(fixation_filtered)){
  #extract start and end of fixation
  fix_start = as.numeric(fixation_filtered[i,2])
  fix_end = as.numeric(fixation_filtered[i,3])
  fixation_frame <- df1 %>% dplyr::filter(.[[1]] >= fix_start, .[[1]]<= fix_end)
  #extract three dimensional ey eray x y z of fixation
  fix_x <- mean(fixation_frame[,2], na.rm=TRUE)
  fix_y <- mean(fixation_frame[,3], na.rm=TRUE)
  fix_z <- mean(fixation_frame[,4], na.rm=TRUE)
  frame1 <- data.frame(fix_x, fix_y, fix_z)
  frame2 <- rbind(frame2,frame1)
}
fixation_filtered <- cbind(fixation_filtered, frame2)
fixation_filtered <- fixation_filtered[,c("fixation_number_row","start","end","dur","fix_x","fix_y","fix_z","x","y")]

#extract the position of fixation
frame5 <-data.frame()
frame6 <-data.frame()
for(i in 1:nrow(fixation_filtered)){
  fix_start = as.numeric(fixation_filtered[i,2])
  fix_end = as.numeric(fixation_filtered[i,3])
  fixation_frame <- df3 %>% dplyr::filter(.[[2]] >= fix_start, .[[2]]<= fix_end)
  #identify all fixation points
  x_series <- fixation_frame[,3]
  y_series <- fixation_frame[,4]
  time_series <- fixation_frame[,2]
  frame5 <- data.frame(time_series, x_series, y_series)
  frame6 <- rbind(frame6,frame5)
}
fixation_points <- frame6

#Calculate angle between current and next fixation points
fix_x=fixation_filtered[5]
fix_y=fixation_filtered[6]
fix_z=fixation_filtered[7]
fix_OcuRad <- acos(((fix_x)*(lead(fix_x))+(fix_y)*(lead(fix_y))+(fix_z)*(lead(fix_z)))/(sqrt((fix_x)^2+(fix_y)^2+(fix_z)^2)*sqrt((lead(fix_x))^2+(lead(fix_y))^2+(lead(fix_z))^2)))
fix_OcuDeg <- (fix_OcuRad*180)/pi
colnames(fix_OcuDeg) <- "OcuDeg"
data_fixation <- cbind (fixation_filtered, fix_OcuDeg)

#Scanpath_length
Scanpath_length = apply(fix_OcuDeg, 2, sum, na.rm =TRUE)

#Saccade
#convert x and y to angle in degrees
angx <- (-(asin(eyeX/(sqrt(eyeX^2+eyeZ^2)))))*180/pi
colnames(angx) <- "angx"
angy <- (-(asin(eyeY/(sqrt(eyeY^2+eyeZ^2)))))*180/pi
colnames(angy) <- "angy"
saccade <- extract_saccades(angx,angy,70,methods = list(method_ek),velocity_function = saccadr::diff_ek,)
#omit those with duration of zero
saccade <- saccade %>% dplyr::filter(Duration > 0)

#extract ocudeg from saccade
frame3 <-data.frame()
frame4 <-data.frame()
for(i in 1:nrow(saccade)){
  sac_start = as.numeric(saccade[i,3])
  sac_end = as.numeric(saccade[i,4])
  sac_fr_length = (sac_end)-(sac_start)
  saccade_frame <- df1 %>% dplyr::filter(.[[8]] >= sac_start, .[[8]]<= sac_end)
  #match the saccade start and end time using FOVE app timestamp
  sac_apptime <- saccade_frame[,1]
  sac_start_apptime <- min(sac_apptime)
  sac_end_apptime <- max(sac_apptime)
  sac_amplitude <- saccade[i,14]
  sac_duration <- sac_end_apptime-sac_start_apptime
  sac_vel_average <- (sac_amplitude)/(sac_duration)
  colnames(sac_vel_average) <- "sac_vel_average"
  frame3 <- data.frame(sac_start, sac_end, sac_start_apptime, sac_end_apptime, sac_duration, sac_amplitude, sac_vel_average)
  frame4 <- rbind(frame4,frame3)
}

#Identify saccades that clash with blinks part1
if(nrow(blink_onlytrue)>0){
  frame4$filter<-rep(FALSE,nrow(frame4))
  for(i in 1:nrow(frame4)){
    for(j in 1:blink_num){
      if(frame4[i,3]>=blink_time[j,1]&frame4[i,3]<=blink_time[j,2]){
        frame4$filter[i] <- TRUE
      }else{
        if(frame4[i,4]>=blink_time[j,1]&frame4[i,4]<=blink_time[j,2]){
          frame4$filter[i] <- TRUE
        }else{
          if(frame4[i,3]<=blink_time[j,1]&frame4[i,4]>=blink_time[j,2]
          ){
            frame4$filter[i] <- TRUE
          }
        }
      }
    }
  }
  #Omit unnecessary saccades
  saccade_filtered1 <- dplyr::filter(frame4, frame4$filter == FALSE)
}else{
  saccade_filtered1 <- frame4
}

#Identify saccades that clash with blinks part2
if(nrow(blink_closedfirst>0)){
  saccade_filtered2 <- dplyr::filter(saccade_filtered1, saccade_filtered1[,3]>max(blink_closedfirst[,1]))
}else{
  saccade_filtered2 <- saccade_filtered1
}

#identify saccades that  are outside range (see line 165) and omit
if(nrow(outside_visual)>0){
  saccade_filtered2$filter_bool <- rep(FALSE, nrow(saccade_filtered2))
  for(i in 1:nrow(saccade_filtered2)){
    for(j in 1:nrow(outside_visual)){
      if(saccade_filtered2[i,3]>=outside_visual[j,2]&saccade_filtered2[i,4]<=outside_visual[j,2]){
      saccade_filtered2$filter_bool[i] <- TRUE
    }else{
    if(saccade_filtered2[i,3]==outside_visual[j,2]|saccade_filtered2[i,4]==outside_visual[j,2]){
      saccade_filtered2$filter_bool[i] <- TRUE
    }
    }
    }
  }
  saccade_filtered2 <- dplyr::filter(saccade_filtered2, saccade_filtered2$filter_bool == FALSE)
}
  
#check data
xy <-data.frame(x,y,angx,angy)
data_check <- mutate(df1, xy)

#plot saccade xy
frame_xy1_start <- data.frame()
frame_xy1_end <- data.frame()
frame_xy2 <- data.frame()
for(i in 1:nrow(saccade_filtered2)){
  frame_xy1_start <- dplyr::filter(data_check, saccade_filtered2[i,3]==data_check$Application.Time)
  frame_xy1_end <- dplyr::filter(data_check, saccade_filtered2[i,4]==data_check$Application.Time)
  frame_xy2 <- rbind(frame_xy2, frame_xy1_start, frame_xy1_end)
}
saccade_graph <- plot(frame_xy2$x, frame_xy2$y, type = "l", xlab="X", ylab="Y", xlim=c(-1.5, 1.5), ylim=c(-1, 1))

#calculate mean saccade duration
saccade_duration <- mean(saccade_filtered2[,5])
#calculate mean saccade amplitude
saccade_amplitude_mean <- mean(saccade_filtered2[,6])
#calculate mean saccade velocity
saccade_velocity_mean <- mean(saccade_filtered2[,7])

#calculate number of saccades
saccade_number <- nrow(saccade_filtered2)

data_saccade <- saccade_filtered2


#hull area
library(grDevices)
library(sp)

#find convex hull of points
frame_xy3 <- frame_xy2[,c(1,10,11)]
colnames(fixation_points) <- c("hull_time", "hull_x", "hull_y")
colnames(frame_xy3) <- c("hull_time", "hull_x", "hull_y")
df4 <- rbind(fixation_points, frame_xy3)
df4_ordered <- df4[order(df4$hull_time),]
df4_ordered %>% dplyr::distinct()
#extract x and y and put them into one data frame called 'hullxy'
hullx <- df4_ordered[,2]
hully <- df4_ordered[,3]
hullxy <- data.frame("x"=hullx, "y"=hully)
hull_indices <- chull(hullxy)
# Add the first point to the end to close the polygon
hull_indices <- c(hull_indices, hull_indices[1])
#Plot the points
pdf("graphhull_building1.pdf")
graph2 <- plot(hullxy, cex=0.5,  xlab="X", ylab="Y", xlim=c(-1.5, 1.5), ylim=c(-1, 1))
lines(hullxy[hull_indices,], col = "red", lwd = 2)
dev.off()
#Extract the coordinates of the convex hull points
hull_points <- hullxy[hull_indices,]
# Calculate the area using the Shoelace formula
hullarea <- 0.5 * abs(sum(hull_points[-nrow(hull_points), 1] * hull_points[-1, 2] - hull_points[-1, 1] * hull_points[-nrow(hull_points), 2]))

#Calculate BCEA
sigma_x <- sd(fixation_filtered$x)
sigma_y <- sd(fixation_filtered$y)
rho <- cor(fixation_filtered$x, fixation_filtered$y, method="pearson")
# Degrees of freedom (2 for a 2D ellipse)
dof <- 2
# Probability corresponding to 1 standard deviation (68.2% or 0.682)
probability <- 0.682
# Calculate the chi-square critical value
chi_square_value <- qchisq(probability, dof)
#BCEA
BCEA <- chi_square_value * pi * sigma_x * sigma_y * sqrt(1 - rho^2)


#summarize data
summary <- data.frame(fixation_number, fixation_duration,  saccade_number, saccade_duration, saccade_amplitude_mean, saccade_velocity_mean, Scanpath_length, blink_num, hullarea, BCEA)
rownames(summary)<- NULL

library(openxlsx)
output_set <- list("summary"=summary, "fixation"=data_fixation, "check"=data_check, "saccade"=data_saccade)
write.xlsx(output_set, "result10.20_building1.xlsx")