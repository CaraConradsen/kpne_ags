# -------------------------------
# Base R Weekly Planner (Top X-axis, Borders, Day numbers)
# -------------------------------

# project data
project_blocks <- fread("C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/analysis/workflow_diagram/ms_plan.csv")
colnames(project_blocks)[1] = "day"
project_blocks$day = as.Date(project_blocks$day)

# Define plotting range
start_date <- as.Date("2025-10-09")  # Monday
end_date   <- as.Date("2025-11-28")  # ~8 weeks
all_days <- seq(start_date, end_date, by="day")

# Map weekdays (Mon=1..Sun=7)
wday_num <- as.POSIXlt(all_days)$wday
wday_num <- ifelse(wday_num==0, 7, wday_num)

# Week numbers for Y-axis
week_num <- as.numeric(factor(format(all_days, "%U")))

# Set up empty plot
par(mar=c(1,0.5,0.5,0.5))
plot(NA, xlim=c(0.5, 6.25), ylim=c(max(week_num)+0.2,0.5),
     xaxt="n", yaxt="n", xlab="", ylab="", main="",
     bty="n")


# Grey background for Sat/Sun
rect(xleft=6-0.5, xright=7+0.5, ybottom=0, ytop=max(week_num)+1,
     col="gray99", border=NA)

# Grid lines
abline(h=0:max(week_num)+0.5, col="gray50", lty=1)
abline(v=1:7-0.5, col="gray50", lty=1)
box()

# Pre-calculate week numbers for all days
week_lookup <- data.frame(
  day = all_days,
  week = as.numeric(factor(format(all_days, "%U")))
)

# X-axis on top
axis(3, at=1:7, labels=c("Mon","Tue","Wed","Thu","Fri","Goals/Notes",""),
     pos = 0.6, tick = FALSE)

# Y-axis: weeks on left
axis(2, at=unique(week_num), pos = 0.55, tick = FALSE,
     labels=unique(week_num), las=1)

# Draw AM/PM boxes with solid borders and day numbers
for(i in 1:nrow(project_blocks)){
  day <- project_blocks$day[i]
  period <- project_blocks$period[i]
  colr <- project_blocks$pcolour[i]
  
  wday <- as.POSIXlt(day)$wday
  wday <- ifelse(wday==0, 7, wday)
  
  # Look up the correct week
  week <- week_lookup$week[week_lookup$day == day]
  
  # AM lower half, PM upper half
  ybottom <- week - ifelse(period=="AM", 0.4, 0)
  ytop    <- week - ifelse(period=="AM", 0, -0.5)
  
  # Draw rectangle with solid border
  rect(xleft=wday-0.49, xright=wday+0.49, ybottom=ybottom, ytop=ytop,
       col=colr, border= NA, lwd=1.5)
  
  # Wrap text with tighter spacing
  lines <- strwrap(project_blocks$label[i], width = 30)
  n_lines <- length(lines)
  
  # Adjust vertical spacing manually (smaller = tighter)
  # line_spacing <- 0.1  # try 0.02 or 0.015 for even tighter
  line_spacing <- 0.1 * (par("usr")[4] - par("usr")[3]) / length(unique(week_num))
  
  mid_y <- (ybottom + ytop) / 2
  start_y <- mid_y + (n_lines / 2 - 0.5) * line_spacing
  
  for (j in seq_along(lines)) {
    y_pos <- start_y - (j - 1) * line_spacing
    text(x = wday - 0.45, y = y_pos, labels = lines[j], cex = 0.65, adj = c(0, 0.5))
  }
  
  
  # Add day number in top-right corner
  rect_width = 0.5 
    margin =  0.01
  
  day_num <- format(day, "%d")
  if (period=="AM"){
    text(x=wday + rect_width - margin, y=ytop + 0.35,
       labels=day_num, cex=0.6, adj=c(1,-8))
    }

}
