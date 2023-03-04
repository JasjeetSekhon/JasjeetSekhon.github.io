### PS 236 Section, October 2008   ###
### Plotting                       ###


## Load the following .RData file from my website
setwd("C:/Users/Charles/R")
load(file="education.RData")

## Let's get a summary of the data
summary(X)

## Let's get a picture of the distribution of incomes
hist(X$inctot)

## Ugh! Close this graph
dev.off()

## Can we make it better?
## See ?par and ?plot for more information
plot(
  # Plot a histogram object with 51 bins
  hist(X$inctot,breaks=50),
  # Use probabilities rather than gross numbers
  freq=FALSE,
  # Turn off the axes
  xaxt="n", yaxt="n",
  # Give the x-axis a label
  xlab="Individual income",
  # Give the y-axis a label
  ylab="Frequency",
  # Give it a title (with a break using \n)
  main="Histogram of income \n 2000 Census, 0.1% Sample")
  # Add an x (1) axis
  axis(1,
    # With tick marks at intervals of $200K between $0 and $1M
    at=seq(0,1000000, 200000),
    # Label the tick marks
    labels=c("0", "200K", "400K", "600K", "800K", "1M"))
    
# What if we want a smoothed representation
plot(
  # Plot the density of the income distribution
  density(X$inctot),
  xaxt="n", yaxt="n",
  xlab="Individual income",
  ylab="Frequency",
  main="Histogram of income \n 2000 Census, 0.1% Sample")
  axis(1,
    at=seq(0,1000000, 200000),
    labels=c("0", "200K", "400K", "600K", "800K", "1M"))
    
## Suppose that we care about just the middle class, defined as those with
##  income between the first and third quartiles
quantile(X$inctot)

MC <- X[X$inctot > 8400 & X$inctot < 86275,]

plot(density(MC$inctot),
  # Set the limits of the x-axis to correspond with our sample limiting
  xlim=c(8400,86275),
  xaxt="n", yaxt="n",
  xlab="Individual income",
  ylab="Frequency",
  main="Histogram of middle class income \n (individuals between first and third quartiles)")
  axis(1, at=seq(10000,85000,25000),
    labels=c("10K", "35K", "60K", "85K"))
    
## Let's see the density of incomes for whites overlayed on that for blacks
plot(density(MC$inctot[MC$race=="Black/Negro"]),
  xlim=c(8400,86275),
  xaxt="n", yaxt="n",
  xlab="Individual income",
  ylab="Frequency",
  main="Histogram of middle class income \n (individuals between first and third quartiles)",
  sub="Density for whites in red, for blacks in black")
  axis(1, at=seq(10000,85000,25000),
    labels=c("10K", "35K", "60K", "85K"))
  lines(density(MC$inctot[MC$race=="White"]), col="red")
  
## But we print in black and white! Let's use a dashed line type for whites
plot(density(MC$inctot[MC$race=="Black/Negro"]),
  xlim=c(8400,86275),
  xaxt="n", yaxt="n",
  xlab="Individual income",
  ylab="Frequency",
  main="Histogram of middle class income \n (individuals between first and third quartiles)",
  sub="Density for whites dashed, solid for blacks")
  axis(1, at=seq(10000,85000,25000),
    labels=c("10K", "35K", "60K", "85K"))
  lines(density(MC$inctot[MC$race=="White"]), lty=2)
  
## Let's have these graphs next to one another
par(
  # Use 'mfrow=c(a,b)' to create an a x b array of graphs filled row-by-row
  # Use 'mfcol' to fill by columns (BTW, 'mf' stands for Multi-Frame)
  mfrow=c(2,2))
  plot(density(MC$inctot[MC$race=="White"]))
  plot(density(MC$inctot[MC$race=="Black/Negro"]))
  ## Plot a histogram for Japanese and Chinese---'|' is an 'or' operator
  plot(density(MC$inctot[MC$race=="Japanese" | MC$race=="Chinese"]))
  plot(density(MC$inctot[MC$race=="American Indian or Alaska Native"]))
  
## Knowing what they look like, let's clean it up by changing titles and
##   using the same scales
par(mfrow=c(2,2))
  plot(density(MC$inctot[MC$race=="White"]),
    xaxt="n", yaxt="n", main="Whites", xlim=c(8400,86275),xlab="", ylab="",
    ylim=c(0,.00004))
      axis(1, at=seq(10000,85000,25000),
        labels=c("10K", "35K", "60K", "85K"))
  plot(density(MC$inctot[MC$race=="Black/Negro"]),
    xaxt="n", yaxt="n", main="Blacks", xlim=c(8400,86275),xlab="", ylab="",
    ylim=c(0,.00004))
      axis(1, at=seq(10000,85000,25000),
        labels=c("10K", "35K", "60K", "85K"))
  plot(density(MC$inctot[MC$race=="Japanese" | MC$race=="Chinese"]),
    xaxt="n", yaxt="n", main="Asian", xlim=c(8400,86275),xlab="", ylab="",
    ylim=c(0,.00004))
      axis(1, at=seq(10000,85000,25000),
        labels=c("10K", "35K", "60K", "85K"))
  plot(density(MC$inctot[MC$race=="American Indian or Alaska Native"]),
    xaxt="n", yaxt="n", main="Amer. Ind/AK Native", xlim=c(8400,86275),xlab="", ylab="",
    ylim=c(0,.00004))
      axis(1, at=seq(10000,85000,25000),
        labels=c("10K", "35K", "60K", "85K"))
        
## What happens if we plot a factor by some continuous variable
plot(MC$educrec, MC$inctot)

## Let's get rid of the empty NA factor
educ <- MC$educrec[,drop=TRUE]

plot(educ, MC$inctot, xaxt="n", yaxt="n",
  xlab="Years of education",
  ylab="Income")
axis(1, at=c(1:9), labels=c("None", "1-4", "5-8", "9", "10", "11", "12",
  "13-15", "16+"))
axis(2, at=seq(10000,85000,15000), labels=c("10K", "25K", "40K", "55K", "70K",
  "85K"))


## Let's get those box plots by race
## Let's define 4 races
race <- ifelse(MC$race == "White", "White", NA)
race <- ifelse(MC$race == "Black/Negro", "Black", race)
race <- ifelse(MC$race=="Japanese" | MC$race=="Chinese", "Asian", race)
race <- ifelse(is.na(race), "Other", race)
race <- as.factor(race)

## Load the 'lattice' package
library(lattice)

## Do a box (and whisker---BW) plot of income by education level for each race
bwplot(educ ~ MC$inctot | race,
  main="Box plot of income by education by race",
  ylab="Years of education", xlab="Income",
  scales=list(x = list(at=seq(10000,85000,25000),
    labels=c("10K", "35K", "60K", "85K")),
    y = list(labels=c("None", "1-4", "5-8", "9", "10", "11", "12",
    "13-15", "16+"))))
    
## Let's make 'educ' numeric
years <- ifelse(educ == "None or preschool", 0, NA)
years <- ifelse(educ == "Grade 1, 2, 3, or 4", 2.5, years)
years <- ifelse(educ == "Grade 5, 6, 7, or 8", 6.5, years)
years <- ifelse(educ == "Grade 9", 9, years)
years <- ifelse(educ == "Grade 10", 10, years)
years <- ifelse(educ == "Grade 11", 11, years)
years <- ifelse(educ == "Grade 12", 12, years)
years <- ifelse(educ == "1 to 3 years of college", 14, years)
years <- ifelse(educ == "4+ years of college", 16, years)

## Let's do a regression of income on education
reg <- lm(MC$inctot ~ years)

## Let's make a scatterplot (not the best data setup for this---sorry)
plot(years ~ MC$inctot)
  ## Let's add the regression line
  abline(reg)
  
## QQ-plots compare the distributions of either:
##   1. Two empirical distributions 'qqplot(x,y)'
##   2. An empirical distribution against a theoretical distribution,
##      e.g., normal 'qqnorm(x)'
## Let's do the latter and compare the distribution of income to the normal
##   distribution. This lets us assess how closely the distribution of income
##   is to normal.
## Before running it: What do you think it will look like (the normal
##   distribution is on the x-axis and income is on the y-axis)
## First, recall the density of income
plot(density(MC$inctot))

qqnorm(MC$inctot)
  ## Add a line showing equal quantile distributions
  qqline(MC$inctot)
  
## It's even starker using the full-sample incomes
qqnorm(X$inctot)
  qqline(X$inctot)
  
## Let's save a plot
## If you are using LaTeX, it's best to save a PostScript format file
## Start the PostScript driver and name the file
postscript("plot.eps",
  # Specify the width and the height of the image that you are saving
  width=5, height=5,
  # The following line uses the same fonts used in LaTeX
  encoding="TeXtext.enc", family="ComputerModern")
  
## For non-LaTeXers, the best way to save an image is via PDF
# pdf(file="plot.pdf")

## Now, whatever you put on your graphics device will be saved.
bwplot(educ ~ MC$inctot | race,
  main="Box plot of income by education by race",
  ylab="Years of education", xlab="Income",
  scales=list(x = list(at=seq(10000,85000,25000),
    labels=c("10K", "35K", "60K", "85K")),
    y = list(labels=c("None", "1-4", "5-8", "9", "10", "11", "12",
    "13-15", "16+"))))
## Once you've put everything on the graph, turn the device and driver off
dev.off()