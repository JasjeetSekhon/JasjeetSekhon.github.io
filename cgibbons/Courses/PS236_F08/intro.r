### PS 236 Section, September 2008 ###
### Basic R code ###

## We can create our own data
## There are a few ways to store a variable
a <- 10
a
b = 6
b

## Let's make some vectors
## c() joins the objects inside to make a vector
a <- c(1,2,3,4,5,6,7,8,9,10)
a

## Here's an easier way to accomplish the same thing
a1 <- c(1:10)
a1

## How long is our vector?
length(a)

## Note that this does NOT work for vectors:
dim(a)

## Are they the same?
a == a1
sum(a == a1)

## Here's another way
b <- seq(1,10)
b

## We could count by 2's
b2 <- seq(1,10,2)
b2

## Does a equal b2?
a == b2

## But are all the numbers in b2 also in a?
b2 %in% a

## What percent of a is in b2?
sum(a %in% b2)/length(a)

## Let's multiply a and a1
a * a1
## This is like the dot product in matrix algebra

## What happens here? Remember, b2 is shorter than a
a
b2
a * b2
## b2 is recycled until it is the same length as a

## We can create an identity vector
d <- rep(1,10)
d

## What if we want just one element from the vector?
a[5]

## Or a few?
b[c(2:5)]

## What about the last five?
b[ (length(b) - 5):length(b) ]

## Let's use a STATA file
##   A 3% random sample of the 2006 American Community Survey by the Census
##   Bureau from the IPUMS database
## First, load the 'foreign' package
## If you don't have this package, R should download it for you
library(foreign)

## Now, set the working directory --- note the direction of the slashes
## BE SURE TO SET TO YOUR COMPUTER'S WORKING DIRECTORY (where you stored the
##   data file).
setwd("C:/Users/Charles/R")

## Load the data and store it as 'X' ##
X <- read.dta("ownership.dta")

## What if we had forgotten how to use 'read.dta'?
?read.dta
help(read.dta)

## Let's see what we have
dim(X)

## So how many variables do our data have?
dim(X)[2]

## How many observations?
dim(X)[1]

## What do the data look like?
summary(X)

## Notice anything?

## The year is 2006 for every observation, so let's delete that variable
X <- X[ , -c(1)]
## Note that we select observations (rows) before the comma and variables
##   (columns) after the comma
## Alternatively, we could have written:
##   X <- X[ , c(2:5)]
##   X <- X[ , c(2:dim(X)[2])]

## Let's name what we have
state <- X[,1]    # The state
own <- X[,2]      # An indicator of home ownership
income <- X[,3]   # Household income
povLine <- X[,4]  # Income as a percentage of the poverty line

## How many households have incomes below 0?
sum(income < 0)

## Let's set all incomes below 0 to 0
income <- ifelse(income < 0, 0, income)
## First, you provide a logical statement, then provide the value of the
##   variable when that statement is true, followed by the value when it is not
sum(income < 0)

## Let's make a new variable indicating whether a household's income is above
##   or below the poverty line
poverty <- ifelse(povLine < 100, 1, 0)

## What percentage of families have incomes below the poverty line?
sum(poverty)/length(poverty)


## The 'class' of an object greatly impacts what you can do with it
## A 'data frame' is a collection of vectors of anything and of possibly
##   differing lengths
class(X)

## A vector can be numeric
class(income)

## Or integers
class(povLine)

## Let's make a character vector of numbers
char <- c("1", "2", "3", "4", "5")
class(char)

## We can turn char into a numeric vector
char <- as.numeric(char)
char
class(char)

## Or back
char <- as.character(char)
char
class(char)

## A vector can also be a factor
## A factor takes on a limited set of integer values or character expressions
##   that represents something else
class(state)

## What 'levels' or values of the factor are included?
levels(state)

## Let's create a new data frame, Y, that includes only CA observations
Y <- X[state == "California", ]

## Now what are the summary statistics?
summary(Y)

## Let's extract the state FIPS vector from Y
yState <- Y[,1]
yState[1:10]

## What class is yFIPS?
## It should be a factor with only one level, "06"
class(yState)
levels(yState)

## Notice that R does NOT discard unused levels by default!
## This has killed me before
yState <- yState[, drop=TRUE]
levels(yState)

## Let's say that we only want CA data --- we can remove the rest
poverty[1:10]

rm(X, state, own, income, povLine, poverty)

poverty[1:10]

## What's left?
ls()

## Do our variables have names?
dimnames(Y)
## The first set are the names of the rows and the second the names of the
##   columns

## Here, we have a 'list' object type --- it's like a vector, where each element
##   can be a vector itself
## To get the second object of a list, we use the following '[[ ]]' command
dimnames(Y)[[2]]

## Now it's just a vector --- let's get the name of the first column
dimnames(Y)[[2]][1]

## Let's get rid of the state FIPS variable
Y <- Y[, -c(1)]

## And name the rest
dimnames(Y)[[2]] <- c("own", "income", "povRatio")
dimnames(Y)[[2]]

## Let's fix these data as we did above, using a different method

## There is no 'own' variable
own[1:10]

## Let's 'attach' the dim names of Y
attach(Y)

own[1:10]

## Now alter the data
income <- ifelse(income < 0, 0, income)
poverty <- ifelse(povRatio < 100, 1, 0)
ownInd <- ifelse(own %in% c("Owned free and clear", "Owned with mortgage or loan"),
  1, NA)
ownInd <- ifelse(own == "With cash rent", 0, ownInd)
ownInd <- as.factor(ownInd)

detach(Y)

# Let's create a data frame with these variables
CaData <- data.frame(income, poverty, ownInd)

## Let's save CaData in a more accessible format, comma separated values
write.csv(CaData, file="CAdata.csv")






