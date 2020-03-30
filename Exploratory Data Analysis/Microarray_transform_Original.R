library(affy)
library(affyQCReport)
#downloading a data set 
#Downloading a micro array data set and set it as you directory
setwd("/Users/shreeyagarg/Desktop/Files/GSE29797_RAW")
info = ReadAffy()
info
gdata = as.data.frame(justRMA())
gdata = t(gdata)
summary(info) #summary of data
dim(info) #dimensions of data
View(info) #view data
str(info)
class(info) # packages involved
typeof(info)
#Presentation of data
image(info[,1])
hist(info,main="Histogram") # histogram of data
## EXPLORARTORY DATA ANALYSIS-
#QQ PLOT FOR G DATA as exploratory data analysis
qqnorm(gdata,main="QQ Plot")
qqline(gdata,main="QQ LINE PLOT",col='red')
#Clustering analysis 
wss <- (nrow(gdata)-1)*sum(apply(gdata,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(gdata, 
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

# K-Means Cluster Analysis
fit <- kmeans(gdata, 5) # 5 cluster solution
# get cluster means 
aggregate(gdata,by=list(fit$cluster),FUN=mean)
# append cluster assignment
mydata <- data.frame(gdata, fit$cluster)
#exploratory data analysis of correlation plot
correlationPlot(info)
#exploratory data analysis of boxplot
RNA_deg = AffyRNAdeg(info)
View(summaryAffyRNAdeg(RNA_deg))
#Transformation
#data normalisation
rma = rma(info)
matexp = exprs(rma)
par(mfrow = c(1,2))
#PreProcessing and Transformation
boxplot(log2(exprs(info)), main = "Before normalization",
        ylab = "log2(Intensity)")
boxplot(matexp, main = "After normalization")
#exploratory data analysis OF rna Degradation Plot
plotAffyRNAdeg(RNA_deg)
par(mfrow = c(1,1))
## T test-Data Analysis
dataset1 = t(gdata[1:7,1])
dataset2 = t(gdata[8:16,1])
t.test.gene.1 = t.test(dataset1, dataset2, "two.sided")

## EXPLORATORY DATA ANALYSIS
## Scatterplots
setwd("/Users/shreeyagarg/Desktop/Files")
result <- read.table("FPKM_COUNT.txt", header=TRUE)
head(result)
plot(result,main="Scatter Plot") #Scatter Plot
#Scatter Plot by taking combination of two datas one by one
plot(result$FPKM,main="Scatter Plot of FPKM")
plot(result$count, result$FPKM,main="Scatter Plot of FPKM vs count")
pairs(result, panel = panel.smooth,main="All Scatter Plots") # all scatter plots
## Density Plot and Barplot
barplot(result$count, horiz = F,main="BarPlot" ,ylab = "count") # Barplot Vertical
barplot(result$count, horiz = T,main="Horizontal BarPlot",ylab = "count") #Barplot horizontal
plot(density(result$FPKM, na.rm = T)) #Density Plot
rug(result$FPKM) # desity plot with rug
## Line plots
plot(result$count, type = 'l',main="Line Plot Continuous") #continuous line plot
plot(result$count, type = 'b',main="Line Plot with dots") #line plot with dots


