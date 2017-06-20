
library(clusterSim)
library(ggplot2)
library(reshape2)

##  Information on R's hypergeometric functions
#?Hypergeometric

#helper function to aggregate successes 
sums <- function(samples, df){
  sums<- NULL
  for (i in 1:(ncol(df))){
    res <- c()
    for (j in 0:samples){
      res <- c(res, sum(df[ ,i] == j))
    }
    sums <- cbind(sums, res)
  }
  return(sums)
}


##  hg_sims simulates a draw of size @sample from a hypergeometric population of size @N @trials times
##  for each combination of successes and failures in @N (i.e. if @N = 10 then combinations are 
##  (0 successes, 10 failures), (1 successes, 9 failures), ..., (10 successes, 0 failures)).
##  @N - total population size
##  @samples - the number of samples drawn 
##  @trials - the number of times to perform independent random draws
##  @returns a dataframe witn N+1 columns (b/c 0 is included) and samples+1 rows (O included)
hg_sims<-function(samples, N, trials){
  df <- NULL
  for (i in 0:N){
    df <- cbind(df, rhyper(trials, i, N-i, samples))
  }
  return(as.data.frame(df))
}


##  plotting function
##  df is a dataframe of normalized aggregate successes from simulation
pdf_plot <- function(df, N){
  sample_size = ncol(df)-1
  test <- df
  test['successes'] <- 0:N
  melted = melt(test, id.vars="successes")
  ggplot() + 
    geom_line(data=melted, aes(x=successes, y=value, group=variable, color = variable), size=1) +
    labs(title=paste("Sample Size of", sample_size, " "), x =paste("Successes out of ", toString(N), sep=" "), y = "Probability") +
    guides(color=guide_legend(title="Successes\nin Sample"))
}

draw_dist_plot <-function(successes, df, N){
  df <- df[c(successes, 'count')]
  melted <- melt(df, id.vars='count')
  print(melted)

  ggplot(data=melted, aes(x=count, y=value, fill=variable)) +
    geom_bar(stat="identity", position=position_dodge()) +
    labs(title=paste(toString(successes), ' successes out of ', N), x='Successes in sample', y='Density') +
    stat_smooth(aes(color=variable), method='auto', se = FALSE)
  
}

## Set desired parameters
samples <- 6
N <- 50
trials <- 50000

##  run the simulation
raw_data = hg_sims(samples, N, trials)

##  aggregate the data
aggregated_data = sums(samples, raw_data)
colnames(aggregated_data) <- 0:N
rownames(aggregated_data) <- 0:samples

##  normalize the results by successes in N to get a pdf
agg_data_normed <- data.Normalization(aggregated_data, type = "n10", normalization = "row")

##  plot the pdf's of the draw successes.  Data parameter is transposed and ensured to be a dataframe
pdf_plot(as.data.frame(t(agg_data_normed)), N)


##  Plot entire draw for x successes in N
x=c(25)
##  first normalize by sample to get pmf and add successes column
samples_normed <- as.data.frame(data.Normalization(aggregated_data, type = "n10", normalization = "col"))
samples_normed['count'] = c(0:(nrow(samples_normed)-1))
##  plot the distribution over successes in draws for specified successes in N
draw_dist_plot(x, samples_normed, N)
##  print probability successes out of N greater then specified # of successes (i.e. prob )



right_tail_plot <- function(df, sample_successes, pop_successes){
  sample_successes <- sample_successes +1
  df1 <- as.data.frame(t(df))
  df1['successes'] <- c(0:(nrow(df1)-1))
  shade <- df1[(pop_successes+1):nrow(df1),]
  print(sum(shade[,sample_successes]))
  ggplot(data=df1, aes(x=successes)) + geom_line(aes(y=df1[,sample_successes]), color='blue') + 
    geom_area(data=shade, aes(x=pop_successes:(nrow(df1)-1), y=shade[,sample_successes]), fill='blue') +
    labs(title=paste(sample_successes,' successes in sample of', toString(nrow(df)-1)), x='Successes in population', 
         y='Probability')
}

sample_successes <- 3
pop_successes <- 25

right_tail_plot(agg_data_normed, sample_successes, pop_successes)

?Hypergeometric
for (i in c(0:6)){
  print(i)
  print(dhyper(i, 25, 25, 6))
}
for (i in c(0:6)){
  print(i)
  print(phyper(i, 25, 25, 6))
}
