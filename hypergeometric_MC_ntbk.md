Hypergeometric Monte Carlo
================

This is a Monte Carlo sampling of the hypergeometric distribution

First load the libraries we will use

``` r
library(clusterSim)
library(ggplot2)
library(reshape2)
```

Information on R's hypergeometric distribution function can be found using this command

``` r
#?Hypergeometric
```

A function to Monte Carlo sample from a hypergeometric distribution

``` r
##  hg_mc - a function that monte carlo samples from a hypergeometric distribution
##  For each x in N, where N is the population size and x is the true number of successes in N, 
##  this function will draw a random sample of specified size (@sample) a specified number of times(@trials). 
##  @N - total population size
##  @samples - the number of samples drawn 
##  @trials - the number of times to perform independent random draws
##  @returns a dataframe witn N+1 columns (b/c 0 is included) and @samples+1 rows (O included)

hg_mc<-function(samples, N, trials){
  df <- NULL
  for (i in 0:N){
    df <- cbind(df, rhyper(trials, i, N-i, samples))
  }
  return(as.data.frame(df))
}
```

A helper function to aggregate the raw data from the simulations

``` r
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
```

First we set the parameters we want. Here we set a population of 50, sample size of 6, and 10,000 trials.

``` r
## Set desired parameters for the simulation
samples <- 6  # sample size
N <- 50 # population size 
trials <- 10000 # number of trials
```

We now run the simulation and aggregate the raw data

``` r
##  run the simulation
raw_data = hg_mc(samples, N, trials)

##  aggregate the data and rename the rows and cols
aggregated_data = sums(samples, raw_data)
colnames(aggregated_data) <- 0:N
rownames(aggregated_data) <- 0:samples
```

Let's look at the aggregated data. Each row is the observed number of successes in the sample. Each column is the true number of successes in the population.

``` r
print(aggregated_data)
```

    ##       0    1    2    3    4    5    6    7    8    9   10   11   12   13
    ## 0 10000 8752 7731 6780 5907 5180 4445 3830 3355 2814 2362 2036 1720 1474
    ## 1     0 1248 2153 2897 3476 3733 4082 4279 4296 4257 4159 4019 3830 3557
    ## 2     0    0  116  316  581  991 1294 1616 1917 2301 2660 2849 3014 3224
    ## 3     0    0    0    7   36   92  168  257  399  561  712  933 1211 1425
    ## 4     0    0    0    0    0    4   11   17   31   66  102  156  204  287
    ## 5     0    0    0    0    0    0    0    1    2    1    5    7   21   32
    ## 6     0    0    0    0    0    0    0    0    0    0    0    0    0    1
    ##     14   15   16   17   18   19   20   21   22   23   24   25   26   27
    ## 0 1203 1052  847  682  618  484  365  312  226  187  149   97   92   59
    ## 1 3289 3117 2779 2522 2410 2029 1789 1516 1416 1207  946  840  671  586
    ## 2 3379 3345 3473 3568 3389 3453 3260 3244 3006 2818 2555 2402 2124 1956
    ## 3 1672 1893 2145 2289 2473 2656 2928 3006 3185 3308 3450 3329 3408 3314
    ## 4  411  515  655  801  939 1135 1342 1525 1662 1891 2098 2368 2529 2723
    ## 5   44   76   97  131  161  228  299  365  450  510  735  848 1049 1171
    ## 6    2    2    4    7   10   15   17   32   55   79   67  116  127  191
    ##     28   29   30   31   32   33   34   35   36   37   38   39   40   41
    ## 0   44   35   29   26   13    7    2    5    2    1    0    0    0    0
    ## 1  479  330  321  220  186  140   85   81   48   37   23   11    6    3
    ## 2 1710 1539 1340 1144  938  818  664  515  409  271  206  154  101   69
    ## 3 3119 3126 2861 2720 2554 2299 2167 1884 1671 1379 1182  961  757  602
    ## 4 3001 3064 3327 3453 3486 3558 3478 3445 3376 3282 3097 2858 2574 2363
    ## 5 1413 1606 1784 1998 2253 2483 2742 3025 3300 3588 3721 4009 4143 4117
    ## 6  234  300  338  439  570  695  862 1045 1194 1442 1771 2007 2419 2846
    ##     42   43   44   45   46   47   48   49    50
    ## 0    0    0    0    0    0    0    0    0     0
    ## 1    0    1    0    0    0    0    0    0     0
    ## 2   36   25    7    2    1    0    0    0     0
    ## 3  392  232  162   91   46    8    0    0     0
    ## 4 1930 1609 1383  919  623  354  112    0     0
    ## 5 4375 4312 4099 3823 3458 2893 2167 1208     0
    ## 6 3267 3821 4349 5165 5872 6745 7721 8792 10000

Let's create some plots to help visualize the data First we write a function to plot the the distribution of the observed successes in a sample over the true number of successes in the population.

``` r
##  plotting function
##  df is a dataframe of normalized aggregate successes from simulation
pdf_plot <- function(df, N){
  sample_size = ncol(df)-1
  test <- df
  test['successes'] <- 0:N
  melted = melt(test, id.vars="successes")
  ggplot() + 
    geom_line(data=melted, aes(x=successes, y=value, group=variable, color = variable), size=1) +
    labs(title=paste("Sample Size of", sample_size, " "), x =paste("True number of successes in population size ", toString(N), sep=" "), y = "Probability density") +
    guides(color=guide_legend(title="Successes\nin Sample"))
}
```

Call the plotting function

``` r
##  normalize the results by successes in N to get a pdf
agg_data_normed <- data.Normalization(aggregated_data, type = "n10", normalization = "row")

##  plot the pdf's of the draw successes.  Data parameter is transposed and ensured to be a dataframe
pdf_plot(as.data.frame(t(agg_data_normed)), N)
```

![](hypergeometric_MC_ntbk_files/figure-markdown_github/unnamed-chunk-9-1.png) This next plotting function is to examine the distribution of sample successes for a specific number of successes in the population. Multiple population successes are allowed for comparison.

``` r
##  To plot the distribution of successes in samples for a specified number of successes in the population
##  @successes - the number of successes in the population
##  @df - a data frame containing the data
##  @N - the population size
 
draw_dist_plot <-function(successes, df, N){
  sample_size = nrow(df) -1
  df <- df[c(successes, 'count')]
  melted <- melt(df, id.vars='count')
  ggplot(data=melted, aes(x=count, y=value, fill=variable)) +
    geom_bar(stat="identity", position=position_dodge()) +
    labs(title=paste('Probability of successes in sample when true population successes are \n ', toString(successes), ' out of ', N), x=paste('Successes in sample size of ', sample_size), y='Probability Density') +
    stat_smooth(aes(color=variable), method='auto', se = FALSE)
}
```

We first set the population successes that we want to look at. Here we choose 20 and 30.

``` r
##  Plot entire sample distribution for x successes in N. Choose x (can choose more than one to compare)
successes=c(20, 30)
```

We normalize the data and call the plotting function

``` r
##  first normalize by sample to get pmf and add successes column
samples_normed <- as.data.frame(data.Normalization(aggregated_data, type = "n10", normalization = "col"))
samples_normed['count'] = c(0:(nrow(samples_normed)-1))
##  plot the distribution over successes in draws for specified successes in N
draw_dist_plot(successes, samples_normed, N)
```

![](hypergeometric_MC_ntbk_files/figure-markdown_github/unnamed-chunk-12-1.png) If we want to examine a specific number of successes in a sample and find the likelihood of the true number of population successes we can plot the sample and sum the area under the curve. This function does that for the right tail of the distribution.

``` r
right_tail_plot <- function(df, sample_successes, pop_successes){
  sample_successes <- sample_successes +1
  df1 <- as.data.frame(t(df))
  df1['successes'] <- c(0:(nrow(df1)-1))
  shade <- df1[(pop_successes+1):nrow(df1),]
  prob <- 100* sum(shade[,sample_successes])
  ggplot(data=df1, aes(x=successes)) + geom_line(aes(y=df1[,sample_successes]), color='blue') + 
    geom_area(data=shade, aes(x=pop_successes:(nrow(df1)-1), y=shade[,sample_successes]), fill='blue') +
    labs(title=paste('Curve: ', sample_successes-1,' successes in sample of', toString(nrow(df)-1), '\nShaded Area: P[True Population Successes > ',
                     pop_successes, '] = ', prob, '%'), x='Successes in population', 
         y='Probability')
}
```

We then call the function based on the parameters we want to see. Here we are looking at 4 successes in the sample size of 6 where we want to know the probability that the true number of successes in the population is greater than 25

``` r
##  parameters for plotting function
sample_successes <- 4 # number of successes in sample
pop_successes <- 25 # probability that true population successes are greater than this when observing sample_successes

right_tail_plot(agg_data_normed, sample_successes, pop_successes)
```

![](hypergeometric_MC_ntbk_files/figure-markdown_github/unnamed-chunk-14-1.png)
