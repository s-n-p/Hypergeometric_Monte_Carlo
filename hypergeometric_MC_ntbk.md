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
    ## 0 10000 8796 7781 6726 5968 5159 4486 3818 3310 2768 2406 1970 1781 1414
    ## 1     0 1204 2091 2956 3376 3875 4071 4297 4308 4280 4163 4048 3773 3620
    ## 2     0    0  128  305  621  891 1264 1603 1948 2335 2556 2861 3028 3197
    ## 3     0    0    0   13   35   72  173  269  396  561  777  963 1183 1442
    ## 4     0    0    0    0    0    3    6   13   35   55   91  141  217  300
    ## 5     0    0    0    0    0    0    0    0    3    1    7   17   18   27
    ## 6     0    0    0    0    0    0    0    0    0    0    0    0    0    0
    ##     14   15   16   17   18   19   20   21   22   23   24   25   26   27
    ## 0 1207 1021  841  721  568  472  357  284  230  169  148  117   79   68
    ## 1 3316 3004 2788 2530 2342 2011 1795 1591 1324 1182  976  849  697  628
    ## 2 3420 3438 3511 3449 3508 3414 3256 3095 3011 2847 2589 2392 2152 1959
    ## 3 1631 1916 2094 2349 2486 2758 2886 3093 3159 3188 3364 3326 3255 3239
    ## 4  372  557  669  806  929 1084 1379 1542 1798 1969 2166 2393 2652 2699
    ## 5   52   58   92  141  150  245  291  359  437  582  678  830 1029 1216
    ## 6    2    6    5    4   17   16   36   36   41   63   79   93  136  191
    ##     28   29   30   31   32   33   34   35   36   37   38   39   40   41
    ## 0   41   28   29   13   16    5    6    5    0    1    0    1    0    0
    ## 1  506  369  305  226  180  143   86   47   41   28   26   10    9    2
    ## 2 1703 1528 1285 1141  949  795  647  520  387  306  236  155  100   71
    ## 3 3185 3044 3006 2791 2533 2321 2103 1883 1579 1409 1190  917  777  518
    ## 4 2994 3134 3202 3384 3511 3546 3549 3524 3378 3265 3092 2876 2583 2313
    ## 5 1315 1610 1781 2000 2238 2477 2766 3010 3380 3548 3737 3990 4191 4215
    ## 6  256  287  392  445  573  713  843 1011 1235 1443 1719 2051 2340 2881
    ##     42   43   44   45   46   47   48   49    50
    ## 0    0    0    0    0    0    0    0    0     0
    ## 1    1    1    0    0    0    0    0    0     0
    ## 2   35   13    6    2    0    0    0    0     0
    ## 3  400  294  153   98   45   12    0    0     0
    ## 4 1880 1665 1323  960  635  380  113    0     0
    ## 5 4414 4256 4094 3817 3459 2930 2125 1212     0
    ## 6 3270 3771 4424 5123 5861 6678 7762 8788 10000

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

![](hypergeometric_MC_ntbk_files/figure-markdown_github/unnamed-chunk-9-1.png)

This next plotting function is to examine the distribution of sample successes for a specific number of successes in the population. Multiple population successes are allowed for comparison.

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
##  Plot entire sample distribution for successes in N.
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

We then call the function based on the parameters we want to see. Here we are looking at 4 successes in the sample size of 6 where we want to know the probability that the true number of successes in the population is greater than 25. For a sample size or 6 with 4 successes, the probability that the true number of successes in the population of 50 is greater than 25 is roughly 80%

``` r
##  parameters for plotting function
sample_successes <- 4 # number of successes in sample
pop_successes <- 25 # probability that true population successes are greater than this when observing sample_successes

right_tail_plot(agg_data_normed, sample_successes, pop_successes)
```

![](hypergeometric_MC_ntbk_files/figure-markdown_github/unnamed-chunk-14-1.png)
