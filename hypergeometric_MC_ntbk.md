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
    ## 0 10000 8804 7724 6777 5908 5138 4477 3868 3254 2892 2454 2012 1781 1447
    ## 1     0 1196 2164 2877 3416 3821 4053 4161 4250 4113 4168 3951 3817 3556
    ## 2     0    0  112  333  634  944 1313 1670 2036 2341 2540 2944 3044 3232
    ## 3     0    0    0   13   42   94  151  278  411  584  736  946 1131 1431
    ## 4     0    0    0    0    0    3    6   23   47   67   95  143  219  307
    ## 5     0    0    0    0    0    0    0    0    2    3    7    4    8   26
    ## 6     0    0    0    0    0    0    0    0    0    0    0    0    0    1
    ##     14   15   16   17   18   19   20   21   22   23   24   25   26   27
    ## 0 1214 1051  862  668  549  445  391  287  215  173  150  111   96   56
    ## 1 3365 3132 2821 2568 2364 2033 1826 1486 1414 1114  992  851  737  568
    ## 2 3325 3324 3452 3446 3452 3485 3224 3207 2939 2845 2614 2306 2114 2006
    ## 3 1662 1905 2121 2388 2475 2698 2912 3068 3167 3231 3381 3415 3345 3230
    ## 4  372  507  661  798  969 1106 1329 1555 1716 1990 2111 2370 2584 2864
    ## 5   60   75   79  123  177  215  298  369  504  566  664  819  982 1098
    ## 6    2    6    4    9   14   18   20   28   45   81   88  128  142  178
    ##     28   29   30   31   32   33   34   35   36   37   38   39   40   41
    ## 0   44   37   28   21   13    8    4    2    1    2    0    0    0    0
    ## 1  431  356  300  242  171  123  114   69   38   28   19   12   10    1
    ## 2 1708 1488 1332 1135  955  781  675  505  399  309  218  158   95   49
    ## 3 3226 3095 2924 2810 2547 2297 2070 1838 1663 1425 1164  958  745  584
    ## 4 3033 3131 3262 3301 3465 3490 3517 3478 3411 3178 3111 2871 2608 2292
    ## 5 1330 1597 1758 2021 2282 2544 2752 3072 3264 3597 3754 3921 4150 4284
    ## 6  228  296  396  470  567  757  868 1036 1224 1461 1734 2080 2392 2790
    ##     42   43   44   45   46   47   48   49    50
    ## 0    0    0    0    0    0    0    0    0     0
    ## 1    2    1    0    0    0    0    0    0     0
    ## 2   47   22   10    5    1    0    0    0     0
    ## 3  447  288  177  102   34   13    0    0     0
    ## 4 1999 1625 1303  908  624  310  123    0     0
    ## 5 4254 4175 4109 3813 3454 2861 2159 1162     0
    ## 6 3251 3889 4401 5172 5887 6816 7718 8838 10000

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

![](hypergeometric_MC_ntbk_files/figure-markdown_github/unnamed-chunk-12-1.png)

If we want to examine a specific number of successes in a sample and find the likelihood of the true number of population successes we can plot the sample and sum the area under the curve. This function does that for the right tail of the distribution.

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
