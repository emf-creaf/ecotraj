# Distance metrics for trajectory resemblance

## 1. Introduction

### 1.1 About this vignette

In this vignette you will learn the differences between three distance
coefficients available for assessing trajectory resemblance. We use
small datasets where trajectories occur in a space of two dimensions, so
that geometric calculations can be followed more easily. First of all,
we load `ecotraj`:

``` r
library(ecotraj)
```

    ## Loading required package: Rcpp

### 1.2 The three distance metrics

Let $`T_1`$ and $`T_2`$ be two trajectories to be compared. The first
distance metric is **Segment Path Distance (SPD)** (Besse et al. 2016),
defined as the average of the distances between each *point* composing
$`T_1`$ and $`T_2`$:

``` math
D_{SP}(T_1, T_2) = \frac{1}{n}\sum_{i=1}^{n}{D_{pt}(x_{1i}, T_2)}
```

where $`D_{pt}`$ is the *distance between a point a a trajectory*. SPD
is not symmetric so it needs to be symmetrized:

``` math
D_{SSP}(T_1, T_2) = \frac{D_{SP}(T_1, T_2) + D_{SP}(T_2, T_1)}{2}
```

SPD is appropriate to compare the *location* and *shape* of
trajectories, but is not sensitive to trajectory *direction*. For this
reason, De Cáceres et al. (2019) introduced the **Directed Segment Path
Dissimilarity (DSPD)**, defined as the average of the distance between
each *directed segment* of $`T_1`$ and $`T_2`$:

``` math
D_{DSP}(T_1, T_2) = \frac{1}{n-1}\sum_{i=1}^{n-1}{D_{DS}(S_{1i}, T_2)}
```

where $`D_{SP}`$ is the *distance between a segment and a trajectory*.
As before, DSPD is not symmetric so it needs to be symmetrized:

``` math
D_{SDSP}(T_1, T_2) = \frac{D_{DSP}(T_1, T_2) + D_{DSP}(T_2, T_1)}{2}
```

DSPD is an appropriate metric to compare the *location*, *shape* and
*direction* of trajectories. Nevertheless, the metric does not allow
taking into account differences in trajectory *speed*, because it does
not use the information regarding the time of observations (only the
survey order).

If $`T_1`$ and $`T_2`$ represent the dynamics of two sites that have
been surveyed synchronously (i.e., if $`n = m`$ and $`t_{11} = t_{21}`$;
$`t_{12}=t_{22}`$;… ; $`t_{1n} = t_{2n}`$), a straightforward way of
comparing them is to calculate the average across surveys of
dissimilarity between the two sites, i.e. the mean of the sequence
$`\{d(x_{11}, x_{21}), \, d(x_{12}, x_{22}), \, \dots,\, d(x_{1n}, x_{2n})\}`$.
For a more general solution the **Time-Sensitive Path Distance (TSPD)**
is the average of distances between each *observation* in $`T_1`$ and
$`T_2`$:

``` math
D_{TSP}(T_1, T_2) = \frac{1}{n}\sum_{i=1}^{n}{D_{ot}(\{x_{1i}, t_{1i}\}, T_2)}
```

where $`D_{ot}`$ is the *distance between an observation and a
trajectory*. $`D_{ot}`$ is the calculated as the distance between
$`x_{1i}`$ and the point in $`T_2`$ corresponding to time $`t_{1i}`$,
which may need to be interpolated if does not correspond to any value in
$`\{t_{21}, t_{22}, \dots,t_{2m}\}`$. If $`t_{1i}`$ is beyond the time
boundaries of $`T_2`$, then the distance to the closest time point is
taken. As before, TSPD is not symmetric so it needs to be symmetrized:

``` math
D_{STSP}(T_1, T_2) = \frac{D_{TSP}(T_1, T_2) + D_{TSP}(T_2, T_1)}{2}
```

TSPD is sensitive to differences in *location*, *shape*, *direction* and
*speed*, as will be illustrated in the following examples.

## 2. Linear trajectories

Let us first compare the behavior of the three distance metrics for
comparisons between linear trajectories. In all cases, the reference
trajectory is composed of three linear segments.

### 2.1 Oposed linear trajectories

We compare first a linear trajectory with its oposed one, i.e. a
trajectory going in the exact oposite sense.

``` r
sites <- c("1","1","1","1","2","2","2", "2")
times <- c(0,1,2,3,0,1,2,3)
  
xy<-matrix(0, nrow=8, ncol=2)
xy[2,2]<-1
xy[3,2]<-2
xy[4,2]<-3
xy[5,2]<-3
xy[6,2]<-2
xy[7,2]<-1
xy[8,2]<-0

x <- defineTrajectories(dist(xy), sites = sites, times = times)
```

We can display the two (overlapping) trajectories using:

``` r
trajectoryPCoA(x, 
               traj.colors = c("black", "red", "blue"), lwd = 2,
               time.labels = TRUE)
```

![](TrajectoryDistanceMetrics_files/figure-html/unnamed-chunk-2-1.png)

The two trajectories have the same lengths and speeds:

``` r
trajectoryLengths(x)
```

    ##   S1 S2 S3 Path
    ## 1  1  1  1    3
    ## 2  1  1  1    3

``` r
trajectorySpeeds(x)
```

    ##   S1 S2 S3 Path
    ## 1  1  1  1    1
    ## 2  1  1  1    1

When we examine trajectory (temporal) shifts we see the oposing
character:

``` r
trajectoryShifts(x)
```

    ##   reference site survey time timeRef shift
    ## 1         1    2      1    0       3     3
    ## 2         1    2      2    1       2     1
    ## 3         1    2      3    2       1    -1
    ## 4         1    2      4    3       0    -3
    ## 5         2    1      1    0       3     3
    ## 6         2    1      2    1       2     1
    ## 7         2    1      3    2       1    -1
    ## 8         2    1      4    3       0    -3

Calculating SPD yields zero dissimilarity, because the distance does not
take into account differences in direction:

``` r
trajectoryDistances(x, distance.type = "SPD")
```

    ##   1
    ## 2 0

The other two dissimilarity metrics do yield non-zero values:

``` r
trajectoryDistances(x, distance.type = "DSPD")
```

    ##   1
    ## 2 1

``` r
trajectoryDistances(x, distance.type = "TSPD")
```

    ##   1
    ## 2 2

### 2.2 Equal pathways and speeds but different number of segments

Here we compare three trajectories with the same linear pathway and
speed. They only differ in the number of segments used to describe them:

``` r
sites <- c("1","1","1","1","2","2","2","3","3")
times <- c(0,1,2,3,0,1.5,3,0,3)
  
xy<-matrix(0, nrow=9, ncol=2)
xy[2,2]<-1
xy[3,2]<-2
xy[4,2]<-3
xy[6,2]<-1.5
xy[7,2]<-3
xy[9,2]<-3

x <- defineTrajectories(dist(xy), sites = sites, times = times)
```

We plot the three trajectories in separate panels for clarity:

``` r
par(mfrow=c(3,1))
trajectoryPCoA(subsetTrajectories(x,"1"), 
               traj.colors = c("black"), lwd = 2,
               time.labels = TRUE)
trajectoryPCoA(subsetTrajectories(x,"2"), 
               traj.colors = c("red"), lwd = 2,
               time.labels = TRUE)
trajectoryPCoA(subsetTrajectories(x,"3"),
               traj.colors = c("blue"), lwd = 2,
               time.labels = TRUE)
```

![](TrajectoryDistanceMetrics_files/figure-html/unnamed-chunk-8-1.png)
Note that reversals may occur because of PCoA eigen analysis. But
together the trajectories look like:

``` r
trajectoryPCoA(x, 
               traj.colors = c("black", "red", "blue"), lwd = 2,
               time.labels = TRUE)
```

![](TrajectoryDistanceMetrics_files/figure-html/unnamed-chunk-9-1.png)

We can check that the three trajectories have the same total length and
average speed:

``` r
trajectoryLengths(x)
```

    ##    S1  S2 S3 Path
    ## 1 1.0 1.0  1    3
    ## 2 1.5 1.5 NA    3
    ## 3 3.0  NA NA    3

``` r
trajectorySpeeds(x)
```

    ##   S1 S2 S3 Path
    ## 1  1  1  1    1
    ## 2  1  1 NA    1
    ## 3  1 NA NA    1

There are no temporal shifts between the trajectories:

``` r
trajectoryShifts(x)
```

    ##    reference site survey time timeRef shift
    ## 1          1    2      1  0.0     0.0     0
    ## 2          1    2      2  1.5     1.5     0
    ## 3          1    2      3  3.0     3.0     0
    ## 4          1    3      1  0.0     0.0     0
    ## 5          1    3      2  3.0     3.0     0
    ## 6          2    1      1  0.0     0.0     0
    ## 7          2    1      2  1.0     1.0     0
    ## 8          2    1      3  2.0     2.0     0
    ## 9          2    1      4  3.0     3.0     0
    ## 10         2    3      1  0.0     0.0     0
    ## 11         2    3      2  3.0     3.0     0
    ## 12         3    1      1  0.0      NA    NA
    ## 13         3    1      2  1.0      NA    NA
    ## 14         3    1      3  2.0      NA    NA
    ## 15         3    1      4  3.0      NA    NA
    ## 16         3    2      1  0.0      NA    NA
    ## 17         3    2      2  1.5     1.5     0
    ## 18         3    2      3  3.0      NA    NA

Here SPD yields zero distance, because the three trajectories have the
same shape:

``` r
trajectoryDistances(x, distance.type = "SPD")
```

    ##   1 2
    ## 2 0  
    ## 3 0 0

Since it is defined by means of distances between directed segments,
DSPD seems to be affected by the different segmentation of trajectories,
so that it yields non-zero values:

``` r
trajectoryDistances(x, distance.type = "DSPD")
```

    ##           1         2
    ## 2 0.5833333          
    ## 3 1.3333333 1.5000000

In contrast, TSPD yields again zero distance values, because the
trajectories do not differ in neither speed or shape.

``` r
trajectoryDistances(x, distance.type = "TSPD")
```

    ##   1 2
    ## 2 0  
    ## 3 0 0

### 2.3 Equal pathways but different speeds

In this example the three trajectories have the same segments and
pathways, but they differ in the speed of changes:

``` r
sites <- c("1","1","1","1","2","2","2","2","3","3","3","3")
times <- c(0,0.5,1,1.5,0,1,2,3,0,2,4,6)
  
xy<-matrix(0, nrow=12, ncol=2)
xy[2,2]<-1
xy[3,2]<-2
xy[4,2]<-3
xy[5:8,2] <- xy[1:4,2]
xy[9:12,2] <- xy[1:4,2]

x <- defineTrajectories(dist(xy), sites = sites, times = times)
```

Again we use separate plots to show the differences in speed:

``` r
par(mfrow=c(3,1))
trajectoryPCoA(subsetTrajectories(x,"1"), 
               traj.colors = c("black"), lwd = 2,
               time.labels = TRUE)
trajectoryPCoA(subsetTrajectories(x,"2"), 
               traj.colors = c("red"), lwd = 2,
               time.labels = TRUE)
trajectoryPCoA(subsetTrajectories(x,"3"), 
               traj.colors = c("blue"), lwd = 2,
               time.labels = TRUE)
```

![](TrajectoryDistanceMetrics_files/figure-html/unnamed-chunk-17-1.png)
We can check that there are no differences in segment or total path
lengths, but they indeed differ in trajectory speed:

``` r
trajectoryLengths(x)
```

    ##   S1 S2 S3 Path
    ## 1  1  1  1    3
    ## 2  1  1  1    3
    ## 3  1  1  1    3

``` r
trajectorySpeeds(x)
```

    ##    S1  S2  S3 Path
    ## 1 2.0 2.0 2.0  2.0
    ## 2 1.0 1.0 1.0  1.0
    ## 3 0.5 0.5 0.5  0.5

Differences in speed also lead to temporal shifts between trajectories:

``` r
trajectoryShifts(x)
```

    ##    reference site survey time timeRef shift
    ## 1          1    2      1  0.0     0.0   0.0
    ## 2          1    2      2  1.0     0.5  -0.5
    ## 3          1    2      3  2.0     1.0  -1.0
    ## 4          1    2      4  3.0     1.5  -1.5
    ## 5          1    3      1  0.0     0.0   0.0
    ## 6          1    3      2  2.0     0.5  -1.5
    ## 7          1    3      3  4.0     1.0  -3.0
    ## 8          1    3      4  6.0     1.5  -4.5
    ## 9          2    1      1  0.0     0.0   0.0
    ## 10         2    1      2  0.5     1.0   0.5
    ## 11         2    1      3  1.0     2.0   1.0
    ## 12         2    1      4  1.5     3.0   1.5
    ## 13         2    3      1  0.0     0.0   0.0
    ## 14         2    3      2  2.0     1.0  -1.0
    ## 15         2    3      3  4.0     2.0  -2.0
    ## 16         2    3      4  6.0     3.0  -3.0
    ## 17         3    1      1  0.0     0.0   0.0
    ## 18         3    1      2  0.5     2.0   1.5
    ## 19         3    1      3  1.0     4.0   3.0
    ## 20         3    1      4  1.5     6.0   4.5
    ## 21         3    2      1  0.0     0.0   0.0
    ## 22         3    2      2  1.0     2.0   1.0
    ## 23         3    2      3  2.0     4.0   2.0
    ## 24         3    2      4  3.0     6.0   3.0

If we calculate distances using SPD, the distance metric does not detect
the differences in speed and tells us that the trajectories are equal:

``` r
trajectoryDistances(x, distance.type = "SPD")
```

    ##   1 2
    ## 2 0  
    ## 3 0 0

And the same happens with DSPD:

``` r
trajectoryDistances(x, distance.type = "DSPD")
```

    ##   1 2
    ## 2 0  
    ## 3 0 0

It is only when we apply TSPD that we can observe differences between
trajectories:

``` r
trajectoryDistances(x, distance.type = "TSPD")
```

    ##        1      2
    ## 2 0.6250       
    ## 3 0.9375 0.6250

The distance between the first and the third trajectory is largest
because their difference in speed is also largest.

### 2.4 Space-shifted trajectories

Let us now evaluate a case where trajectories are the same but have been
displaced in one dimension:

``` r
sites <- c("1","1","1","1","2","2","2","2","3","3","3","3")
times <- c(1,2,3,4,1,2,3,4,1,2,3,4)

xy<-matrix(0, nrow=12, ncol=2)
xy[2,2]<-1
xy[3,2]<-2
xy[4,2]<-3
# States are all shifted half unit with respect to site "1"
xy[5:8,2] <- xy[1:4,2] + 0.5 
# States are all shifted one unit with respect to site "1"
xy[9:12,2] <- xy[1:4,2] + 1.0  

x <- defineTrajectories(dist(xy), sites = sites, times = times)
```

We use a single plot, though not very clear, to display the three
trajectories:

``` r
trajectoryPCoA(x, 
               traj.colors = c("black", "red", "blue"), lwd = 2,
               time.labels = TRUE)
```

![](TrajectoryDistanceMetrics_files/figure-html/unnamed-chunk-24-1.png)
In this case differences do not exist in terms of lengths nor speeds:

``` r
trajectoryLengths(x)
```

    ##   S1 S2 S3 Path
    ## 1  1  1  1    3
    ## 2  1  1  1    3
    ## 3  1  1  1    3

``` r
trajectorySpeeds(x)
```

    ##   S1 S2 S3 Path
    ## 1  1  1  1    1
    ## 2  1  1  1    1
    ## 3  1  1  1    1

But (temporal) shifts reflect the spatial ones:

``` r
trajectoryShifts(x)
```

    ##    reference site survey time timeRef shift
    ## 1          1    2      1    1     1.5   0.5
    ## 2          1    2      2    2     2.5   0.5
    ## 3          1    2      3    3     3.5   0.5
    ## 4          1    2      4    4      NA    NA
    ## 5          1    3      1    1     2.0   1.0
    ## 6          1    3      2    2     3.0   1.0
    ## 7          1    3      3    3     4.0   1.0
    ## 8          1    3      4    4      NA    NA
    ## 9          2    1      1    1      NA    NA
    ## 10         2    1      2    2     1.5  -0.5
    ## 11         2    1      3    3     2.5  -0.5
    ## 12         2    1      4    4     3.5  -0.5
    ## 13         2    3      1    1     1.5   0.5
    ## 14         2    3      2    2     2.5   0.5
    ## 15         2    3      3    3     3.5   0.5
    ## 16         2    3      4    4      NA    NA
    ## 17         3    1      1    1      NA    NA
    ## 18         3    1      2    2     1.0  -1.0
    ## 19         3    1      3    3     2.0  -1.0
    ## 20         3    1      4    4     3.0  -1.0
    ## 21         3    2      1    1      NA    NA
    ## 22         3    2      2    2     1.5  -0.5
    ## 23         3    2      3    3     2.5  -0.5
    ## 24         3    2      4    4     3.5  -0.5

In this case, all three metrics are responsive to differences in
trajectory location:

``` r
trajectoryDistances(x, distance.type = "SPD")
```

    ##       1     2
    ## 2 0.125      
    ## 3 0.250 0.125

``` r
trajectoryDistances(x, distance.type = "DSPD")
```

    ##           1         2
    ## 2 0.5000000          
    ## 3 0.3333333 0.5000000

``` r
trajectoryDistances(x, distance.type = "TSPD")
```

    ##     1   2
    ## 2 0.5    
    ## 3 1.0 0.5

### 2.5 Space-expanded trajectories

In this example, the three linear trajectories are surveyed the same
times but they differ in total path length due to differences in
trajectory speed.

``` r
sites <- c("1","1","1","1","2","2","2","2","3","3","3","3")
times <- c(0,1,2,3,0,1,2,3,0,1,2,3)
  
xy<-matrix(0, nrow=12, ncol=2)
xy[2,2]<-1
xy[3,2]<-2
xy[4,2]<-3
xy[5:8,2] <- xy[1:4,2]*1.5
xy[9:12,2] <- xy[1:4,2]*2

x <- defineTrajectories(dist(xy), sites = sites, times = times)
```

We draw the three (overlapping) trajectories:

``` r
trajectoryPCoA(x, 
               traj.colors = c("black", "red", "blue"), lwd = 2,
               time.labels = TRUE)
```

![](TrajectoryDistanceMetrics_files/figure-html/unnamed-chunk-29-1.png)

In this case both lengths and speeds are different between trajectories:

``` r
trajectoryLengths(x)
```

    ##    S1  S2  S3 Path
    ## 1 1.0 1.0 1.0  3.0
    ## 2 1.5 1.5 1.5  4.5
    ## 3 2.0 2.0 2.0  6.0

``` r
trajectorySpeeds(x)
```

    ##    S1  S2  S3 Path
    ## 1 1.0 1.0 1.0  1.0
    ## 2 1.5 1.5 1.5  1.5
    ## 3 2.0 2.0 2.0  2.0

This is also translated to trajectory shifts:

``` r
trajectoryShifts(x)
```

    ##    reference site survey time   timeRef      shift
    ## 1          1    2      1    0 0.0000000  0.0000000
    ## 2          1    2      2    1 1.5000000  0.5000000
    ## 3          1    2      3    2 3.0000000  1.0000000
    ## 4          1    2      4    3        NA         NA
    ## 5          1    3      1    0 0.0000000  0.0000000
    ## 6          1    3      2    1 2.0000000  1.0000000
    ## 7          1    3      3    2        NA         NA
    ## 8          1    3      4    3        NA         NA
    ## 9          2    1      1    0 0.0000000  0.0000000
    ## 10         2    1      2    1 0.6666667 -0.3333333
    ## 11         2    1      3    2 1.3333333 -0.6666667
    ## 12         2    1      4    3 2.0000000 -1.0000000
    ## 13         2    3      1    0 0.0000000  0.0000000
    ## 14         2    3      2    1 1.3333333  0.3333333
    ## 15         2    3      3    2 2.6666667  0.6666667
    ## 16         2    3      4    3        NA         NA
    ## 17         3    1      1    0 0.0000000  0.0000000
    ## 18         3    1      2    1 0.5000000 -0.5000000
    ## 19         3    1      3    2 1.0000000 -1.0000000
    ## 20         3    1      4    3 1.5000000 -1.5000000
    ## 21         3    2      1    0 0.0000000  0.0000000
    ## 22         3    2      2    1 0.7500000 -0.2500000
    ## 23         3    2      3    2 1.5000000 -0.5000000
    ## 24         3    2      4    3 2.2500000 -0.7500000

Since the trajectories differ in length, this is captured by all three
metrics:

``` r
trajectoryDistances(x, distance.type = "SPD")
```

    ##        1      2
    ## 2 0.1875       
    ## 3 0.5000 0.1875

``` r
trajectoryDistances(x, distance.type = "DSPD")
```

    ##           1         2
    ## 2 0.7500000          
    ## 3 1.3333333 0.9166667

``` r
trajectoryDistances(x, distance.type = "TSPD")
```

    ##      1    2
    ## 2 0.75     
    ## 3 1.50 0.75

## 3. Curved trajectories

In this second set of examples we examine the behavior of the metrics
when comparing trajectories that are not always linear.

### 3.1 Constant speed

Here the three trajectories have the same length and speed, but
trajectory 2 and 3 are progressively curved:

``` r
sites <- c("1","1","1","1","2","2","2","2","3","3","3","3")
surveys <- c(1,2,3,4,1,2,3,4,1,2,3,4)

xy<-matrix(0, nrow=12, ncol=2)
xy[2,2]<-1
xy[3,2]<-2
xy[4,2]<-3
xy[5:6,2] <- xy[1:2,2]
xy[7,1]<-0+sqrt(0.5)
xy[7,2]<-1+sqrt(0.5)
xy[8,2]<-xy[7,2]
xy[8,1]<-xy[7,1]+1
xy[9:10,2] <- xy[1:2,2]
xy[11,1] <- 1.0
xy[11,2] <- 1.0
xy[12,1] <- 1.0
xy[12,2] <- 0.0

x <- defineTrajectories(dist(xy), sites = sites, times = times)
```

``` r
trajectoryPCoA(x, 
               traj.colors = c("black", "red", "blue"), lwd = 2,
               time.labels = TRUE)
```

    ## Warning in cmdscale(d, eig = TRUE, add = TRUE, k = nrow(as.matrix(d)) - : only
    ## 10 of the first 11 eigenvalues are > 0

![](TrajectoryDistanceMetrics_files/figure-html/unnamed-chunk-34-1.png)

As expected, no differences are found in terms of lengths or speeds:

``` r
trajectoryLengths(x)
```

    ##   S1 S2 S3 Path
    ## 1  1  1  1    3
    ## 2  1  1  1    3
    ## 3  1  1  1    3

``` r
trajectorySpeeds(x)
```

    ##   S1 S2 S3 Path
    ## 1  1  1  1    1
    ## 2  1  1  1    1
    ## 3  1  1  1    1

The three distance metrics are responsive to differences of trajectory
shape:

``` r
trajectoryDistances(x, distance.type = "SPD")
```

    ##           1         2
    ## 2 0.5743683          
    ## 3 0.6250000 0.4267767

``` r
trajectoryDistances(x, distance.type = "DSPD")
```

    ##           1         2
    ## 2 0.7658244          
    ## 3 0.9023689 0.6868867

``` r
trajectoryDistances(x, distance.type = "TSPD")
```

    ##           1         2
    ## 2 0.7267030          
    ## 3 1.1441228 0.6532815

### 3.2 Different speed

This example is similar to the previous one, but here we changed the
survey times, so that the observed trajectory shapes correspond also to
different speeds:

``` r
sites <- c("1","1","1","1","2","2","2","2","3","3","3","3")
times <- c(0,0.5,1,1.5,0,1,2,3,0,2,4,6)

xy<-matrix(0, nrow=12, ncol=2)
xy[2,2]<-1
xy[3,2]<-2
xy[4,2]<-3
xy[5:6,2] <- xy[1:2,2]
xy[7,1]<-0+sqrt(0.5)
xy[7,2]<-1+sqrt(0.5)
xy[8,2]<-xy[7,2]
xy[8,1]<-xy[7,1]+1
xy[9:10,2] <- xy[1:2,2]
xy[11,1] <- 1.0
xy[11,2] <- 1.0
xy[12,1] <- 1.0
xy[12,2] <- 0.0

x <- defineTrajectories(dist(xy), sites = sites, times = times)
```

The trajectory plot looks as before, except for the time labels:

``` r
trajectoryPCoA(x, 
               traj.colors = c("black", "red", "blue"), lwd = 2,
               time.labels = TRUE)
```

    ## Warning in cmdscale(d, eig = TRUE, add = TRUE, k = nrow(as.matrix(d)) - : only
    ## 10 of the first 11 eigenvalues are > 0

![](TrajectoryDistanceMetrics_files/figure-html/unnamed-chunk-38-1.png)

In this case trajectories differ in speed but not length:

``` r
trajectoryLengths(x)
```

    ##   S1 S2 S3 Path
    ## 1  1  1  1    3
    ## 2  1  1  1    3
    ## 3  1  1  1    3

``` r
trajectorySpeeds(x)
```

    ##    S1  S2  S3 Path
    ## 1 2.0 2.0 2.0  2.0
    ## 2 1.0 1.0 1.0  1.0
    ## 3 0.5 0.5 0.5  0.5

The three metrics detect differences in shape, as before.

``` r
trajectoryDistances(x, distance.type = "SPD")
```

    ##           1         2
    ## 2 0.5743683          
    ## 3 0.6250000 0.4267767

``` r
trajectoryDistances(x, distance.type = "DSPD")
```

    ##           1         2
    ## 2 0.7658244          
    ## 3 0.9023689 0.6868867

``` r
trajectoryDistances(x, distance.type = "TSPD")
```

    ##           1         2
    ## 2 0.9748813          
    ## 3 1.4872932 0.8433407

However, note that the values of SPD and DSPD are exactly the same to
those of the previous example, whereas TSPD yields higher distance
values because of the differences in trajectory speed.

## 4. References

- Besse, P., Guillouet, B., Loubes, J.-M. & François, R. (2016). Review
  and perspective for distance based trajectory clustering. IEEE Trans.
  Intell. Transp. Syst., 17, 3306–3317.

- De Cáceres M, Coll L, Legendre P, Allen RB, Wiser SK, Fortin MJ,
  Condit R & Hubbell S. (2019). Trajectory analysis in community
  ecology. Ecological Monographs 89, e01350.
