# mosumsr

Methods for data segmentation under a sparse regression model. See

Reference


## Installation

package `mosumsr` installable via

```
devtools::install_github("https://github.com/Dom-Owens-UoB/sparse_mosum", subdir = "mosumsr")
```
 

## Usage

We can simulate from a piecewise sparse regression model via
```
set.seed(111)
dat <- mosumsr.sim(500, 50, q = 2, kappa = 4)
```

Identify change points:
```
out <- mosumsr(dat$X, dat$y, 100, do.scale = FALSE)
```

Multiscale:

```
out <- mosumsr.multiscale(dat$X, dat$y, c(50,100,150), do.scale = FALSE)
```

Using cross-validation:
```
cv <- mosumsr.cv(dat$X, dat$y, 50, do.scale = FALSE)
cv.ms <- mosumsr.multiscale.cv(dat$X, dat$y, c(50,100,150), do.scale = FALSE)
```