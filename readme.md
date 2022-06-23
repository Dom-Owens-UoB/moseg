# moseg

Methods for data segmentation under a sparse regression model. See

Reference


## Installation

package `moseg` installable via

```
devtools::install_github("https://github.com/Dom-Owens-UoB/moseg")
```
 

## Usage

We can simulate from a piecewise sparse regression model via
```
set.seed(111)
dat <- moseg.sim(500, 50, q = 2, kappa = 4)
```

Identify change points:
```
out <- moseg(dat$X, dat$y, 100, do.scale = FALSE)
```

Multiscale:

```
out <- moseg.ms(dat$X, dat$y, c(50,100,150), do.scale = FALSE)
```

Using cross-validation:
```
cv <- moseg.cv(dat$X, dat$y, 50, do.scale = FALSE)
cv.ms <- moseg.ms.cv(dat$X, dat$y, c(50,100,150), do.scale = FALSE)
```
