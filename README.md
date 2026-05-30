
# BayesRMST

`BayesRMST` is an R package for Bayesian analysis of Restricted Mean Survival Time (RMST) models using Stan.

## Prerequisites

Before installing `BayesRMST`, users must ensure that:

1. A working C++ toolchain/compiler is installed.
2. The `rstan` package is successfully installed and configured.

Since `BayesRMST` relies on Stan for Bayesian computation, proper compiler and `rstan` installation are required.

---

## Step 1: Verify Compiler Installation

### Windows

Install **Rtools** corresponding to your R version:

https://cran.r-project.org/bin/windows/Rtools/

After installation, verify that the compiler is available:

```r
pkgbuild::has_build_tools(debug = TRUE)
```

If the result is `TRUE`, your compiler setup is ready.

### macOS

Install Apple's Command Line Tools:

```bash
xcode-select --install
```

### Linux

Install the GNU build tools. For Ubuntu/Debian:

```bash
sudo apt update
sudo apt install build-essential
```

---

## Step 2: Install RStan

Install the required packages:

```r
install.packages(c("StanHeaders", "rstan"), dependencies = TRUE)
```

Load `rstan` and verify the installation:

```r
library(rstan)

example(stan_model, package = "rstan", run.dontrun = TRUE)
```

You can also check the Stan configuration:

```r
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
```

For detailed platform-specific instructions, consult the RStan installation guide:

https://mc-stan.org/users/interfaces/rstan

---

## Step 3: Install devtools

If you do not already have `devtools` installed:

```r
install.packages("devtools")
```

---

## Step 4: Install BayesRMST from GitHub

Install the development version directly from GitHub:

```r
devtools::install_github("ksbakar/BayesRMST")
```

---

## Step 5: Load the Package

```r
library(BayesRMST)
```

Verify the installation:

```r
packageVersion("BayesRMST")
```

---

## Troubleshooting

### RStan Compilation Errors

If Stan models fail to compile:

1. Confirm that your compiler toolchain is installed correctly.
2. Restart R after installing Rtools (Windows) or Command Line Tools (macOS).
3. Update `rstan` and `StanHeaders`:

```r
install.packages(c("StanHeaders", "rstan"))
```

4. Check your compiler configuration:

```r
pkgbuild::has_build_tools(debug = TRUE)
```

### Windows-Specific Issues

Ensure that:

* Rtools is installed.
* Rtools is on the system PATH.
* Your R version matches the installed Rtools version.

---

## Example

```r
library(BayesRMST)

## Example analysis

## Datasets for evaluation of cetuximab in head and neck cancer

data(cetux, package="survextrap")
head(cetux)
y <- cetux[,c("years","d","treat")] # time, event, trt
head(y)

## North Central Cancer Treatment Group: Lung Cancer Data

data(cancer, package="survival")
x <- cancer[,c("time","status","sex")] # time, event, trt
x$sex = x$sex-1; x$status = x$status-1
head(x)

## Acute Myelogenous Leukemia survival data

head(aml) # time, event, trt
aml$x = as.numeric(aml$x)-1
head(aml)

## Bayesian M-spline

out_ms <- survRMST(data=y, modelType="BayesMspline")
print(out_ms)
plot(out_ms)
plot(out_ms, omega=TRUE)

## Bayesian non-parametric (DDP)

out_dp <- survRMST(data=y, modelType="BayesNonPara")
print(out_dp)
plot(out_dp)

## Bayesian parametric

out_para <- survRMST(data=y, modelType="BayesPara")
print(out_para)
plot(out_para)


```

## Reporting Issues

Please report bugs, installation issues, or feature requests through the GitHub Issues page for the repository.

## License

See the `LICENSE` file for licensing information.
