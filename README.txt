Bayesian Multiply Robust Estimaton of the Average Causal Effect

Please follow the following steps to compute the multiply robust projection and calibration estimator (Gochanour et al., 2021+) of
the Average Causal Effect (ACE). Hopefully this code will soon be implemented in an R package to make things even easier. Until then,
please follow the following steps.

(1) Download files

Download R Code.zip to your computer. Unzip the folder and place all files inside your preferred R working directory.

(2) Prepare your datafile

Before applying our method, coerce your data into a dataframe with the following conventions:

Dataframe name should be "dat"

Covariates, which may be continuous or categorical, should be named/renamed x1, x2, ...., xn, where n is
the total number of covariates desired.

Treatment variable should be numeric with the values 0 and 1 and named "r"

Your outcome variables for treatments 1 and 0 should be named "y1" and "y0" and should be numeric.

IMPORTANT: COLUMNS MUST BE IN THE FOLLOWING ORDER:
"x1", "x2", ...., "xn", "r", "y1", "y0"

Counsult "dat.RData" for an example. This data is the NHANES data discussed in Gochanour et al. (2021+), 
processed to follow the above conventions.

(3) Specify your models

Open the Main.R script. Set working directory to source file location.

Modify lines 14-26 of Main.R to reflect your desired candidate models (two regression and two propensity).

(4) Run the script

Run the entire "Main.R" script from top to bottom. If you like, you can set a different C value in the function calls. Here,
C is the desired number of sets of weights (called "T" in the paper).

(5) Examine your results

Your results will be returned to the object "finaltable" when the code finishes running. The results present, for each estimator, the point estimate
psi and its 95% Bayesian credible interval.

NOTE: Results of running the code on the example data may not exactly match those of Gochanour et al. (2021+). This is because generating the weights is
a random process. If you'd like results to be consistent from run to run, please specify a seed.
