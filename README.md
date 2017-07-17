# hadsge
Code for Bayesian estimation of a heterogeneous agent DSGE model (MATLAB).

## Running the code

### Dependencies

1. **Gensys**, available from http://sims.princeton.edu/yftp/gensys/
  
2. **myAD** automatic differentiation package, available from https://github.com/sehyoun/MATLABAutoDiff

Optional:
The MATLAB file ``functions/solve_ss.m`` solves for the model's steady state. I also provide a C implementation of this function which can be compiled as a MEX file. To compile, in MATLAB run
```
mex functions/solve_ss.c
```
in the _hadsge_ directory.

### Running
  
The file ``run.m``  has several sections:

#### Section 1: Solving

The model solution method is given in the ``@InvestmentModel`` class. The code
```
m = InvestmentModel(parameters);
m.solve();
```
computes the model solution where ``parameters`` is a vector of model parameters (see ``@InvestmentModel/InvestmentModel.m`` for a description).

Various aspects of the model solution (both steady-state and dynamic) may be accessed as properties of the InvestmentModel object after the ``solve()`` method has been called. See the source file for details.

The model solution can be split into three parts:

1. Solve for steady state:
``
m.solve_steady_state();
``

2. Compute derivatives:
``
m.compute_derivatives();
``

3. Solve for dynamics: ``m.solve_dynamics();``

This is to facilitate re-solving after changing parameters that don't affect the steady state or don't affect either the steady state or the derivatives.

####  Section 2: Simulating

The  ``@InvestmentModel`` class also provides a method to simulate from the solved model. Gaussian innovations must be provided by the user:

```
T = 200;
innov = randn(3,T);
data = m.simulate(innov);
```

(Note: the measurement equation may be amendede by changing ``m.C`` and ``m.D``.)

### Section 3: Evaluating the posterior

Given some data (e.g. the simulated data ``data`` above), the posterior for a given set of parameters may be calculated. Gaussian innovations must be provided by the user. These are used to compute the covariance matrix of the state equation via simulation.

```
T_for_sim = 1000;
innov_for_sim = randn(3, T_for_sim);

posterior_obj_fn = @(theta) m.evaluate_logposterior_for_parameters(theta, data, innov_for_sim);

log_posterior = posterior_obj_fn(estimation_parameters)
```

where ``estimation_parameters`` is a vector of parameters. Note that the ``estimation_parameters`` vector is shorter than the initializing ``parameters`` vector above. The initializing ``parameters`` vector includes some parameters which are assumed fixed (i.e. $\alpha$, $\beta$, $\delta$). The division of parameters into fixed and for-estimation is easy to modify in the source.

### Section 4: Estimating using Metropolis-Hastings

Given an objective function, the ``metropolis_hastings()`` function executes the Metropolis-Hastings algorithm and writes the output chains (in CSV format without headers) to a file passed as an input argument (see the source for a full description of input arguments).

This section runs parallel chains using MATLAB's parallel toolbox if ``use_parallel_toolbox`` is set to ``true`` or a single chain in serial if not.

### Section 5: Loading and plotting estimation results

The ``@MCMCResults`` class loads the output of a successful run of the ``metropolis_hastings()`` function. The user must pass as an input argument the path to the output files on disk. Various elements of the output may be accessed as properties of the ``@MCMCResults`` object:

```
fraction_to_keep = 0.75;
results = MCMCResults('mh_output', fraction_to_keep)

parameter_means = mean(results.combined(:,1:num_parameters))
```
where ``num_parameters`` is the number of estimated parameters.


## Model

...
