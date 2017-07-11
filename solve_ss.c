/*==========================================================
 * spread.c - example in MATLAB External Interfaces
 *========================================================*/

#include "mex.h"
#include <math.h>

void print_array(mwSize n, double *array) {
    mwSize i;
    for (i = 0; i < n; i++) {
        printf("%f ",array[i]);
    }
    printf("\n");
}

double interp(mwSize n, double *x_grid, double *y_grid, double xx) {
    mwSize i;
    double slope;

    if (xx <= x_grid[0]) {
        return y_grid[0];
    } else {
        for (i = 1; i < n; i++) {
            if (xx < x_grid[i]) {
                slope = (y_grid[i] - y_grid[i-1]) / (x_grid[i] - x_grid[i-1]);
                return y_grid[i-1] + (xx - x_grid[i-1]) * slope;
            }
        }
        return y_grid[n-1];
    }
}

double foc(double kp, double *params, mwSize nz, double *z0, double *Pz_row,
    mwSize nk, double *k0, double k, double *kp_grid)
{
    double alpha, beta, delta, phi, one_minus_delta;
    double lhs, rhs;
    double kpp, kpp_over_kp, prodn_part, adj_part;
    mwSize zp_idx;

    alpha = params[0];
    beta = params[1];
    delta = params[2];
    phi = params[3];

    one_minus_delta = 1.0 - delta;

    lhs = 1 + 2 * phi * (kp/k-(1-delta));

    rhs = 0.0;
    for (zp_idx = 0; zp_idx < nz; zp_idx++) {
        kpp = interp(nk, k0, kp_grid + zp_idx*nk, kp);
        kpp_over_kp = kpp/kp;
        adj_part = phi * (kpp_over_kp*kpp_over_kp - one_minus_delta*one_minus_delta);
        prodn_part = z0[zp_idx] * alpha * pow(kp,alpha-1.0);

        rhs += beta * Pz_row[zp_idx] * (prodn_part + one_minus_delta + adj_part);
    }
    return lhs-rhs;
}

double abs_foc(double kp, double *params, mwSize nz, double *z0, double *Pz_row,
    mwSize nk, double *k0, double k, double *kp_grid) {
    return fabs(foc(kp, params, nz, z0, Pz_row, nk, k0, k, kp_grid));
}

double golden(double x_min, double x_max, double *params, double *options, mwSize nz,
    double *z0, double *Pz_row, mwSize nk, double *k0, double k, double *kp_grid)
{
    mwSize i, max_iter;
    double a, b, c, d, gr, solver_tolerance, fc, fd, err;
    double x_star;

    max_iter = 100;
    solver_tolerance = options[1];
    gr = 1.618;

    a = x_min;
    b = x_max;

    c = b - (b - a) / gr;
    d = a + (b - a) / gr;

    for (i = 0; i < max_iter; i++) {
        fc = abs_foc(c, params, nz, z0, Pz_row, nk, k0, k, kp_grid);
        fd = abs_foc(d, params, nz, z0, Pz_row, nk, k0, k, kp_grid);
        
        if (fc < fd) {
            b = d;
        } else {
            a = c;
        }
        c = b - (b - a) / gr;
        d = a + (b - a) / gr;

        err = fabs(c - d);

        if (err < solver_tolerance) {
            break;
        }
    }
    x_star = (b + a) * 0.5;
    return x_star;
}

void solve_ss(double *params, double *options, mwSize nk, double *k0, mwSize nz, double *z0,
    double *Pzp, double *kp_guess, double *output_ptr)
{
    mwSize k_idx,z_idx,zp_idx,idx;
    mwSize iter;
    mwSize max_iter = 100;
    double abs_diff, max_abs_diff;
    double policy_tolerance;
    double kp[nk*nz], kp_old[nk*nz];
    double kp_star;
    double foc_lb, foc_ub;

    policy_tolerance = options[0];

    /* copy kp_guess into kp */
    for (idx = 0; idx < (nk*nz); idx++) {
        kp[idx] = kp_guess[idx];
    }

    /* iteration */
    for (iter = 0; iter < max_iter; iter++) {
        /* copy kp into kp_old */
        for (idx = 0; idx < (nk*nz); idx++) {
            kp_old[idx] = kp[idx];
        }

        /* loop over all points in state-space */
        max_abs_diff = 0;
        idx = 0;
        for (z_idx = 0; z_idx < nz; z_idx++) {
            for (k_idx = 0; k_idx < nk; k_idx++) {

                /* check lower bound */
                foc_lb = foc(k0[0], params, nz, z0, Pzp+z_idx*nz, nk, k0, k0[k_idx], kp_old);

                if (foc_lb > 0) {
                    kp_star = k0[0];
                } else {
                    /* check upper bound */
                    foc_ub = foc(k0[nk-1], params, nz, z0, Pzp+z_idx*nz, nk, k0, k0[k_idx], kp_old);

                    if (foc_ub < 0) {
                        kp_star = k0[nk-1];
                    } else {
                        /* find interior solution */
                        kp_star = golden(k0[0], k0[nk-1], params, options, nz, z0, Pzp+z_idx*nz, nk, k0, k0[k_idx], kp_old);
                    }
                }

                kp[idx] = kp_star;

                abs_diff = fabs(kp_star - kp_old[idx]);
                if (abs_diff > max_abs_diff) {
                    max_abs_diff = abs_diff;
                }

                idx++;
            }
        }
        if (max_abs_diff < policy_tolerance) {
            break;
        }
    }
    /* printf("%d\n",iter); */

    for (idx = 0; idx < (nk*nz); idx++) {
        output_ptr[idx] = kp[idx];
    }
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    mwSize nk, nz;
    double *params, *k0, *z0, *Pzp, *kp_guess, *options;
    double *output_ptr;

    /* check for proper number of arguments */
    if(nrhs!=6) {
        mexErrMsgIdAndTxt("MyToolbox:spread:nrhs","5 inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:spread:nlhs","1 output required.");
    }

    /* create a pointer to the real data in the input matrix  */
    params = mxGetPr(prhs[0]);
    k0 = mxGetPr(prhs[1]);
    z0 = mxGetPr(prhs[2]);
    Pzp = mxGetPr(prhs[3]);
    kp_guess = mxGetPr(prhs[4]);
    options = mxGetPr(prhs[5]);

    nk = mxGetM(prhs[1]);
    nz = mxGetM(prhs[2]);

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix((mwSize)nk,(mwSize)nz,mxREAL);
    output_ptr = mxGetPr(plhs[0]);

    /* call the computational routine */
    solve_ss(params, options, nk, k0, nz, z0, Pzp, kp_guess, output_ptr);
}
