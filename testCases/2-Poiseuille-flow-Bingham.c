/**
#  Planar Couette flow of Generalized Newtonian Fluid

This code extends the method used in [/sandbox/M1EMN/Exemples/bingham_simple.c](/../sandbox/M1EMN/Exemples/bingham_simple.c)
and generalizes it for any Power Law fluid (using regularization method). Another
difference between the two is that this code calculates the second invariant of
deformation tensor at the face-centers of the cells instead of the cell centers.

## Mathematical Formulations

Unlike the Newtonian fluids, non-Newtonian fluids do not have a linear stress-strain rate relationship.
One way to represent the relationship is using the Generalized Newtonian fluid method:
$$ \tau = \tauy + 2\mu_0D_{ij}^n $$
The fluid is such that
if $\|\tau\| \le  \tauy$ then there is no motion $D_{ij}=0\:\forall\:(i,j)$
if the stress is high enough $\|\tau\| >  \tauy$ then there is motion

*Note:* that $\|\tau\|$ is the modulus defined as the Euclidian norm  $\sqrt{\frac{1}{2}{\tau_{ij} \tau_{ij}}}$.
 It is not $\sqrt{\tau_{11}^2 + \tau_{12}^2}$ as in Balmorth et al. (2006), which is the Frobenius norm.

$D_{ij}$ is the shear strain rate tensor (or the deformation tensor)

$D_{ij}=(u_{i,j}+u_{j,i})/2$: the components in 2D:
$$D_{11}=\frac{\partial u}{\partial x}$$
$$D_{12} =\frac{1}{2}\left(\frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x}\right)$$
$$D_{21} =D_{12} =\frac{1}{2}\left( \frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x}\right)$$
$$D_{22}=\frac{\partial v}{\partial y}$$

In the Euclidian norm we have:
$$\|D\|=\sqrt{\frac{D_{ij}D_{ij}}{2}}$$
The second invariant defined by $D_2=\sqrt{D_{ij}D_{ij}}$ (this is the Frobenius norm)
is given by:
$$D_2^2= D_{ij}D_{ij}= \left( \frac{\partial u}{\partial x}\right)^2 + \left(\frac{\partial v}{\partial y}\right)^2
 +  \frac{1}{2}\left( \frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x}\right)^2$$
and we have obviously $\|D_{ij}\| = D_2/\sqrt{2}$

## Numerical regularization
$$ \tau_{ij} = \tauy\left(\frac{D_{ij}}{\|D_{ij}\|}\right) + 2\mu_0\|D_{ij}\|^{n-1}D_{ij}^n $$
Factorising with $2D_{ij}$ to obtain a equivalent viscosity
$$\tau_{ij} = 2\left(\mu_0 \|D_{ij}\|^{n-1} + \frac{\tauy}{2 \|D_{ij}\|}\right)D_{ij}$$
$$\tau_{ij} = 2 \mu_{eq}D_{ij}$$
$$\mu_{eq} = \mu_0\|D_{ij}\|^{n-1} + \frac{\tauy}{2\|D_{ij}\|}$$
$\mu$ is the min of $\mu_{eq}$ and a large $\mu_{max}$ so that the viscosity does not blow up.
$$ \mu = \text{min}\left(\mu_{eq}, \mu_{max}\right) $$
*Note:* We present here the formulation in Balmforth, he uses $\dot{\gamma}$ which is by his definition $\sqrt{\frac{1}{2}\dot{\gamma_{ij}}\dot{\gamma_{ij}}}$
and as $\dot{\gamma_{ij}}=2 D_{ij}$ then $\dot{\gamma}$ is $\sqrt{2}D_2$, that is why we have a $\sqrt{2}$ in the equations.


##  Exact solution in the proposed case


We look at an unidirectional flow, a pure shear flow  $u(y)$, $v=0$, so
$D_{11}=D_{22}=0$ and $D_{12}=D_{21}=(1/2)(\partial u/\partial y)$.
$$ \tau_{12} = 2\mu D_{12}^n  + \tauy =
  2^{1-n}\mu\left(\frac{\partial u}{\partial y}\right)^n + \tauy $$

Equilibrium between pressure gradient and viscosity (writting $\tau$ for a shorthand of $\tau_{12}$)
$$0=-\frac{\partial p}{\partial x} + \frac{\partial \tau}{\partial y}$$
as there is no stress at the free surface $y=h$, the stress is
$$ \tau = \left(-\frac{\partial p}{\partial x}\right)(h-y)$$
the stress $\tau$ increases from the free surface, as long as $\tau<\tauy$,
we are under the threshold,
so shear is zero: $\frac{\partial u}{\partial y} =0$,  hence velocity is constant, say it is $U$.
Let us define
 $Y=h-\tauy/(-\frac{\partial p}{\partial x})$, where  $\tau=\tauy$.


So :
$$ \left\{\tau<\tauy, \frac{\partial u}{\partial y} = 0,\:\&\:u=U\:\forall\:Y<y<h\right\} $$
Then going down:
$0<y<Y$ we have $\tau = 2^{1-n}\mu\left(\frac{\partial u}{\partial y}\right)^n + \tauy$.
This gives:
$$\tauy + 2^{1-n}\mu\left(\frac{\partial u}{\partial y}\right)^n =  \left(-\frac{\partial p}{\partial x}\right)(Y-y)$$
After some straight forward manipulations:
$$\left(\frac{\partial u}{\partial y}\right) =  \left(\frac{1}{2^{1-n}\mu}\right)^{1/n}\left(-\frac{\partial p}{\partial x}\right)^{1/n}(Y-y)^{1/n}
 = A_{p\mu}(Y-y)^{1/n}$$
and this allows to solve for the velocity profile
$$u = \frac{nA_{p\mu}}{n+1}\left(Y^{\frac{n+1}{n}} - \left(Y-y\right)^{\frac{n+1}{n}}\right)$$
which is indeed zero in $y=0$, and for  $y=Y$, we have the plug flow $u=U$ of value:
$$U= \frac{n}{n+1}A_{p\mu}Y^{\frac{n+1}{n}}$$
$$A_{p\mu} = \left(\frac{1}{2^{1-n}\mu}\right)^{1/n}\left(-\frac{\partial p}{\partial x}\right)^{1/n}$$
For the present case $-\frac{\partial p}{\partial x} = 1, \mu = 1, h = 1$, which gives:
$A_{p\mu} = \frac{1}{2^{(1-n)/n}}, Y = 1 - \tauy$

# Code
*/
#include "navier-stokes/centered.h"

// Global parameters
char file_name[80];
double tauy = 0.0;        // Yield stress
double mu_0 = 1.0;         // Base viscosity
double mu_max = 1000.;     // Maximum viscosity for regularization
double n = 1.0;            // Power law exponent
int max_iter = 1e4;        // Maximum iterations
#define DT_MAX (1e-3)      // Maximum timestep

/**
 * @brief Main entry point for the planar Couette flow simulation.
 *
 * This function initializes the simulation grid and domain parameters (domain size, origin, timestep, and convergence
 * tolerance), and sets the boundary conditions (periodic for right-left, slip at the top, and no-slip at the bottom).
 * It then sequentially configures and runs simulations for four fluid models—Newtonian, Power law, Herschel-Bulkley, and
 * Bingham—by setting the corresponding parameters (yield stress, base viscosity, and power-law exponent) and invoking
 * the simulation routine.
 *
 * The output file name is initially sourced from the first command-line argument and updated for each fluid case.
 *
 * @return int Exit status of the program.
 */
int main(int argc, char const *argv[])
{
  
  sprintf(file_name, "%s", argv[1]);
  
  // Initialize grid and domain
  init_grid(1<<6);
  L0 = 1.0;
  origin(-0.5, -0.5);
  DT = DT_MAX;
  stokes = true;
  TOLERANCE = 1e-5;

  // Set boundary conditions
  // Right-left boundaries are periodic
  periodic(right);
  
  // Slip at the top
  u.t[top] = neumann(0);
  u.n[top] = neumann(0);
  
  // No slip at the bottom
  u.n[bottom] = dirichlet(0);
  u.t[bottom] = dirichlet(0);

/**
 Values of yield stress, viscosity, and coefficient.
 - Newtonian: $\mu_0 = 1.0$; $\tauy = 0.$ and n = 1
 - Power law: $\mu_0 = 1.0$; $\tauy = 0.$ and n = 0.5
 - Herschel-Bulkley: $\mu_0 = 1.0$; $\tauy = 0.25$ and n = 0.5
 - Bingham: $\mu_0 = 1.0$; $\tauy = 0.25$ and n = 1
*/
  for (int i = 0; i < 4; i++) {
    // default parameters
    tauy = 0.0;
    mu_0 = 1.0;
    n = 1.0;
    // Set parameters based on fluid type
    switch (i) {
      case 0: // Newtonian
        tauy = 0.0;
        mu_0 = 1.0;
        n = 1.0;
        sprintf(file_name, "newtonian");
        break;
      case 1: // Power law
        tauy = 0.0;
        mu_0 = 1.0;
        n = 0.5;
        sprintf(file_name, "powerlaw");
        break;
      case 2: // Herschel-Bulkley
        tauy = 0.25;
        mu_0 = 1.0;
        n = 0.5;
        sprintf(file_name, "herschelbulkley");
        break;
      case 3: // Bingham
        tauy = 0.25;
        mu_0 = 1.0;
        n = 1.0;
        sprintf(file_name, "bingham");
        break;
    }
    fprintf(ferr, "Running case %d: tauy = %g, mu_0 = %g, n = %g\n", 
            i, tauy, mu_0, n);
    run();
  }
}

// un is used to search for a stationary solution
scalar un[];
// muv will be used as the face vector for viscosity
face vector muv[];

/**
## Initialization event
*/
event init(t = 0) {
  // Set Non-Newtonian viscosity
  mu = muv;
  
  /**
   Pressure gradient `mdpdx`
   $$-\frac{dp}{dx} = 1 $$
  */
  const face vector mdpdx[] = {1.0, 0.0};
  a = mdpdx;
  
  // Initialize velocity field at rest
  foreach() {
    u.x[] = 0;
    u.y[] = 0;
    un[] = 0;
  }
  
  dump(file = "start");
}

/**
## Monitoring convergence 
We look for a stationary solution by checking changes in velocity field.
*/
event logfile(i += 500; i <= max_iter) {
  double du = change(u.x, un);
  fprintf(ferr, "i = %d: err = %g\n", i, du);
  
  if (i > 0 && du < 1e-6) {
    dump(file = file_name);
    return 1; /* stop */
  }
  
  if (i == max_iter) {
    dump(file = file_name);
  }
}

/**
 * @brief Computes the effective viscosity for a generalized Newtonian fluid at cell faces.
 *
 * This event updates the face viscosity field by calculating a locally adaptive viscosity based
 * on the rate-of-strain tensor. The procedure involves:
 * - Computing deformation tensor components from the velocity field:
 *   - D11 from consecutive horizontal velocity differences.
 *   - D22 from averaged vertical velocity differences.
 *   - D12 from a combination of horizontal and vertical differences.
 * - Calculating the second invariant of the deformation tensor:
 *   $$D_2 = \\frac{\\sqrt{D_{11}^2 + D_{22}^2 + 2\\,D_{12}^2}}{\\Delta},$$
 *   where \\(\\Delta\\) is the grid spacing.
 * - Determining the equivalent viscosity:
 *   $$\\mu_{eq} = \\mu_0\\left(\\frac{D_2}{\\sqrt{2}}\\right)^{n-1} + \\frac{\\tauy}{\\sqrt{2}\\,D_2},$$
 * - Enforcing an upper limit by taking
 *   $$\\mu = \\min(\\mu_{eq}, \\mu_{max}).$$
 *
 * If the deformation invariant \\(D_2\\) is zero, viscosity is set based on the following conditions:
 * - If the yield stress \\(\\tauy\\) is positive or the exponent \\(n < 1\\), the viscosity defaults to \\(\\mu_{max}\\).
 * - Otherwise, if \\(n = 1\\) it is set to \\(\\mu_0\\); for other cases, it is set to 0.
 *
 * The computed viscosity is scaled by the corresponding face metric and stored in the face viscosity field.
 */
event properties(i++) {
  /**
  Implementation of generalized Newtonian viscosity:
  
  The second invariant is $D_2=\sqrt{D_{ij}D_{ij}}$ (Frobenius norm)
  $$D_2^2 = D_{ij}D_{ij} = D_{11}^2 + 2D_{12}^2 + D_{22}^2$$
  
  The equivalent viscosity is:
  $$\mu_{eq} = \mu_0\left(\frac{D_2}{\sqrt{2}}\right)^{n-1} + 
               \frac{\tauy}{\sqrt{2} D_2}$$
  
  Finally: $\mu = \min(\mu_{eq}, \mu_{max})$
  */
  foreach_face() {
    // Calculate deformation tensor components at face centers
    double D11 = (u.x[] - u.x[-1,0]);
    double D22 = ((u.y[0,1] - u.y[0,-1]) + (u.y[-1,1] - u.y[-1,-1])) / 4.0;
    double D12 = 0.5 * (((u.x[0,1] - u.x[0,-1]) + 
              (u.x[-1,1] - u.x[-1,-1])) / 4.0 + (u.y[] - u.y[-1,0]));
    
    // Calculate second invariant
    double D2 = sqrt(sq(D11) + sq(D22) + 2.0 * sq(D12)) / Delta;
    
    // Calculate effective viscosity
    double mu_temp;
    if (D2 > 0.0) {
      double temp = tauy / (sqrt(2.0) * D2) + 
                   mu_0 * exp((n - 1.0) * log(D2 / sqrt(2.0)));
      mu_temp = min(temp, mu_max);
    } else {
      if (tauy > 0.0 || n < 1.0) {
        mu_temp = mu_max;
      } else {
        mu_temp = (n == 1.0 ? mu_0 : 0.0);
      }
    }
    
    // Apply viscosity at face
    muv.x[] = fm.x[] * mu_temp;
  }
}