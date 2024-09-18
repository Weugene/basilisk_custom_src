/**
# Aslam Extrapolations

This module implements the PDE-based Aslam extrapolation
methods ([Aslam 2003](#aslam2004partial)). Using these
functions, a discontinuous field can be extended from the
liquid phase to the gas-phase and vice-versa using
constant and extrapolations.

These methods are expensive, since the PDE must be solved
at steady-state. However, most of the times they are
used just for a few steps in order to extend the field for
a few cells across the interface, without covering the whole
domain.
*/

#include "mapregion.h"
#include "redistance.h"

void vof_to_ls (scalar f, scalar levelset, int imax = 3) {
  double deltamin = L0/(1 << grid->maxdepth);
  foreach()
    levelset[] = -(2.*f[] - 1.)*deltamin*0.75;
#if TREE
  restriction({levelset});
#endif
  redistance (levelset, imax);
}

/**
## Constant Extrapolation

The constant extrapolation of a field *f*, existing only in
a portion of space defined by *ls*, can be extrapolated in
the whole domain solving the following PDE:

$$
\dfrac{\partial f}{\partial t}
+ H(\phi)\hat{\mathbf{n}}\cdot\nabla f = 0
$$

where $H(\phi)$ is the Heaviside function, used to define the
region where the field should be extrapolated and where it
should not be modified:

$$
H(\phi) =
\begin{cases}
  1 & \text{if } \phi > 0,\\
  0 & \text{if } \phi \leq 0.
\end{cases}
$$

while $\hat{\mathbf{n}}$ is the unit normal.
*/

trace
void constant_extrapolation (
  scalar f,                   // field to extrapolate
  scalar ls,                  // level set field
  double cfl = 0.5,           // CFL number
  int nmax,                   // number of maximum time steps
  (const) scalar s = {-1},    // source term, default zero
  (const) scalar c = {-1},    // vof field, optional
  int nl = 0,                 // from which layer of cells (optional, min=0, max=2)
  int nointerface = 0,        // remove the interface from the heaviside (default false)
  int inverse = 0             // the vof field if = 1 if gas (default false)
)
{
  scalar H[];
  vector n[], gf[];

  if (s.i < 0)
      s[] = 0.;

  /**
  We compute the gradients of the level set for
  the calcation of the interface normals.
  */

  gradients ({ls}, {n});

  /**
  We compute the normals and the heaviside function H,
  which is non-null in the region where the field must
  be extrapolated. In case of vof field, the user can
  decide to extrapolate the field from a layer of cells
  which in not adjacent to the interface. This can be
  specified setting the variables *nl*, 0 by default,
  it can be set to 1 or 2. */

  if (c.i > 0)
    mapregion (H, c, nl=nl, nointerface=nointerface, inverse=inverse);
  else
    foreach()
      H[] = (ls[] <= 0.) ? 0. : 1.;

  foreach() {
    double maggf = 0.;
    foreach_dimension()
      maggf += sq (n.x[]);
    maggf = sqrt (maggf);
    foreach_dimension()
      n.x[] /= (maggf + 1.e-10);
  }

  /**
  We solve the PDE inside a loop over the maximum number
  of time steps that, in principle, should ensure to arrive
  at steady-state. In practice, less time-steps will be
  sufficients to extrapolate the fields over the interface. */

  int ts = 0;

  while (ts < nmax) {

    /**
    The gradients of the extrapolated functions are updated
    using an upwind scheme. */

    foreach()
      foreach_dimension()
        gf.x[] = (n.x[] <= 0.) ? (f[1] - f[])/Delta : (f[] - f[-1])/Delta;

    /**
    We solve a single step of the PDE. */

    foreach() {
      double dt = cfl*Delta;
      double nscalargf = 0.;
      foreach_dimension()
        nscalargf += n.x[]*gf.x[];
      f[] -= dt*H[]*(nscalargf - s[]);
    }
    ts++;
  }
}

/**
## Linear Extrapolation

The linear extrapolation of the field *f* along the normal
direction, can be achieved from the solution of the
following PDE:

$$
\dfrac{\partial f}{\partial t}
+ H(\phi)\left(\hat{\mathbf{n}}\cdot\nabla f - f_n \right) = 0
$$

Which is similar to the equation resolved with the constant
extrapolation procedure, except for the source term $f_n$,
which is the directional derivative of $f$ in the normal
direction, defined as:

$$
f_n = \hat{\mathbf{n}}\cdot \nabla f
$$

Since this function is defined only in the region where $f$
is defined, we apply the constant extrapolation in order to
obtain a field $f_n$ defined on the whole domain:

$$
\dfrac{\partial f_n}{\partial t}
+ H(\phi)\hat{\mathbf{n}}\cdot\nabla f_n = 0
$$

Therefore, the linear interpolation requires the solution of
two different extrapolation PDEs, and the same logic applies
if higher order extrapolations have to be solved.
*/

void linear_extrapolation (
  scalar f,                   // field to extrapolate
  scalar ls,                  // level set field
  double cfl,                 // CFL number
  int nmax,                   // number of maximum time steps
  (const) scalar s = {-1},    // source term, default zero
  (const) scalar c = {-1},    // vof field, optional
  int nl = 0,                 // from which layer of cells (optional, min=0, max=2)
  int nointerface = 0,        // remove the interface from the heaviside (default false)
  int inverse = 0             // the vof field if = 1 if gas (default false)
)
{
  scalar H[], fn[];
  vector n[], gf[];

  if (s.i < 0)
      s[] = 0.;

  /**
  We compute the gradients of the level set for
  the calcation of the interface normals, and the
  gradient of *f* for the calculation of the directional
  derivative. */

  gradients ({ls, f}, {n, gf});

  /**
  We compute the normals and the heaviside function H,
  which is non-null in the region where the field must
  be extrapolated. In case of vof field, the user can
  decide to extrapolate the field from a layer of cells
  which in not adjacent to the interface. This can be
  specified setting the variables *nl*, 0 by default,
  it can be set to 1 or 2. */

  if (c.i > 0)
    mapregion (H, c, nl = (nl == 0.) ? 1. : min (nl+1, 2),
        nointerface=nointerface, inverse=inverse);
  else
    foreach()
      H[] = (ls[]+Delta <= 0.) ? 0. : 1.;

  foreach() {
    double maggf = 0.;
    foreach_dimension()
      maggf += sq (n.x[]);
    maggf = sqrt (maggf);
    foreach_dimension()
      n.x[] /= (maggf + 1.e-10);
  }

  /**
  We compute the directional derivative *fn*. */

  foreach()
    foreach_dimension()
      fn[] += n.x[]*gf.x[];

  /**
  We solve the constant extrapolation for extending
  the directional derivative. */

  int ts = 0;

  while (ts < nmax) {

    /**
    The gradients of the extrapolated functions are updated
    using an upwind scheme. */

    foreach()
      foreach_dimension()
        gf.x[] = (n.x[] <= 0.) ? (fn[1] - fn[])/Delta : (fn[] - fn[-1])/Delta;

    /**
    We solve a single step of the PDE. */

    foreach() {
      double dt = cfl*Delta;
      double nscalargf = 0.;
      foreach_dimension()
        nscalargf += n.x[]*gf.x[];
      fn[] -= dt*H[]*(nscalargf - s[]);
    }
    ts++;
  }

  /**
  We solve a constant extrapolation with source equal to
  the directional derivative. */

  if (c.i > 0)
    constant_extrapolation (f, ls, cfl, nmax, fn, c, nl, nointerface);
  else
    constant_extrapolation (f, ls, cfl, nmax, fn);
}

/**
## Linear Extrapolation with constant value on surface

The linear extrapolation of the field *f* along the normal
direction, can be achieved from the solution of the
following PDE:

$$
\dfrac{\partial f}{\partial t}
+ H(\phi)\left(\hat{\mathbf{n}}\cdot\nabla f - f_n \right) = 0
$$

Which is similar to the equation resolved with the constant
extrapolation procedure, except for the source term $f_n$,
which is the directional derivative of $f$ in the normal
direction, defined as:

$$
f_n = \hat{\mathbf{n}}\cdot \nabla f
$$

Since this function is defined only in the region where $f$
is defined, we apply the constant extrapolation in order to
obtain a field $f_n$ defined on the whole domain:

$$
\dfrac{\partial f_n}{\partial t}
+ H(\phi)\hat{\mathbf{n}}\cdot\nabla f_n = 0
$$

Therefore, the linear interpolation requires the solution of
two different extrapolation PDEs, and the same logic applies
if higher order extrapolations have to be solved.
*/

void linear_extrapolation_constant_surface_value (
    scalar f,                   // field to extrapolate
    scalar ls,                  // level set field
    double f_solid,             // value on interface
    double cfl,                 // CFL number
    int nmax,                   // number of maximum time steps
    (const) scalar s = {-1},    // source term, default zero
    (const) scalar c = {-1},    // vof field, optional but given with fs
    (const) scalar fs = {-1},   // vof field, optional but given with c
    int nl = 0,                 // from which layer of cells (optional, min=0, max=2)
    int nointerface = 0,        // remove the interface from the heaviside (default false)
    int inverse = 0             // the vof field if = 1 if gas (default false)
)
{
    assert (!((c.i < 0) ^ (fs.i < 0)));
    scalar H[], fn[];
    vector n[], gf[];

    if (s.i < 0)
        s[] = 0.;

    /**
    We compute the gradients of the level set for
    the calculation of the interface normals, and the
    gradient of *f* for the calculation of the directional
    derivative. */

    gradients ({ls, f}, {n, gf});

    /**
     We compute normal *n*. */
    foreach() {
        double maggf = 0.;
        foreach_dimension()
        maggf += sq (n.x[]);
        maggf = sqrt (maggf);
        foreach_dimension()
        n.x[] /= (maggf + 1.e-10);
    }

    /**
    We compute the directional derivative *fn*. */

    foreach() {
        fn[] = 0;
        foreach_dimension()
        fn[] += n.x[] * gf.x[];
    }

    /**
    We set f_solid on the interface. */
    foreach() {
        if (c[] > F_ERR && c[] < 1.-F_ERR) {
            f[] = f_solid + Delta*(c[] - 0.5)*fn[];
        }
    }
    boundary({f});

    /**
    We compute the normals and the heaviside function H,
    which is non-null in the region where the field must
    be extrapolated. In case of vof field, the user can
    decide to extrapolate the field from a layer of cells
    which in not adjacent to the interface. This can be
    specified setting the variables *nl*, 0 by default,
    it can be set to 1 or 2. */

    if (c.i > 0)
        mapregion (H, c, nl = (nl == 0.) ? 1. : min (nl+1, 2), nointerface=0, inverse=inverse);
    else
        foreach()
            H[] = (ls[]+Delta <= 0.) ? 0. : 1.;

    /**
    We solve the constant extrapolation for extending
    the directional derivative. */

    int ts = 0;

    while (ts < nmax) {

        /**
        The gradients of the extrapolated functions are updated
        using an upwind scheme. */

        foreach()
            foreach_dimension()
                gf.x[] = (n.x[] <= 0.) ? (fn[1] - fn[]) / Delta : (fn[] - fn[-1]) / Delta;

        /**
        We solve a single step of the PDE. */

        foreach() {
            double dt = cfl*Delta;
            double nscalargf = 0.;
            foreach_dimension()
            nscalargf += n.x[]*gf.x[];
            fn[] -= dt*H[]*(nscalargf - s[]);
            fn[] *= min(1 - c[] + fs[], 1); // min(1 - chi + chi_b, 1)
        }
        ts++;
    }

    /**
    We solve a constant extrapolation with source equal to
    the directional derivative.
    Here we don't modify interface values, therefore `nointerface=1`. */

    if (c.i > 0)
        constant_extrapolation (f, ls, cfl, nmax, fn, c, nl, nointerface=1, inverse=inverse);
    else
        constant_extrapolation (f, ls, cfl, nmax, fn, inverse=inverse);
}

// Define a type for a function pointer that takes coordinates and returns a value
typedef double (*solid_function)(double x, double y, double z);
const scalar src_zero = {-1};

void linear_extrapolate_fields (
    scalar * fields,                    // array of fields to extrapolate
    scalar ls,                          // level set field
    solid_function * f_solid_fun_list,  // array of functions returning the solid value based on coordinates
    double cfl,                         // CFL number
    int nmax,                           // number of maximum time steps
    scalar * s_list = NULL,             // array of source terms, corresponding to fields, default zero
    (const) scalar c = {-1},            // vof field, optional but given with fs
    (const) scalar fs = {-1},           // BPM target vof field, optional but given with c
    int nl = 0,                         // from which layer of cells (optional, min=0, max=2)
    int nointerface = 0,                // remove the interface from the heaviside (default false)
    int inverse = 0                     // the vof field if = 1 if gas (default false)
)
{
    int len_fields = list_len(fields);
    assert (!((c.i < 0) ^ (fs.i < 0)));
    scalar H[], fn[];
    vector n[];

    /**
    We compute the gradients of the level set for
    the calculation of the interface normals, and the
    gradient of *f* for the calculation of the directional
    derivative. */
    scalar * fl = NULL;
    vector * gfl = NULL;
    for (scalar f in fields){
        fl = list_append (fl, f);
        vector gf = new vector;
        foreach_dimension() {
            gf.x.gradient = NULL;
        #if TREE
            gf.x.refine = gf.x.prolongation = refine_linear;  // Linear prolongation on refinement
            gf.x.restriction = restriction_volume_average;  // Volume-averaged restriction on coarsening
        #endif
        }
        gfl = vectors_append(gfl, gf);
    }
    fl = list_append (fl, ls);
    gfl = vectors_append(gfl, n);

    gradients (fl, gfl);
    boundary(gfl);

    /**
     We compute normal *n*. */
    foreach() {
        double maggf = 0.;
        foreach_dimension()
            maggf += sq (n.x[]);
        maggf = sqrt (maggf);
        foreach_dimension()
            n.x[] /= (maggf + 1.e-10);
    }
    boundary((scalar *){n});

    /**
    We compute the normals and the heaviside function H,
    which is non-null in the region where the field must
    be extrapolated. In case of vof field, the user can
    decide to extrapolate the field from a layer of cells
    which in not adjacent to the interface. This can be
    specified setting the variables *nl*, 0 by default,
    it can be set to 1 or 2. */

    if (c.i > 0)
        mapregion (H, c, nl = (nl == 0.) ? 1. : min (nl+1, 2), nointerface=0, inverse=inverse);
    else
        foreach()
            H[] = (ls[] + Delta <= 0.) ? 0. : 1.;

    for (int i = 0; i < len_fields; i++){
        scalar f = fields[i];
        vector gf = gfl[i];
        solid_function f_solid_fun = f_solid_fun_list[i];
        // Check if a corresponding source term is provided

        scalar s = (s_list != NULL && s_list[i].i >= 0) ? s_list[i] : src_zero;
        /**
        We compute the directional derivative *fn*. */
        foreach() {
            fn[] = 0;
            foreach_dimension()
                fn[] += n.x[] * gf.x[];
        }
        boundary({fn});

        /**
        We set f_solid on the interface. */
        foreach() {
            if (c[] > F_ERR && c[] < 1. - F_ERR) {
                f[] = f_solid_fun(x, y, z) + Delta * (c[] - 0.5) * fn[];
            }
        }
        boundary({f});

        /**
        We solve the constant extrapolation for extending
        the directional derivative. */

        int ts = 0;

        while (ts<nmax) {

            /**
            The gradients of the extrapolated functions are updated
            using an upwind scheme. */

            foreach()
                foreach_dimension()
                    gf.x[] = (n.x[] <= 0.) ? (fn[1] - fn[]) / Delta : (fn[] - fn[-1]) / Delta;
            boundary((scalar *){gf});

            /**
            We solve a single step of the PDE. */

            foreach() {
                double dt = cfl * Delta;
                double nscalargf = 0.;
                foreach_dimension()
                    nscalargf += n.x[] * gf.x[];
                fn[] -= dt * H[] * (nscalargf - s[]);
                fn[] *= min(1 - c[] + fs[], 1); // min(1 - chi + chi_b, 1)
            }
            boundary({fn});
            ts++;
        }

        /**
        We solve a constant extrapolation with source equal to
        the directional derivative.
        Here we don't modify interface values, therefore `nointerface=1`. */

        if (c.i > 0)
            constant_extrapolation (f, ls, cfl, nmax, fn, c, nl, nointerface = 1, inverse = inverse);
        else
            constant_extrapolation (f, ls, cfl, nmax, fn, inverse = inverse);
    }
    delete ((scalar *) gfl); free (gfl);
}

/**
## References

~~~bib
@article{aslam2004partial,
  title={A partial differential equation approach to multidimensional extrapolation},
  author={Aslam, Tariq D},
  journal={Journal of Computational Physics},
  volume={193},
  number={1},
  pages={349--355},
  year={2004},
  publisher={Elsevier}
}
~~~
*/

