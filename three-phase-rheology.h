/**
# Three-phase interfacial flows

This file helps setup simulations for flows of two fluids separated by
an interface (i.e. immiscible fluids) with solid obstacles. It is typically used in
combination with a [Navier--Stokes solver](navier-stokes/centered.h).

The interface between the fluids is tracked with a Volume-Of-Fluid
method. The volume fraction $f$ and $fs$ leads to next interpretation of averaged value of A =(\rho, \mu, \kappa)
\begin{table}[]
\begin{tabular}{lll}
f & fs & A    \\
1 & 1  & A\_3 \\
0 & 1  & A\_3 \\
1 & 0  & A\_1 \\
0 & 0  & A\_2
\end{tabular}
\end{table}
The above definition of variables leads to a specific definition of the averaged value A.
The densities and dynamic viscosities for fluid 1 and 2 are *rho1*,
*mu1*, *rho2*, *mu2*, respectively, in a solid rho3 and mu3.

To use this module it needs to define the following variables:
- REACTION_MODEL: the reaction model
- HEAT_TRANSFER: the heat transfer model
- IGNORE_SOLID_MU: if defined, the solid phase will be ignored and properties are defined as a function of the liquid phase only.
*/

#include "vof.h"
scalar f[], fs[];
scalar T[];

#if REACTION_MODEL != NO_REACTION_MODEL
    scalar alpha_doc[];
#endif
scalar * interfaces = {f};


double mu0 = 0, rho1 = 1., mu1 = 0., rho2 = 1., mu2 = 0., rho3 = 1., mu3 = 0.;
double kappa1 = 0, kappa2 = 0, kappa3 = 0;  // W/(m*K)
double Cp1 = 0, Cp2 = 0, Cp3 = 0;
int N_smooth = 1;

/**
Auxiliary fields are necessary to define the (variable) specific
volume $\alpha=1/\rho$ as well as the cell-centered density. */

face vector alphav[];
face vector alphamv[];
scalar rhov[];

#ifdef HEAT_TRANSFER
    scalar rhoCpv[];
#endif

event defaults (i = 0) {
    alpha = alphav;
    rho = rhov;
    alpham = alphamv;
    /**
    If the viscosity and conductivity are non-zero, we need to allocate the face-centered
    viscosity and conductivity fields. */

    if (mu1 || mu2)
        mu = new face vector;
#ifdef HEAT_TRANSFER
    if (kappa1 || kappa2)
        kappa = new face vector;
    rhoCp = rhoCpv;
#endif
    /**
    We add the interface to the default display. */

    display ("draw_vof (c = 'f');");
}

/**
 * The following functions are used to define the property mean value: arithmetic (homogeneous) and harmonic.
 */

#ifndef var_hom
    #ifdef IGNORE_SOLID_MU
        #define var_hom(f, fs, A1, A2, A3) ((A2) + ((A1) - (A2))*clamp(f,0.,1.))
    #else
        #define var_hom(f, fs, A1, A2, A3) ((1.0 - clamp(fs,0.,1.))*((A2) + ((A1) - (A2))*clamp(f,0.,1.)) + clamp(fs,0.,1.)*(A3))
    #endif
#endif

#ifndef var_harm
    #ifdef IGNORE_SOLID_MU
        #define var_harm(f, fs, A1, A2, A3) ( (A1)*(A2)/( clamp(f,0.,1.) * ((A2) - (A1)) + (A1) ) )
    #else
        #define var_harm(f, fs, A1, A2, A3) (1.0/(  (1.0 - clamp(fs,0.,1.))*( clamp(f,0.,1.)*(1.0/(A1) - 1.0/(A2)) + 1.0/(A2)) + clamp(fs,0.,1.)/(A3)))
    #endif
#endif

/**
The density, density * heat capacity are defined using arithmetic averages by default.
*/

#ifndef rho
#define rho(f, fs) var_hom(f, fs, rho1, rho2, rho3)
#endif

#ifndef rhoCp
#define rhoCp(f, fs) var_hom(f, fs, rho1*Cp1, rho2*Cp2, rho3*Cp3);
#endif

/**
Meanwhile, the viscosity and thermal conductivity are defined as a harmonic mean by default.
The user can overload these definitions to use other types of averages.
Usually, it is assumed that mu1 is variable, mu2 and mu3 are constant.

In such systems with a large difference in property values, the common
property (effective viscosity or conductivity) is determined by the phase with the
lower value. The harmonic mean is more appropriate because it gives more weight
to the lower value, accurately reflecting how the system behaves.
 */

#ifndef kappa
    #define kappa(f, fs) var_harm(f, fs, kappa1, kappa2, kappa3)
#endif

#include "viscosity_model.h"

/**
We have the option of using some "smearing" of the density/viscosity
jump. */

#ifdef FILTERED
    scalar sf1[], sf2[];
#else
    #define sf1 f
    #define sf2 fs
#endif


event tracer_advection (i++)
{

    /**
    When using smearing of the density jump, we initialise *sf* with the
    vertex-average of *f*. */
#ifdef FILTERED
    scalar sf_s[];
    filter_scalar(f, sf1);
    for (int i_smooth=2; i_smooth<=N_smooth; i_smooth++){
        filter_scalar(sf1, sf_s);
        foreach() sf1[] = sf_s[];
        boundary({sf1});
    }
    filter_scalar(fs, sf2);
    for (int i_smooth=2; i_smooth<=N_smooth; i_smooth++){
        filter_scalar(sf2, sf_s);
        foreach() sf2[] = sf_s[];
        boundary({sf2});
    }

#endif // FILTERED

#ifndef sf1
    filter_scalar(f, sf1);
#endif // !sf1

#ifndef sf2
    filter_scalar(fs, sf2);
#endif // !sf2
#if TREE
    sf1.prolongation = refine_bilinear;
    sf2.prolongation = refine_bilinear;
    sf1.dirty = true; // boundary conditions need to be updated
    sf2.dirty = true; // boundary conditions need to be updated
#endif

}

#include "fractions.h"

event properties (i++) {
    foreach_face() {
        double ff1 = face_value (sf1, 0); // liquid fraction on face
        double ff2 = face_value (sf2, 0); // solid fraction on face
        alphav.x[] = fm.x[]/rho(ff1, ff2);
        alphamv.x[] =  alphav.x[] /(1.0 + ff2*dt/eta_s); // for modified Chorin's method
        if (mu1 || mu2) {
            face vector muv = mu;
            double Tf = face_value (T, 0); // temperature on face
#if REACTION_MODEL != NO_REACTION_MODEL
            double alpha_doc_f = face_value (alpha_doc, 0); // degree of cure on face
            muv.x[] = fm.x[]*mu(ff1, ff2, alpha_doc_f, Tf);
#else
            muv.x[] = fm.x[]*mu(ff1, ff2, 0, Tf);
#endif
        }
#ifdef HEAT_TRANSFER
        if (kappa1 || kappa2) {
            face vector kappav = kappa;
            kappav.x[] = fm.x[]*kappa(ff1, ff2);
        }
#endif
    }
    foreach(){
        rhov[] = cm[] * rho(sf1[], sf2[]);
#ifdef HEAT_TRANSFER
        rhoCpv[] = cm[] * rhoCp(sf1[], sf2[]);
#endif
    }

#if TREE
    sf1.prolongation = fraction_refine; //after changing we restore prolongation operator
    sf2.prolongation = fraction_refine;
    sf1.dirty = true; // boundary conditions need to be updated
    sf2.dirty = true; // boundary conditions need to be updated
#endif
}
