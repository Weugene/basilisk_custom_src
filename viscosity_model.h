/**
# Variable rheology models
 $$\mu = \mu_1 \exp(\frac{E_\eta}{RT})(\frac{\alpha_{gel}}{\alpha_{gel}-\alpha})^f(\alpha, T)$$
 $$f(\alpha, T) = A + B \alpha$$

 $$\mu = \mu_1 \exp(\frac{E_\eta}{RT}+\chi \alpha^2)$$
 **/
#define VISCOSITY_MODEL_CONSTANT 1
#define VISCOSITY_MODEL_DUSI 2
#define VISCOSITY_MODEL_GELATION 3

// User must define VISCOSITY_MODEL
#ifndef VISCOSITY_MODEL
    #error "VISCOSITY_MODEL is not defined. Please define it as one of the allowed values."
#endif

// Check if REACTION_MODEL is one of the allowed values
#if (VISCOSITY_MODEL != VISCOSITY_MODEL_CONSTANT) && \
    (VISCOSITY_MODEL != VISCOSITY_MODEL_DUSI) && \
    (VISCOSITY_MODEL != VISCOSITY_MODEL_GELATION)
    #error "Invalid VISCOSITY_MODEL value. It must be one of the following: 1 (VISCOSITY_MODEL_CONSTANT), 2 (VISCOSITY_MODEL_DUSI), 3 (VISCOSITY_MODEL_GELATION)."
#endif


double mu_eff = 0; // ????

#ifndef muf1
    #if VISCOSITY_MODEL == VISCOSITY_MODEL_CONSTANT
        #define muf1(alpha_doc, T) (mu1)
    #elif VISCOSITY_MODEL == VISCOSITY_MODEL_DUSI
        double Eeta_by_Rg = 0.1; //Kelvin
        double chi = 1;
        #define muf1(alpha_doc, T) (mu0*exp(Eeta_by_Rg/(T) + chi*(alpha_doc)))
    #elif VISCOSITY_MODEL == VISCOSITY_MODEL_GELATION
        double Eeta_by_Rg = 0.1; //Kelvin
        double chi = 1;
        double alpha_gel = 0.8;
        #define fpol(alpha_doc, T) ((A)*alpha_doc + (B))
        #define muf1(alpha_doc, T) (mu0*exp(Eeta_by_Rg/(T))*pow(alpha_gel/(alpha_gel-(alpha_doc)), fpol((alpha_doc), (T))))
    #endif
#endif

#ifndef mu
    #define mu(f, fs, alpha_doc, T) var_harm(f, fs, muf1(alpha_doc, T), mu2, mu3)
//slow convergence
//#define mu(f, fs, alpha_doc, T) ((1.0 - clamp(fs,0.,1.))*(mu2 + (muf1(alpha_doc, T) - mu2)*clamp(f,0.,1.)) + clamp(fs,0.,1.)*mu3)
//#define mu(f, fs, alpha_doc, T) ((fs == 0)*(mu2 + (muf1(alpha_doc, T) - mu2)*clamp(f,0.,1.)) + clamp(fs,0.,1.)*mu3)
//#define mu(f, fs, alpha_doc, T) (1.0/(  (fs == 0)*( clamp(f,0.,1.)*(1.0/muf1(alpha_doc, T) - 1.0/mu2) + 1.0/mu2) + clamp(fs,0.,1.)/mu3))
// #define mu(f, fs, alpha_doc, T) (1.0/(  (1.0 - clamp(fs,0.,1.))*( clamp(f,0.,1.)*(1.0/muf1(alpha_doc, T) - 1.0/mu2) + 1.0/mu2) + clamp(fs,0.,1.)/mu3))
//#define mu(f, fs, alpha_doc, T) ((1 - clamp(fs,0.,1.))/( clamp(f,0.,1.)*(1.0/muf1(alpha_doc, T) - 1.0/mu2) + 1.0/mu2) + clamp(fs,0.,1.)*mu3 )
#endif