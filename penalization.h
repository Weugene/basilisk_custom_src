const vector zerocf[] = {0.,0.,0.};
#ifdef BRINKMAN_PENALIZATION
    #define fbp (fs[])
    extern scalar fs;
    extern face vector fs_face;
    double eta_s = 1e-15, nu_s = 1, lambda_slip = 0;
    double m_bp = 0;
    (const) scalar a_br = unity, b_br = unity; // useful for Robin BC
    (const) vector U_solid = zerocf;
    #if BRINKMAN_PENALIZATION == 1 //Dirichlet BC
        (const) vector target_U = zerocf, n_sol = zerocf;
        (const) face vector target_Uf = zerof;
        #define PLUS_CONSTANT_BRINKMAN_RHS + (fbp*(target_U.x[])/eta_s)
        #define PLUS_VARIABLE_BRINKMAN_RHS - (fbp*(u.x[]/eta_s))
        #define PLUS_NUMERATOR_BRINKMAN    + 0
        #define PLUS_DENOMINATOR_BRINKMAN  + fbp*dt*(sq(Delta)/eta_s)
    #elif BRINKMAN_PENALIZATION == 2 //Neumann BC
        vector target_U[], n_sol[];
        #define CALC_GRAD
        #define PLUS_CONSTANT_BRINKMAN_RHS  + (fbp*(target_U.x[])/eta_s)
        #define PLUS_VARIABLE_BRINKMAN_RHS  + fbp*dt*( - (scalar_a_by_b(n_sol, grad_u) )/eta_s)
        #define PLUS_NUMERATOR_BRINKMAN   + fbp*dt*( - (0.5*Delta*(m_scalar_a_by_b(n_sol, u)) - sq(Delta)*target_U.x[])/eta_s)
        #define PLUS_DENOMINATOR_BRINKMAN + 0
    #elif BRINKMAN_PENALIZATION == 3 //Robin BC
        vector target_U[], n_sol[];
        #define CALC_GRAD
        #define PLUS_CONSTANT_BRINKMAN_RHS  + (fbp*(target_U.x[])/eta_s)
        #define PLUS_VARIABLE_BRINKMAN_RHS + fbp*dt*( - (a_br[]*u.x[] + b_br[]*scalar_a_by_b(n_sol, grad_u) )/eta_s)
        #define PLUS_NUMERATOR_BRINKMAN   + fbp*dt*( - (0.5*Delta*b_br[]*(m_scalar_a_by_b(n_sol, u)) - sq(Delta)*target_U.x[])/eta_s)
        #define PLUS_DENOMINATOR_BRINKMAN + fbp*dt*(a_br[]*sq(Delta)/eta_s)
    #elif BRINKMAN_PENALIZATION == 4 //No penetration & slip BC
        vector target_U[], n_sol[]; //here target_U will be recalcilated each time step
        #define CALC_GRAD_U_TAU
        #define PLUS_CONSTANT_BRINKMAN_RHS + (fbp*(target_U.x[])/eta_s)
        #define PLUS_VARIABLE_BRINKMAN_RHS - (fbp*(u.x[]/eta_s ))
        #define PLUS_NUMERATOR_BRINKMAN    + 0
        #define PLUS_DENOMINATOR_BRINKMAN  + fbp*dt*(sq(Delta)/eta_s)
//    #elif BRINKMAN_PENALIZATION == 5//Ideal slip, no friction
//        #define CALC_GRAD
//        #define PLUS_BRINKMAN_RHS         + fbp*dt*(nu_s*(u.x[1] - 2*u.x[] +u.x[-1])/sq(Delta) - n_sol.x[]*scalar_a_by_b(n_sol, umUs)/eta_s)
//        #define PLUS_NUMERATOR_BRINKMAN   + fbp*dt*(nu_s*(u.x[1] + u.x[-1]) - sq(Delta)*n_sol.x[]*(scalar_a_by_b(n_sol, umUs) - n_sol.x[]*u.x[])/eta_s)
//        #define PLUS_DENOMINATOR_BRINKMAN + fbp*dt*(sq(n_sol.x[]*Delta)/eta_s + 2*nu_s)
    #else
        #define BRINKMAN_PENALIZATION_ERROR_BC 1
    #endif
    #undef SEPS
    #define SEPS 1e-15
    #ifdef DEBUG_BRINKMAN_PENALIZATION
        vector dbp[], total_rhs[], residual_of_u[];
        vector utau[], grad_utau_n[];
        #define gradun grad_utau_n.x[]
    #else
        double gradun;
    #endif
#else
    #define fbp 0
    double gradun;
    #define PLUS_CONSTANT_BRINKMAN_RHS  0
    #define PLUS_VARIABLE_BRINKMAN_RHS 0
    #define PLUS_NUMERATOR_BRINKMAN 0
    #define PLUS_DENOMINATOR_BRINKMAN 0
#endif

#define PLUS_BRINKMAN_RHS  (PLUS_VARIABLE_BRINKMAN_RHS + PLUS_CONSTANT_BRINKMAN_RHS) //- fbp*( (u.x[] - target_U.x[])/eta_s )

struct Brinkman {
    vector u;
    face vector uf;
    scalar rho;
    double dt;
};


double give_etas(double m_bp, double mindelta, double nu_min){
    return sq(m_bp * mindelta) / nu_min;
}

double give_mbp(double eta_s, double mindelta, double nu_min){
    return sqrt(eta_s * nu_min) / mindelta;
}

/**
 * The Brinkman penalization method is used to simulate the flow around solid obstacles.
 * The method is based on the introduction of a penalization term in the momentum equation.
 * The penalization term depends on the kinematic viscosity of the fluid $\nu=\mu_1/\rho_1$ ,
 * the penalization coefficient eta_s, and the penalization parameter m_bp - the number of mesh cells
 * used to resolve the Brinkman penalization layer $\sqrt{\eta_s \nu}$.
 */
void set_penalization_parameters (face vector mu, scalar rho, double new_m_bp, double new_eta_s){
    int maxlevel = grid->maxdepth;
    double mindelta = L0 / (1 << maxlevel);
    double nu_min = 1e+10;
    foreach( reduction(min:nu_min) ){
        double nu = norm(mu) / rho[];
        if (nu < nu_min) nu_min = nu;
    }

    if (nu_min > SEPS) {
        if (fabs(new_m_bp) > 0) { // m_bp has higher priority
            m_bp = new_m_bp;
            eta_s = give_etas(m_bp, mindelta, nu_min);
        } else if (fabs(new_m_bp) == 0 && fabs(eta_s) < SEPS) { // nothing is set, m_bp = 1, eta_s(m_bp)
            m_bp = 1;
            eta_s = give_etas(m_bp, mindelta, nu_min);
        } else { // only eta is set
            eta_s = new_eta_s;
            m_bp = give_mbp(eta_s, mindelta, nu_min);
        }
        fprintf(
            ferr,
            "Brinkman penalization params for u: eta_s=%g, m_bp=%g, minDelta=%g, nu_min=%g\n",
            eta_s, m_bp, mindelta, nu_min
        );
    }
}
