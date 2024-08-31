
#if BRINKMAN_PENALIZATION == 4
void calc_target_U(const vector u, vector target_U, const vector normal){
    if (!is_constant(U_solid.x)) foreach() foreach_dimension() target_U.x[] = U_solid.x[];
    if (fabs(lambda_slip) > 0.) {
        double ubyn;
        #ifndef DEBUG_BRINKMAN_PENALIZATION
        vector utau[]; // otherwise utau will be defined globally
        #endif
        if (!is_constant(U_solid.x)) foreach() foreach_dimension() u.x[] -= U_solid.x[];
        //    if (!is_constant(U_solid.x)) fprintf(ferr, "U_solid.x");
        foreach() {
            ubyn = 0;
            foreach_dimension() ubyn += u.x[]*normal.x[];
            foreach_dimension() utau.x[] = u.x[] - ubyn*normal.x[];
        }
        if (!is_constant(U_solid.x)) foreach() foreach_dimension() u.x[] += U_solid.x[];

        foreach() {
            if (0 < fs[] && fs[] < 1) { //See here!
				coord gradun;
                foreach_dimension() {
					gradun.x = (n_sol.x[] > 0) ? (utau.x[2] - utau.x[1])/Delta : (utau.x[-1] - utau.x[-2])/Delta;
                    gradun = ((utau.x[-1] - utau.x[-2])*0.5*(fs[-1] + fs[-2]) + (utau.x[2] -  utau.x[1])*0.5*(fs[1] + fs[2]))*n_sol.x[] / Delta
                            #if dimension > 1
                            + ((utau.x[0,-1] - utau.x[0,-2])*0.5*(fs[0,-1] + fs[0,-2]) + (utau.x[0,2] - utau.x[0,1])*0.5*(fs[0,1] + fs[0,2]))*n_sol.y[]/Delta
                            #endif
                            #if dimension > 2
                            + ((utau.x[0, 0, -1] - utau.x[0, 0, -2])*0.5*(fs[0, 0, -1] + fs[0, 0, -2]) + (utau.x[0, 0, 2] - utau.x[0, 0, 1])*0.5*(fs[0, 0, 2] + fs[0, 0, 1]))*n_sol.z[] / Delta
                            #endif
                            ;
                    target_U.x[] += lambda_slip*gradun;
                }
                //            if (target_U.x[]) fprintf(ferr, "Ut = %g, U_s = %g, l=%g gr=%g\n", target_U.x[], U_solid.x[], lambda_slip, gradun);
            }else{
                foreach_dimension() target_U.x[] = u.x[];
            }
        }
    }
}
#endif
void brinkman_correction_u (vector u, double dt){
#if BRINKMAN_PENALIZATION == 4
    if (!is_constant(target_U.x)) calc_target_U(u, target_U, n_sol);
#endif
    foreach() {
        foreach_dimension(){
            u.x[] = (u.x[] + (fbp*dt/eta_s)*target_U.x[])/(1. + fbp*dt/eta_s);
        }
    }
    boundary ((scalar *){u});
}
static int i_bpm=0;
void brinkman_correction_uf (face vector uf){

    //sticky way. Everything near 1 cell is fixed at the surface (or has velocity of the solid)
#if STICKY_SOLID == 1
    if (i_bpm==0) {fprintf(ferr, "Sticky solid. uf.x is corrected even for a partial solid inclusion.\n"); i_bpm++;}
    foreach_face()
        if ((uf.x[] - target_Uf.x[]) && (fs_face.x[] > 0))
            uf.x[] = target_Uf.x[];
    //Not sticky way. Only cells inside of solid has velocity target_Uf
#elif NO_STICKY_SOLID == 1
    if (i_bpm==0) {fprintf(ferr, "No sticky solid. uf.x is corrected only inside solids.\n"); i_bpm++;}
    foreach_face()
    if ((uf.x[] - target_Uf.x[]) && (fs_face.x[] <= 0))
        uf.x[] = target_Uf.x[];
#else //LINEAR_STICKINESS
    if (i_bpm==0) {fprintf(ferr, "Linear stickiness. uf.x is corrected proportionally to the fs.\n"); i_bpm++;}
    foreach_face() {
        uf.x[] = (1 - fs_face.x[])*uf.x[] + fs_face.x[]*target_Uf.x[];
    }
#endif
}