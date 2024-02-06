#ifndef SMALL_EPS
    #define SMALL_EPS 1.e-30
#endif

void MinMaxValues(scalar * list, double * arr_eps) {// for each scalar min and max
    double arr[10][2], small_val = 1e-10;
    int ilist = 0;
    for (scalar s in list) {
        double mina= HUGE, maxa= -HUGE;
        foreach( reduction(min:mina) reduction(max:maxa) ){
            if (fabs(s[]) < mina) mina = fabs(s[]);
            if (fabs(s[]) > maxa) maxa = fabs(s[]);
        }
#if _MPI
        MPI_Allreduce (MPI_IN_PLACE, &mina, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce (MPI_IN_PLACE, &maxa, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
        arr[ilist][0] = mina;
        arr[ilist][1] = maxa;
        ilist++;
//        fprintf(stderr, "arr for i=%d", ilist);
    }
    int i = 0;
    for (scalar s in list){

#if EPS_MAXA == 1
        if (arr[i][1] > small_val){
            arr_eps[i] *=arr[i][1];
        }
#elif EPS_MAXA == 2
        if (arr[i][1] - arr[i][0] > small_val){
            arr_eps[i] *= arr[i][1] - arr[i][0];
        }
#else
        if (0.5*(arr[i][0] + arr[i][1]) > small_val){
            arr_eps[i] *= 0.5*(arr[i][0] + arr[i][1]);
        }
#endif
        else{
            arr_eps[i] *= 1;
        }
#ifdef DEBUG_MINMAXVALUES
        fprintf(stderr, "MinMaxValues: name=%s, min=%g, max=%g, eps=%g\n", s.name, arr[i][0], arr[i][1], arr_eps[i]);
#endif
        i++;
    }
}

int count_cells(double t, int i){
    int tnc = 0, maxlev = 0;
    foreach( reduction(+:tnc) reduction(max:maxlev) ){
        tnc++;
        if (level > maxlev) maxlev = level;
    }
#if _MPI
    int nc = 0;
    foreach(serial, noauto){
        nc++;
    }
    int rank, h_len;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Get_processor_name(hostname, &h_len);
    printf("i %d t %g hostname %s rank %d num cells %d total num cells %d compression rate %g\n", i, t, hostname, rank, nc, tnc, pow(2, dimension*maxlev)/tnc);
#else
    printf("i %d t %g total num cells %d\n", i, t, tnc);
#endif
    fflush(stdout);
    return tnc;
}
// statistical values inside cells with liquid
stats statsf_weugene (scalar f, scalar fs)
{
    double min = 1e100, max = -1e100, sum = 0., sum2 = 0., volume = 0.;
    foreach(reduction(+:sum) reduction(+:sum2) reduction(+:volume)
    reduction(max:max) reduction(min:min)){
        double dvr = dv()*(1. - fs[]);
        double val = f[]*(1. - fs[]);
        volume += dvr;
        sum    += f[]*dvr;
        sum2   += sq(f[])*dvr;
        if (val > max) max = val;
        if (val < min) min = val;
    }
	if (volume > 0.){
    	sum /= volume; sum2 /= volume;
	}
	sum2 -= sq(sum);
	fprintf(ferr, "***: %g %g\n", sum, sum2);
	stats s;
	s.min = min, s.max = max, s.sum = sum, s.volume = volume; //modified by Weugene
	s.stddev = sum2 > 0. ? sqrt(sum2) : 0.;
    return s;
}

// statistical values inside pure liquid
stats statsf_weugene2 (scalar f, scalar fs)
{
    double min = 1e100, max = -1e100, sum = 0., sum2 = 0., volume = 0.;
    foreach(reduction(+:sum) reduction(+:sum2) reduction(+:volume)
    reduction(max:max) reduction(min:min))
    if (fs[] == 0.) {
        double dvr = dv()*(1. - fs[]);
        volume += dvr;
        sum    += dvr*f[];
        sum2   += dvr*sq(f[]);
        if (f[] > max) max = f[];
        if (f[] < min) min = f[];
    }
	if (volume > 0.){
    	sum /= volume; sum2 /= volume;
	}
	sum2 -= sq(sum);
    fprintf(ferr, "sum=%g\n", sum);
    stats s;
    s.min = min, s.max = max, s.sum = sum, s.volume = volume;
    return s;
}

norm normf_weugene (scalar f, scalar fs)
{
    double avg = 0., rms = 0., max = 0., volume = 0.;
    foreach(reduction(max:max) reduction(+:avg)
    reduction(+:rms) reduction(+:volume)){
        double dvr = dv()*(1. - fs[]);
        double v = fabs(f[])*(1. - fs[]);
        if (v > max) max = v;
        volume += dvr;
        avg    += dvr*v;
        rms    += dvr*sq(v);
    }
    norm n;
    n.avg = volume ? avg/volume : 0.;
    n.rms = volume ? sqrt(rms/volume) : 0.;
    n.max = max;
    n.volume = volume;
    return n;
}

double change_weugene (scalar s, scalar sn, scalar fs)
{
    double max = 0.;
    foreach(reduction(max:max)) {
        double ds = fabs (s[] - sn[])*(1. - fs[]);
        if (ds > max)
            max = ds;
        sn[] = s[];
    }
    return max;
}
/**
 * Smoothing function
 * f insput scalar field
 * sf - output smoothed scalar field
 */
void filter_scalar(scalar f, scalar sf){
#if dimension <= 2
    foreach()
    sf[] = (4.*f[] +
            2.*(f[0,1] + f[0,-1] + f[1,0] + f[-1,0]) +
            f[-1,-1] + f[1,-1] + f[1,1] + f[-1,1])/16.;
#else // dimension == 3
    foreach()
        sf[] = (8.*f[] +
            4.*(f[-1] + f[1] + f[0,1] + f[0,-1] + f[0,0,1] + f[0,0,-1]) +
            2.*(f[-1,1] + f[-1,0,1] + f[-1,0,-1] + f[-1,-1] +
                f[0,1,1] + f[0,1,-1] + f[0,-1,1] + f[0,-1,-1] +
                f[1,1] + f[1,0,1] + f[1,-1] + f[1,0,-1]) +
            f[1,-1,1] + f[-1,1,1] + f[-1,1,-1] + f[1,1,1] +
            f[1,1,-1] + f[-1,-1,-1] + f[1,-1,-1] + f[-1,-1,1])/64.;
#endif
#if TREE
    sf.prolongation = refine_bilinear;
    sf.dirty = true; // boundary conditions need to be updated
#endif
}

void filter_scalar_N_times(scalar f, scalar sf, int N_smooth){
    scalar sf_s[];
    filter_scalar(f, sf);
    for (int i_smooth=2; i_smooth<=N_smooth; i_smooth++){
        filter_scalar(sf, sf_s);
        foreach() sf[] = sf_s[];
    }
}

/**
Calculate scalar from face vector. compute viscosity in a cell
*/

void calc_scalar_from_face(const face vector vf, scalar vs){
    foreach() {
        double vsum = 0;
        foreach_dimension() {
            vsum += vf.x[] + vf.x[1];
        }
        vs[] = vsum/(2.0*dimension);
    }
    boundary((scalar *){vs});
}

/**
A function to rescale normals so that they are unit vectors w.r.t. the
2-norm (by default, the 1-norm is adopted for efficiency purposes). */
coord normalize_coord(coord n){
    double nn = SMALL_EPS;
    foreach_dimension() nn += sq(n.x); // (sqrt(sq(nf.x[]) + sq(nf.y[])))
    nn = sqrt(nn);
    foreach_dimension() n.x /= nn;
    return n;
}

/**
 * Compute sign function for a given (x,y) for cylinder obstacles.
 * Inner area has sign > 0, outer one has < 0.
 * @param x input x coordinate
 * @param y input y coordinate
 * @param xc x coordinate of a cylinder center
 * @param yc y coordinate of a cylinder center
 * @param size radius of a cylinder
 * @return sign function
 */
double sign_function_of_cylinder_obstacle(double x, double y, double xc, double yc, double size){
    return sq(size) - sq(x  - xc) - sq(y - yc); // inner distance > 0, outer distance < 0
}


/**
 * Compute sign function for a given (x,y,z) for a bubble in 2D space.
 * Inner area has sign < 0, outer one has > 0.
 * @param x input x coordinate
 * @param y input y coordinate
 * @param xc x coordinate of a bubble center
 * @param yc y coordinate of a bubble center
 * @param size radius of a bubble
 * @return sign function
 */
double sign_function_of_bubble(double x, double y, double xc, double yc, double size){
    return sq(x  - xc) + sq(y - yc) - sq(size); // inner distance > 0, outer distance < 0
}

/**
 * Compute sign function for a given (x,y,z) for a bubble.
 * Inner area has sign < 0, outer one has > 0.
 * @param x input x coordinate
 * @param y input y coordinate
 * @param z input z coordinate
 * @param xc x coordinate of a bubble center
 * @param yc y coordinate of a bubble center
 * @param zc z coordinate of a bubble center
 * @param size radius of a bubble
 * @return sign function
 */
double sign_function_of_bubble(double x, double y, double z, double xc, double yc, double zc, double size){
    return sq(x  - xc) + sq(y - yc) + sq(z - zc) - sq(size); // inner distance > 0, outer distance < 0
}

/**
 * This function compute scalar field from analytical representation of bunch of cylinders.
 * @param f volume fraction
 * @param ns number of cylinders
 * @param generate_cylinders function to generate cylinders centers (xc, yc) and radii (R)
 * @param obstacle_pattern distance function for each cylinder
 */
void compute_volume_fraction (
        scalar f,
        const int ns,
        void (* generate_cylinders) (double xc[], double yc[], double R[], const int ns),
        double (*distance_function) (double x, double y, double xc, double yc, double size)
)
{
    face vector face_f[];
    vertex scalar phi[];
    double xc[ns], yc[ns], R[ns];
    generate_cylinders(xc, yc, R, ns);

    foreach_vertex() {
        phi[] = HUGE;
        /**
        Since the medium is periodic, we need to take into account all
        the disk images using periodic symmetries. */

        for (double xp = -L0; xp <= L0; xp += L0)
            for (double yp = -L0; yp <= L0; yp += L0)
                for (int i = 0; i < ns; i++)
                    for (int i = 0; i < ns; i++)
                        phi[] = intersection (phi[], distance_function(x, y, xc[i] - xp, yc[i] - yp, R[i]));
//		phi[] = -phi[];
    }
    fractions (phi, f, face_f);
}

/**
 * function to compute vector for each cell of scalar field `f`.
 * @param f scalar field
 * @param normal_vector computed normals
 */
void compute_normal(const scalar f, vector normal_vector) {
    foreach () {
        coord n = interface_normal(point, f);
        foreach_dimension() { normal_vector.x[] = n.x; }
    }
}

/**
 * Compute theoretical norms of cylinder on the interface. Out of interface the norm is {0,0,0}.
 * @param f input scalar field
 * @param f_normal output vector field of normals
 * @param ns number of cylinders
 * @param generate_cylinders function to generate cylinders centers (xc, yc) and radii (R)
 */
void calc_norms_theoretical(
        const scalar f,
        vector f_normal,
        const int ns,
        void (* generate_cylinders) (double xc[], double yc[], double R[], const int ns),
){
    double xc[ns], yc[ns], R[ns];
    generate_cylinders(xc, yc, R, ns);
    foreach() {
        if (f[] > 0 && f[] < 1){
            for (int i = 0; i < ns; i++) {
                double magn = sqrt(sq(x - xc[i]) + sq(y - yc[i])) + 1e-15;
                f_normal.x[] = (x - xc[i]) / magn;
                f_normal.y[] = (y - yc[i]) / magn;
                f_normal.z[] = 0;
            }
        }else{
            foreach_dimension() f_normal.x[] = 0;
        }
    }
}