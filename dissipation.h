//
// Created by Weugene on 01.04.2022.
//

struct Dissipation {
    scalar dis;
    vector u;
    face vector mu;
    bool dump_dis; // nullify the source term
};

/**
This function compute source term due to viscous heat generation: 
$q_v = \mu\left( \grad \mathbf{u} + (\grad \mathbf{u})^T\right) : \grad \mathbf{u}$.
The viscous dissipation term $\Phi$ in a three-dimensional incompressible fluid, expressed in LaTeX code, is:
\Phi = \mu \left( 
    2\left(\frac{\partial u}{\partial x}\right)^2 + 
    2\left(\frac{\partial v}{\partial y}\right)^2 + 
    2\left(\frac{\partial w}{\partial z}\right)^2 + 
    \left(\frac{\partial u}{\partial y} + \frac{\partial v}{\partial x}\right)^2 + 
    \left(\frac{\partial u}{\partial z} + \frac{\partial w}{\partial x}\right)^2 + 
    \left(\frac{\partial v}{\partial z} + \frac{\partial w}{\partial y}\right)^2 
\right)
Note: by default it nullifies the source term. Also it works only in case of incompressibility. 
For compressible case it needs to add 
$-\frac23\mu\left( \frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} + \frac{\partial w}{\partial z}\right)$
*/
void dissipation (scalar dis, (const) vector u, (const) face vector mu, bool dump_dis = true)
{
#if TREE
    /* conservative coarse/fine discretisation (2nd order) */
    if (dump_dis){
        foreach () {
            dis[] = 0.;
        }
    }

    foreach_dimension() {
        face vector taux[];
        foreach_face(x) {
            taux.x[] =  2.*mu.x[]*sq((u.x[] - u.x[-1])/Delta);
        }
        #if dimension > 1
        /**
        $\mu_y \left( \frac{\partial u}{\partial y}\right) \left( \frac{\partial u}{\partial y} + \frac{\partial v}{\partial x}\right)$
        mu.y[]*((u.x[] - u.x[0,-1])/Delta) * ((u.x[] - u.x[0,-1])/Delta + ((u.y[1,-1] + u.y[1,0])/2 - (u.y[-1,-1] + u.y[-1,0])/2)/(2*Delta))
        */
        foreach_face(y)
            taux.y[] = mu.y[]*(u.x[] - u.x[0,-1] +
                   (u.y[1,-1] + u.y[1,0])/4. -
                   (u.y[-1,-1] + u.y[-1,0])/4.)
                  *(u.x[] - u.x[0,-1])/sq(Delta);
        #endif
        #if dimension > 2
        /**
        $\mu_z \left( \frac{\partial u}{\partial z}\right) \left( \frac{\partial u}{\partial z} + \frac{\partial w}{\partial x}\right)$
        mu.y[]*((u.x[] - u.x[0,0,-1])/Delta) * ((u.x[] - u.x[0,0,-1])/Delta + ((u.y[1,0,-1] + u.y[1,0,0])/2 - (u.y[-1,0,-1] + u.y[-1,0,0])/2)/(2*Delta))
        */
        foreach_face(z)
            taux.z[] = mu.z[]*(u.x[] - u.x[0,0,-1] +
                   (u.z[1,0,-1] + u.z[1,0,0])/4. -
                   (u.z[-1,0,-1] + u.z[-1,0,0])/4.)
                  *(u.x[] - u.x[0,0,-1])/sq(Delta);
        #endif
        foreach () {
            // average values laying on faces into cell center
            double d = 0;
            foreach_dimension()
                d += (taux.x[1] + taux.x[]);
            dis[] += d/(2.*dimension); // average over all faces
        }
    }
#else
/*   /\* "naive" discretisation (only 1st order on trees) *\/ */
    foreach () {
        if (dump_dis)
            dis[] = 0.;
        foreach_dimension(){
            dis[] += ((mu.x[1,0]*sq(u.x[1] - u.x[])
                  + mu.x[]*sq(u.x[] - u.x[-1])
#if dimension > 1
                  + mu.y[0,1]*(u.x[0,1] - u.x[] +
                (u.y[1,0] + u.y[1,1])/4. -
                (u.y[-1,0] + u.y[-1,1])/4.)*(u.x[0,1] - u.x[])/2.
                + mu.y[]*(u.x[] - u.x[0,-1] +
                (u.y[1,-1] + u.y[1,0])/4. -
                (u.y[-1,-1] + u.y[-1,0])/4.)*(u.x[] - u.x[0,-1])/2.
#endif
#if dimension > 2
                   + mu.z[0,0,1]*(u.x[0,0,1] - u.x[] +
                (u.z[1,0,0] + u.z[1,0,1])/4. -
                (u.z[-1,0,0] + u.z[-1,0,1])/4.)*(u.x[0,0,1] - u.x[])/2.
                - mu.z[]*(u.x[] - u.x[0,0,-1] +
                (u.z[1,0,-1] + u.z[1,0,0])/4. -
                (u.z[-1,0,-1] + u.z[-1,0,0])/4.)*(u.x[] - u.x[0,0,-1])/2.
#endif
                ) )/sq(Delta);
      }
    }
#endif
}