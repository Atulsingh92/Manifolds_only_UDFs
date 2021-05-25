#include "udf.h"
#include <math.h>

#define in1 36
#define in2 35
#define in3 34
#define in4 33
#define out1 51
#define out2 50
#define out3 49
#define out4 48

#define frelax 0.0001

#define L 672
#define Dia 6
#define mu 2.18e-5
#define PI 3.14159
#define density 0.9462
/*
Step 1: integrate the mass flux at every outlet and set value on corresponding BC
Step 2: integrate the static pressure at every inlet and set value on corresponding BC
*/

DEFINE_PROFILE(massflux_inlet, thread, position)
{
	/* Integrates mass flux through upstream outflow interface
	   Sets the mass flux profile on the downstream interface
	 */

	int ID; /* the ID of zone with upstream boundary */
	face_t f;
	real NV_VEC(area);

	real Amagnitude, Flux, mflux_glob;
	real mdot_sum, mdot_glob;
	real A_sum, A_glob;

	face_t f_upstream;
	Domain* domain_upstream;          /* domain is declared as a variable   */
	domain_upstream = Get_Domain(1);  /* returns fluid domain pointer       */

	/* Find the mass flux through the upstream interface outlet */
	/* set the ID of the upstream outlet interface surface */

	int zone_ID = THREAD_ID(thread);
		if (zone_ID == in1)
	{
		ID = out1;
	}
	else if (zone_ID == in2)
	{
		ID = out2;
	}
	else if (zone_ID == in3)
	{
		ID = out3;
	}
	else if (zone_ID == in4)
	{
		ID = out4;
	}
	else
	{
		ID = 10000;
	}

	if (ID != 10000)
	{
		Thread* thread_upstream = Lookup_Thread(domain_upstream, ID);
		Flux  = 0.0;
		Amagnitude = 0.0;
		begin_f_loop(f_upstream, thread_upstream)
			if PRINCIPAL_FACE_P(f_upstream, thread_upstream)
			{
				F_AREA(area, f_upstream, thread_upstream);
				Amagnitude += NV_MAG(area);
				Flux += F_FLUX(f_upstream, thread_upstream);
			}
		end_f_loop(f_upstream, thread_upstream)

		mdot_glob = PRF_GRSUM1(Flux*Amagnitude);
		A_glob = PRF_GRSUM1(Amagnitude);
		mflux_glob = mdot_glob / A_glob; 
	}
	else
	{
		mflux_glob = 0.0;
	}

	begin_f_loop(f, thread)
	{	
		F_PROFILE(f, thread, position) = mflux_glob;
	}
	end_f_loop(f, thread)
	Message0("massflowrate at %d is %g\n", zone_ID, mflux_glob);
}


/*
double calulateReynolds(Domain *domain , int id)
{
	real uc, vc, wc, reynolds;
    Thread *t;
    int cells = 0;
    cell_t c;
    face_t f;
    real Amag = 0.0;
    real Amag_sum;
    real vmag_glob, flux=0.0;
    t = Lookup_Thread(domain, id);
    real NV_VEC(area);

    real vmag, vmag_sum, vmag_total, rho, mu_local;
    
    begin_f_loop(f,t)
    {
		
        uc = F_U(f,t); 
        vc = F_V(f,t); 
        wc = F_W(f,t);     
        rho = density;
        vmag += sqrt(pow(uc,2) + pow(vc,2) + pow(wc,2));
        flux += F_FLUX(f,t);
        cells += 1;
        Amag += NV_MAG(area);
    }
    end_f_loop(f,t)
       
   Amag_sum  = PRF_GRSUM1(Amag);
   vmag_sum = PRF_GRSUM1(vmag*Amag_sum);
    
    /*vmag_glob = vmag_sum / Amag_sum;
    vmag_glob = PRF_GRSUM1(flux*Amag_sum)/(Amag_sum*rho* PI *0.25 * pow(Dia,2));
    
    reynolds = rho * vmag_glob * Dia / mu ;

    return reynolds;
    Message0("vmag_glob is %g\n", vmag_glob);
}
*/

DEFINE_PROFILE(pressure_outlet, thread, position)
{
	/* Averages static pressure on downstream interface (in:*) inflow BC and applies it at the upstream outflow interface  (out:*)
	 */

	int ID; /* the ID of zone with downstream boundary */
	face_t f;
	real NV_VEC(area);

	real Amagnitude, pressure, press_avg;
	real PA_sum, PA_glob;
	real A_sum, A_glob, velocity_glob, density_sum, density_glob ;
	real mdot_sum, mdot_glob, Flux, mflux_glob , reynolds, velocity, velocity_sum, friction_factor;

	face_t f_downstream;
	Domain* domain_downstream;         
	domain_downstream = Get_Domain(1); 
	
	/* set the ID of the downstream interface */
	int zone_ID = THREAD_ID(thread); /* applied to */

		if (zone_ID == out1)
	{
		ID = in1;
	}
	else if (zone_ID == out2)
	{
		ID = in2;
	}
	else if (zone_ID == out3)
	{
		ID = in3;
	}
	else if (zone_ID == out4)
	{
		ID = in4;
	}
	else
	{
		ID = 10000;
	}

	if (ID!=10000)
	{		
		Thread* thread_downstream = Lookup_Thread(domain_downstream, ID);
		Amagnitude = 0.0; 
		Flux = 0.0;
		pressure = 0.0;
		
		begin_f_loop(f_downstream, thread_downstream)
			if PRINCIPAL_FACE_P(f_downstream, thread_downstream)
			{
				F_AREA(area, f_downstream, thread_downstream);
				Amagnitude += NV_MAG(area);
				Flux += F_FLUX(f_downstream, thread_downstream);
				pressure +=  F_P(f_downstream, thread_downstream); /*static pressure*/
				/*velocity_squared += pow(F_U(f_downstream, thread_downstream),2)+pow(F_V(f_downstream, thread_downstream),2)+pow(F_W(f_downstream, thread_downstream),2) * pow(NV_MAG(area),2);*/
			}
		end_f_loop(f_downstream, thread_downstream)

		/* If executing in parallel we need a global summation across processors: */
		mdot_glob = PRF_GRSUM1(Flux*Amagnitude);
		
		A_glob = PI*0.25*pow(Dia,2); /*PRF_GRSUM1(Amagnitude);*/

		mflux_glob = mdot_glob / A_glob; 

		PA_glob = PRF_GRSUM1(pressure*Amagnitude);

		velocity_glob = mflux_glob / (density*A_glob);

		reynolds = density * fabs(velocity_glob) * Dia / mu;  /*calulateReynolds(domain_downstream,ID); */

		friction_factor = 1.0/pow((1.8*log10(reynolds)-1.5),2);

		press_avg = (PA_glob/A_glob) + (0.5*density*pow(velocity_glob,2)) - (friction_factor*(L/Dia)*(0.5*density*pow(velocity_glob,2)))*frelax;
	}
	else
	{
		press_avg =  0.0;
	}
		begin_f_loop(f, thread)
		{			
	        /*velocity = mflux_glob/ (density*(PI*0.25*pow(Dia,2)));*/
			F_PROFILE(f, thread, position) = press_avg;
		}
		end_f_loop(f, thread)	
	Message0("\n");
	/*Message0("area of %d is %g\n", zone_ID, A_glob);  Area is correct */
	Message0("massflowrate at in%d: is %g\n", zone_ID, mflux_glob);
	Message0("reynolds at %d is %g\n", zone_ID, reynolds);
	Message0("ff at at %d is %g\n", zone_ID, friction_factor);
	Message0("velocity %d is %g\n", zone_ID, velocity_glob);
	Message0("calculated_p at %d is %g\n", zone_ID, press_avg);
	Message0("\n");
}
