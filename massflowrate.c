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

#define frelax 0.00011

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
	/*real Area[ND_ND];*/
	real  Amagnitude, Flux, mflux_glob;
	real mdot_sum, mdot_glob;
	real A_sum, A_glob, dens1, A;

	Amagnitude = 0.0;
	dens1 = 0.0;
	A_glob = 0.0;
	dens1 = 0.0;
	mdot_glob = 0.0;
	
	face_t f_upstream;
	Domain* domain_upstream;         
	domain_upstream = Get_Domain(1);  
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
		begin_f_loop(f_upstream, thread_upstream)
			if PRINCIPAL_FACE_P(f_upstream, thread_upstream)
			{
				F_AREA(area, f_upstream, thread_upstream);
				Amagnitude = NV_MAG(area);
				Flux = F_FLUX(f_upstream, thread_upstream);
				mdot_sum += Flux;
				A_sum += Amagnitude;
				dens1 += F_R(f_upstream, thread_upstream)*NV_MAG(area);
			}
		end_f_loop(f_upstream, thread_upstream)
		
	A_glob =PRF_GRSUM1(A_sum);
	dens1 = PRF_GRSUM1(dens1);
	dens1  = dens1 / A_glob ;
	mdot_glob = PRF_GRSUM1(mdot_sum); 	/* Equivalent to PRF_GRSUM1(Flux*A_glob) */		
	mflux_glob =  mdot_glob / A_glob;
	}
	else
	{
		mdot_glob = 0.0;
	}
	
	begin_f_loop(f, thread)
	{
		
		F_PROFILE(f, thread, position) = mdot_glob;	
	}
	end_f_loop(f, thread)
	Message0("\n----------------------------------------------\n");	
	Message0("massflowrate @ RHS id %d is %g\n", zone_ID, mdot_glob);	
}


/*DEFINE_PROFILE(pressure_outlet, thread, position)*/
/*{*/
/*	/* Averages static pressure on downstream interface (in:*) inflow BC and applies it at the upstream outflow interface  (out:*)*/
/*	 */
/*	int ID; /* the ID of zone with downstream boundary */
/*	face_t f;*/
/*	real NV_VEC(area);*/
/**/
/*	/*real Area[ND_ND];*/
/*	real Amagnitude, press_avg;*/
/*	real A_sum, A_glob, A;*/
/*	real mdot_sum, Flux, mflux_glob , reynolds, velocity, friction_factor;*/
/*	real press, dens, Density;*/
/**/
/*	mdot_sum=0.0;*/
/*	press = 0.0;*/
/*	A_sum = 0.0;	*/
/*	dens = 0.0;*/
/*	A_glob = 0.0;*/
/*	Flux = 0.0;*/
/*	press=0.0;*/
/*	Density = 0.0;*/
/**/
/*	face_t f_downstream;*/
/*	Domain* domain_downstream;         */
/*	domain_downstream = Get_Domain(1); */
	/**/
/*	/* set the ID of the downstream interface */
/*	int zone_ID = THREAD_ID(thread); /* applied to */
/*		if (zone_ID == out1)*/
/*	{*/
/*		ID = in1;*/
/*	}*/
/*	else if (zone_ID == out2)*/
/*	{*/
/*		ID = in2;*/
/*	}*/
/*	else if (zone_ID == out3)*/
/*	{*/
/*		ID = in3;*/
/*	}*/
/*	else if (zone_ID == out4)*/
/*	{*/
/*		ID = in4;*/
/*	}*/
/*	else*/
/*	{*/
/*		ID = 10000;*/
/*	}*/
/*	if (ID!=10000)*/
/*	{		*/
/*		Thread* thread_downstream = Lookup_Thread(domain_downstream, ID);*/
/*		Amagnitude = 0.0; */
/*		Flux = 0.0;*/
/*		press = 0.0;*/
/*		A_sum = 0.0;	*/
/*		dens = 0.0;*/
/*		begin_f_loop(f_downstream, thread_downstream)*/
/*			if PRINCIPAL_FACE_P(f_downstream, thread_downstream)*/
/*			{	*/
/*				F_AREA(area, f_downstream, thread_downstream); /*Outputs Area vector*/
/*				Amagnitude = NV_MAG(area);*/
/*				A_sum += Amagnitude;*/
/*				Flux = F_FLUX(f_downstream, thread_downstream);*/
/*				mdot_sum += Flux;*/
/*				/*temp = F_T(f_downstream,thread_downstream);*/
/*				cells = cells + 1.0;*/
/*				press +=  F_P(f_downstream, thread_downstream)*NV_MAG(area);*/
/*				dens += F_R(f_downstream, thread_downstream)*NV_MAG(area);*/
/*			}*/
/*		end_f_loop(f_downstream, thread_downstream)*/
/*		A = PI*0.25*Dia*Dia*1e-6;*/
/*		A_glob = PRF_GRSUM1(A_sum); /*2.7e-5*/
/*		Flux  = PRF_GRSUM1(mdot_sum);*/
/*		press = PRF_GRSUM1(press);*/
/*		Density = PRF_GRSUM1(dens);*/
/*		/*temp_sum = PRF_GRSUM1(temp);*/
/*		Density = Density / A_glob; /*Is correct!*/
/*		/*bulk_temp = temp_sum / cells;*/
/*		reynolds = 4*fabs(Flux)/(PI*1e-3*Dia*mu);	*/
/*		velocity = (reynolds*mu)/(Dia*1e-3*Density);*/
/*		/*velocity  = Flux /(Density * A_glob);*/
/*		friction_factor = 1/(pow((1.8*log10(reynolds)-1.5),2));*/
/*		press_avg = (press/A_glob)+ (0.5*Density*pow(velocity,2)) - (friction_factor*(L/Dia)*0.5*Density*pow(velocity,2));*/
/*	}*/
/*	else*/
/*	{*/
/*		press_avg =  0.0;*/
/*	}*/
		/**/
/*		begin_f_loop(f, thread)*/
/*		{			*/
/*			F_PROFILE(f, thread, position) = press_avg ; /*(force/A_glob)-(friction_factor*(L/Dia)*0.5*Density*pow(velocity,2))*frelax;*/
/*		}*/
/*		end_f_loop(f, thread)	*/
/**/
/*	Message0("\n---------------------------------------------\n");*/
/*	Message0("denisty is %g\n", Density);*/
/*	Message0("Area is %g\n", A_glob);*/
/*	Message0("Area manual is %g\n", A);*/
/*	/*Message0("Bulk temperature is %g \n", bulk_temp);*/
/*	Message0("massflowrate at out%d: is %g \n", zone_ID, Flux);*/
/*	Message0("reynolds at %d is %g \n", zone_ID, reynolds);*/
/*	Message0("velocity %d is %g \n", zone_ID, velocity);*/
/*	Message0("ff at at %d is %.8g\n", zone_ID, friction_factor);*/
/*	Message0("calculated_p at %d is %g \n", zone_ID, press_avg);*/
/*	Message0("\n");*/
/*} */