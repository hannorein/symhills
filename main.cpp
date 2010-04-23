/*******************************************************************************
 * Symplectic integrator for Hill's equartions
 * based on http://arxiv.org/abs/0908.2269
 *
 * Hanno Rein 2010
 * University of Cambridge
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see http://www.gnu.org/licenses/
 *******************************************************************************/


#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <fstream>
#include <string>
using namespace std;

#define G 	6.6743e-11
#define OMEGA 	0.00013143546
#define A	130000000.0
#define dt	1.0e0

void kick1(double* r, double* v, double* a, double* Py){
	v[0] += -0.5 *dt * (OMEGA*OMEGA*r[0] - a[0]);
	(*Py)  = v[1]     + 2.0*OMEGA * r[0]     + 0.5 * dt * a[1];
	v[0] += dt*OMEGA*(*Py);
	v[1]  = (*Py)     - OMEGA * r[0]     - OMEGA * (r[0] + dt * v[0]);
}

void drift(double* r, double* v){
	r[0] += dt * v[0];
	r[1] += dt * v[1];
}

void accel(double* r1, double* r2, double* a, double mas){
	double _r = sqrt(
		 (r1[0]-r2[0])*(r1[0]-r2[0])
		+(r1[1]-r2[1])*(r1[1]-r2[1])
	);
	a[0] = -G*mas/_r/_r/_r*(r1[0]-r2[0]);
	a[1] = -G*mas/_r/_r/_r*(r1[1]-r2[1]);
}

void kick2(double* r, double* v, double* a, double* Py){
	
	v[0] += dt * OMEGA * (*Py);
	v[0] += - 0.5 * dt * (OMEGA * OMEGA * r[0] - a[0]);
	v[1]  = (*Py) - 2.0 * OMEGA * r[0] + 0.5 * dt * a[1];
}

int main (int argc, char * argv[]) {
	double p_r[2];
	double p_v[2],p_Py;
	double p_a[2];

	// Setting up test problem
	// Particle on eccentric orbit
	p_r[0]=1000.;
	p_r[1]=100.;
	p_v[0]=1.0;	// set to zero for circular orbit
	p_v[1]=-1.5*OMEGA*p_r[0];
	p_a[0]=0.;
	p_a[1]=0.;
	p_Py  = p_v[1]     + 2.0*OMEGA * p_r[0];

	// Integrate for 10 orbits
	double t_max = 10.*2.0*M_PI/OMEGA;
	while (t_max>0){
		kick1(p_r,p_v,p_a,&p_Py);

		drift(p_r,p_v);

		//Acceleration forces
		//accel(p_r,m_r,p_a,m_m);

		kick2(p_r,p_v,p_a,&p_Py);

		t_max -= dt;
		// Output positions
		cout << p_r[0] <<  "\t" << p_r[1] << endl;
	}
			
	return 0;
}

