#include <iostream>
#include <random>
#include <vector>
#include <future>
#include <thread>
#include <atomic>
#include <math.h>
#include <iomanip>
#include "stb_image_write.h"
#include "vec3.h"
#include "ray.h"
#include "hitablelist.h"
#include "sphere.h"
#include "camera.h"
#include "material.h"


ray integration_step(const ray& r, float h);
//ray integration_step_rk_adapt(const ray& r);
vec3 f(const ray& r);
vec3 f(vec3 a, vec3 b);
vec3 f(vec3 org);
double* potential(double u[]);
double* rk4vec(double t0, double u0[], double dt, double* func(double u[]));
//This is the function that will be changed based on new coordinates.


/*

vec3 color(const ray& r, hitable *world, int depth) {
	hit_record rec;
	if (world->hit(r, 0.001f, FLT_MAX, rec)){
		ray scattered;
		vec3 attenuation;
		if ((depth < 50) && (rec.mat_ptr->scatter(r, rec, attenuation, scattered))) {
			return attenuation * color(scattered, world, depth + 1);
		}
		else {
			return vec3(0.0f, 0.0f, 0.0f);
		}
	}
	
	//Doesn't hit hitable in the world
	vec3 unit_direction = unit_vector(r.direction());
	float t = .5f * (unit_direction.y() + 1.0f);
	return (1.0f - t) * vec3(1.0f, 1.0f, 1.0f) + t * vec3(.5f, .7f, 1.0f);
}
*/
//This basically gets the color of 1 ray each step.
//Need to change depth, to stepsize.

vec3 color(const ray& r, hitable* world, int depth) {
	hit_record rec;
	ray r_c = r;
	const int MAX_INT_STEPS = 10000;
	float h = .1;
	float v = 1;
	int cur_steps = 0;
	while (cur_steps < MAX_INT_STEPS && r_c.origin().length() <= 20.0f) {
		cur_steps++;



		




		//Do step of the light ray
		ray r_n = integration_step(r_c, h);
		//check if there was an object that was hit
		//TODO: Need a way to get the distance between these two vectors as the ray of reflection.
		// These two rays are assumed to be so close to each other that they are basically linear.

		if ((vec3(0.0f, 0.0f, 1.0f) - r_n.origin()).length() < .5)
			return vec3(0.0f, 0.0f, 0.0f);


		if (world->hit(r_c.origin(), r_n.origin(), FLT_MAX, rec)) {
			ray scattered;
			vec3 attenuation;
			if ((depth < 50) && (rec.mat_ptr->scatter(r_c, rec, attenuation, scattered))) {
				return attenuation * color(scattered, world, depth + 1);
			}
			else {
				return vec3(0.0f, 0.0f, 0.0f);
			}
		}
		else {
			r_c = r_n;
		}

	}
	//std::cout << "Reaching Integration Limit" << r_c.origin().length() << " Max Depthd: " << cur_steps << std::endl;
	vec3 unit_direction = unit_vector(r.direction());
	float t = .5f * (unit_direction.y() + 1.0f);
	return (1.0f - t) * vec3(1.0f, 1.0f, 1.0f) + t * vec3(.5f, .7f, 1.0f);
}

/*
ray integration_step(const ray& r, float h) {
	//Euler First Order Integration
	vec3 x_p1 = h * r.direction() + r.origin();
	vec3 v_p1 = h * f(r) + r.direction();

	v_p1.make_unit_vector();
	return ray(x_p1, v_p1);
}


ray integration_step(const ray& r, float h = 1e-3) {
	const float TOL = 1e-5;
	float err = 2 * TOL;
	vec3 new_position, yn;

	while (err > TOL) {

		vec3 k1 = h * f(r);
		vec3 k2 = h * f(r.origin() + (1.0f / 5.0f) * h * vec3(1.0f, 1.0f, 1.0f), r.direction() + (1.0f / 5.0f) * k1);
		vec3 k3 = h * f(r.origin() + (3.0f / 10.0f) * h * vec3(1.0f, 1.0f, 1.0f), r.direction() + ((3.0f / 40.0f) * k1) + ((9.0f / 40.0f) * k2));
		vec3 k4 = h * f(r.origin() + (3.0f / 5.0f) * h * vec3(1.0f, 1.0f, 1.0f), r.direction() + ((3.0f / 10.0f) * k1) - ((9.0f / 10.0f) * k2) + ((6.0f / 5.0f) * k3));
		vec3 k5 = h * f(r.origin() + (1.0f / 1.0f) * h * vec3(1.0f, 1.0f, 1.0f), r.direction() - ((11.0f / 54.0f) * k1) + ((5.0f / 2.0f) * k2) - ((70.0f / 27.0f) * k3) + ((35.0f / 27.0f) * k4));
		vec3 k6 = h * f(r.origin() + (7.0f / 8) * h * vec3(1.0f, 1.0f, 1.0f), r.direction() + ((1631.0f / 55296) * k1) + ((175.0f / 512) * k2) + ((575.0f / 13824) * k3) + ((44275.0f / 110592) * k4) + ((253.0f / 4096) * k5));
		vec3 dy4 = ((37.0f / 378) * k1) + ((250.0f / 621) * k3) + ((125.0f / 594) * k4) + ((512.0f / 1771) * k6);
		vec3 dy5 = ((2825.0f / 27648) * k1) + ((18575.0f / 48384) * k3) + ((13525.0f / 55296) * k4) + ((277.0f / 14336) * k5) + ((1.0f / 4) * k6);
		err = 1e-2 * TOL + (dy4 - dy5).length();
		h = 0.8 * h * pow((TOL * h / err), (1.0f / 4));
		yn = r.direction() + dy4;
		
		new_position = yn * h * r.origin();

	}
	yn.make_unit_vector();
	return ray(new_position, yn);

}
*/

ray integration_step(const ray& r, float h = 1) {
	double* u = new double[6];
	u[0] = r.origin()[0];
	u[1] = r.origin()[1];
	u[2] = r.origin()[2];
	u[3] = r.direction()[0];
	u[4] = r.direction()[1];
	u[5] = r.direction()[2];

	double* ans;
	ans = rk4vec(0, u, h, potential);

	vec3 vel(ans[3], ans[4], ans[5]);

	vel.make_unit_vector();

	ray ret_ray(vec3(ans[0], ans[1], ans[2]), vel);

	delete[] ans;
	delete[] u;
	return ret_ray;
}


double* rk4vec(double t0, double u0[], double dt, double* func(double u[]))
{
	double* f0;
	double* f1;
	double* f2;
	double* f3;
	int i;
	double t1;
	double t2;
	double t3;
	double* u;
	double* u1;
	double* u2;
	double* u3;
	//
	//  Get four sample values of the derivative.
	//
	f0 = func(u0);

	t1 = t0 + dt / 2.0;
	u1 = new double[6];
	for (i = 0; i < 6; i++)
	{
		u1[i] = u0[i] + dt * f0[i] / 2.0;
	}
	f1 = func(u1);

	t2 = t0 + dt / 2.0;
	u2 = new double[6];
	for (i = 0; i < 6; i++)
	{
		u2[i] = u0[i] + dt * f1[i] / 2.0;
	}
	f2 = func(u2);

	t3 = t0 + dt;
	u3 = new double[6];
	for (i = 0; i < 6; i++)
	{
		u3[i] = u0[i] + dt * f2[i];
	}
	f3 = func(u3);
	//
	//  Combine them to estimate the solution.
	//
	u = new double[6];
	for (i = 0; i < 6; i++)
	{
		u[i] = u0[i] + dt * (f0[i] + 2.0 * f1[i] + 2.0 * f2[i] + f3[i]) / 6.0;
	}
	//
	//  Free memory.
	//
	delete[] f0;
	delete[] f1;
	delete[] f2;
	delete[] f3;
	delete[] u1;
	delete[] u2;
	delete[] u3;

	return u;
}
// Zero potential means that there will be no alteration to the final image

vec3 f(const ray& r) {
	return vec3(0, 0, 0);
}

vec3 f(vec3 org, vec3 dir) {
	return vec3(0, 0, 0);
}

vec3 f(vec3 org) {

	return vec3(0, 0, 0);
}

/*
double* potential(double u[]) {
	double* ans = new double[6];

	// du/dt = v(t)
	ans[0] = u[3];
	ans[1] = u[4];
	ans[2] = u[5];
	// dv/dt = f(u)

	ans[3] = 0;
	ans[4] = 0;
	ans[5] = 0;
	return ans;
}
*/
double* potential(double u[]) {

	double* ans = new double[6];

	// du/dt = v(t)
	ans[0] = u[3];
	ans[1] = u[4];
	ans[2] = u[5];
	const double GM = 0.25;
	const double c = 1;
	//dv/dt = f(u)
	
	double x, y, z, dx, dy, dz, r, theta, phi, dr, dtheta, dphi, d2rdt2;
	//Setting up local variables
	x = u[0] - 0;
	z = u[1] - 0;
	y = u[2] + 1;
	dx = u[3];
	dz = u[4];
	dy = u[5];

	r = pow(x * x + y * y + z * z, .5);
	theta = atan(y / (x + 1e-9)); // Protection against x ~ 0
	phi = acos(z / (r + 1e-9)); // Protection against r ~ 0
	dr = (.5 * pow(x * x + y * y + z * z, -.5)) * (2 * x * dx + 2 * y * dy + 2 * z * dz);
	dtheta = (x * dy - y * dx) / (x * x + y * y + 1e-9);
	dphi = (z * dr - r * dz) / (r * r * pow(1 - (z * z) / (r * r), .5) + 1e-9);

	d2rdt2 = -GM * ((1 - (2 * GM) / (r * c * c)) / (1 - (3 * GM) / (r * c * c))) * (r / (r * r * r) - 1) / (r * r);



	ans[3] = d2rdt2 * cos(theta) * sin(phi) - 2 * dr * (sin(theta) * sin(phi) * dtheta - cos(theta) * cos(phi) * dphi) - 2 * r * sin(theta) * cos(phi) * dtheta * dphi;
	ans[5] = d2rdt2 * sin(theta) * sin(phi) + 2 * dr * (cos(theta) * sin(phi) * dtheta + sin(theta) * cos(phi) * dphi) + 2 * r * cos(theta) * cos(phi) * dtheta * dphi;
	ans[4] = d2rdt2 * cos(phi) - 2 * dr * dphi * sin(phi);

	//std::cout << std::fixed << std::setprecision(3) << "Pos: " << x << " " << y << " " << z <<  " " << r << " " << phi << " " << theta << std::endl;


	return ans;
}




//Vector potential function
/*
vec3 f(const ray& r) {
	float rd = (r.origin() - vec3(0.0f, 0.0f, -1.0f)).length();
	return -(r.origin() - vec3(0.0f, 0.0f, -1.0f)) / rd * (2 * rd * rd * rd - 3 * rd * rd + 1);
}

vec3 f(vec3 org, vec3 dir) {
	float rd = (org - vec3(0.0f, 0.0f, -1.0f)).length();
	return -(org - vec3(0.0f, 0.0f, -1.0f)) / rd * (2 * rd * rd * rd - 3 * rd * rd + 1);
}

*/




int main() {
	
	int numX = 400;
	int numY = 200;
	int numSamples = 50;
	int channels = 3;
	
	float aspectRatio = float(numX) / float(numY);
	int stride = numX*channels;
	size_t size = (size_t)numX * numY * 3;
	hitable *list[3];
	list[0] = new sphere(vec3(0.0f, 0.0f, -1.0f), .5f, new lambertian(vec3(0.8f, 0.3f, 0.3f)));
	list[1] = new sphere(vec3(0.0f, -100.5f, -1.0f), 100.0f, new lambertian(vec3(0.8f, 0.8f, 0.0f)));
	list[2] = new sphere(vec3(1.0f, 0.0f, -1.0f), 0.5f, new metal(vec3(0.8f, 0.6f, 0.2f), .5f));
	//list[3] = new sphere(vec3(-1.0f, 0.0f, -1.0f), .5f, new dielectric(1.5));

	hitable* world = new hitable_list(list, 3);
	camera cam;

	char* data = new char[size];

	//New thread logic buidling


	std::size_t max = numX * numY*channels;
	std::size_t cores = std::thread::hardware_concurrency();

	volatile std::atomic<std::size_t> count(0);
	std::vector<std::future<void>> future_vector;

	
	while (cores--) {
		future_vector.emplace_back(
			std::async([=, &world, &count, &cam](){
				while (true) {
					std::size_t index = count++;
					index = index * 3;
					if (index >= max)
						break;
					std::size_t x = (index/3) % numX;
					std::size_t y = (index/3) / numX;

					vec3 col(0.0f, 0.0f, 0.0f);

					for (int s = 0; s < numSamples; s++) {
						float u = float(x + get_uniform_rand()) / float(numX);
						float v = float(y + get_uniform_rand()) / float(numY);
						ray r = cam.get_ray(u, v);
						//vec3 p = r.point_at_parameter(2.0f);
						col += color(r, world, 0);
					}
					col /= float(numSamples);

					col = vec3(sqrt(col[0]), sqrt(col[1]), sqrt(col[2]));
					char red = (char)(col[0] * 255.99f);
					char green = (char)(col[1] * 255.99f);
					char blue = (char)(col[2] * 255.99f);

					data[channels * x + (numY - y - 1) * numX * channels] = red;
					data[channels * x + (numY - y - 1) * numX * channels + 1] = green;
					data[channels * x + (numY - y - 1) * numX * channels + 2] = blue;
					std::cout << index/3 << std::endl;
				}
			}	
		));
	}
















	/*
	for (int y = 0; y < numY; y++) {
		for (int x = 0; x < numX; x++) {

			vec3 col(0.0f, 0.0f, 0.0f);

			for (int s = 0; s < numSamples; s++) {
				float u = float(x + get_uniform_rand()) / float(numX);
				float v = float(y + get_uniform_rand()) / float(numY);


				ray r = cam.get_ray(u, v);
				//vec3 p = r.point_at_parameter(2.0f);
				col += color(r, world, 0);


			}
			col /= float(numSamples);

			col = vec3(sqrt(col[0]), sqrt(col[1]), sqrt(col[2]));
			char red = (char)(col[0] * 255.99f);
			char green = (char)(col[1] * 255.99f);
			char blue = (char)(col[2] * 255.99f);

			data[channels*x + (numY-y-1) * numX * channels] = red;
			data[channels*x + (numY-y-1) * numX * channels + 1] = green;
			data[channels*x + (numY-y-1) * numX * channels + 2] = blue;
			std::cout << int(float(y * numX + x) / (numY * numX) * 100) << std::endl;
		}
	}

	
	*/

	//Wait for all threads to finish their tasks before reading data. 
	for (int i = 0; i < future_vector.size(); i++)
		future_vector[i].get();

	int result = stbi_write_png("runge_kutta_black_hole.png", numX, numY, channels, data, stride);





	//Wrap up the data.
	delete[] data;



}