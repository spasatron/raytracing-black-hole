#include <iostream>
#include <random>
#include <vector>
#include <future>
#include <thread>
#include <atomic>
#include "stb_image_write.h"
#include "vec3.h"
#include "ray.h"
#include "hitablelist.h"
#include "sphere.h"
#include "camera.h"
#include "material.h"


ray integration_step(const ray& r, float h);
vec3 f(const ray& r);

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
	const int MAX_INT_STEPS = 1000;
	float h = .0001;
	float v = 1;
	int cur_steps = 0;
	while (cur_steps < MAX_INT_STEPS) {
		cur_steps++;
		//Do step of the light ray
		r_c = integration_step(r_c, h);
		//check if there was an object that was hit
		if (world->hit(r, 0.00001f, .0001f, rec)) {
			ray scattered;
			vec3 attenuation;
			if ((depth < 50) && (rec.mat_ptr->scatter(r, rec, attenuation, scattered))) {
				return attenuation * color(scattered, world, depth + 1);
			}
			else {
				return vec3(0.0f, 0.0f, 0.0f);
			}
		}

	}
	vec3 unit_direction = unit_vector(r.direction());
	float t = .5f * (unit_direction.y() + 1.0f);
	return (1.0f - t) * vec3(1.0f, 1.0f, 1.0f) + t * vec3(.5f, .7f, 1.0f);
}

ray integration_step(const ray& r, float h) {
	//Euler First Order Integration
	vec3 new_position = h * r.direction() + r.origin();
	vec3 new_velocity = (h * f(r) + r.direction());
	new_velocity.make_unit_vector();
	return ray(new_position, new_velocity);
}
//Vector potential function
vec3 f(const ray& r) {
	float rd = (r.origin() - vec3(0.0f, 0.0f, -1.0f)).length();
	return -(r.origin() - vec3(0.0f, 0.0f, -1.0f)) / rd * (2 * rd * rd * rd - 3 * rd * rd + 1);
}








int main() {
	
	int numX = 300;
	int numY = 150;
	int numSamples = 1000;
	int channels = 3;
	
	float aspectRatio = float(numX) / float(numY);
	int stride = numX*channels;
	size_t size = (size_t)numX * numY * 3;
	hitable *list[4];
	list[0] = new sphere(vec3(0.0f, 0.0f, -1.0f), .5f, new lambertian(vec3(0.8f, 0.3f, 0.3f)));
	list[1] = new sphere(vec3(0.0f, -100.5f, -1.0f), 100.0f, new lambertian(vec3(0.8f, 0.8f, 0.0f)));
	list[2] = new sphere(vec3(1.0f, 0.0f, -1.0f), 0.5f, new metal(vec3(0.8f, 0.6f, 0.2f), .5f));
	list[3] = new sphere(vec3(-1.0f, 0.0f, -1.0f), .5f, new dielectric(1.5));

	hitable* world = new hitable_list(list, 4);
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

					data[index] = red;
					data[index + 1] = green;
					data[index + 2] = blue;
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

	int result = stbi_write_png("mulitthreading.png", numX, numY, channels, data, stride);





	//Wrap up the data.
	delete[] data;



}