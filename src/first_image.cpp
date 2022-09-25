#include <iostream>
#include <random>
#include "stb_image_write.h"
#include "vec3.h"
#include "ray.h"
#include "hitablelist.h"
#include "sphere.h"
#include "camera.h"
#include "material.h"





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










int main() {
	
	int numX = 1200;
	int numY = 600;
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

	for (int y = 0; y < numY; y++) {
		for (int x = 0; x < numX; x++) {

			vec3 col(0.0f, 0.0f, 0.0f);

			for (int s = 0; s < numSamples; s++) {
				float u = float(x + get_uniform_rand()) / float(numX);
				float v = float(y + get_uniform_rand()) / float(numY);


				ray r = cam.get_ray(u, v);
				vec3 p = r.point_at_parameter(2.0f);
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
		}
	}


	int result = stbi_write_png("output.png", numX, numY, channels, data, stride);





	//Wrap up the data.
	delete[] data;



}