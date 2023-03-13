#ifndef SPHEREH
#define SPHEREH

#include "hitable.h"

class sphere : public hitable {
public:
	sphere() {}
	sphere(vec3 cen, float r, material* m) : center(cen), radius(r), mtr_ptr(m) {};

	virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec) const;
	virtual bool hit(const vec3& p1, const vec3& p2, float t_max, hit_record& rec) const;

	vec3 center;
	float radius = 1;
	material* mtr_ptr = nullptr;
};

bool sphere::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
	vec3 oc = r.origin() - center;
	float a = dot(r.direction(), r.direction());
	float b = dot(oc, r.direction());
	float c = dot(oc, oc) - radius * radius;
	float discriminant = b * b - a * c;


	if (discriminant > 0) {
		float temp = (-b - sqrt(b * b - a * c)) / a;
		if (temp < t_max && temp > t_min) {
			rec.t = temp;
			rec.p = r.point_at_parameter(rec.t);
			rec.normal = (rec.p - center) / radius;
			rec.mat_ptr = mtr_ptr;
			return true;
		}
		temp = (-b + sqrt(b * b - a * c)) / a;
		if (temp< t_max && temp > t_min) {
			rec.t = temp;
			rec.p = r.point_at_parameter(rec.t);
			rec.normal = (rec.p - center) / radius;
			rec.mat_ptr = mtr_ptr;
			return true;
		}
	}
	return false;
}

bool sphere::hit(const vec3& p1, const vec3& p2, float t_max, hit_record& rec) const {
	float t_min = 1e-3;
	vec3 oc = p1 - center;
	float tol = 1e-6;
	float a = dot(p2 - p1, p2 - p1);
	float b = dot(p2 - p1, oc);
	float c = dot(oc, oc) - radius * radius;
	float discriminant = b * b - a * c;
	if (discriminant > 0) {
		float temp = (-b - sqrt(b * b - a * c)) / a;
		if ((temp <= t_max) && (temp > t_min)) {
			rec.t = temp;
			rec.p = (p2 - p1) * temp + p1;
			rec.normal = (rec.p - center) / radius;
			rec.mat_ptr = mtr_ptr;
			return true;
		}
		//I'm not really sure why this soluition would ever be pickked
		temp = (-b + sqrt(b * b - a * c)) / a;
		if ((temp < t_max) && (temp > t_min)) {
			rec.t = temp;
			rec.p = (p2 - p1) * temp + p1;
			rec.normal = (rec.p - center) / radius;
			rec.mat_ptr = mtr_ptr;
			return true;
		}
		
	}
	return false;
}






#endif
