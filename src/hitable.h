#ifndef HITABLEH
#define HITABLEH
#include "ray.h"


class material;

struct hit_record{
	float t;
	vec3 p;
	vec3 normal;
	material* mat_ptr;
};



class hitable {
public:
	virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec) const = 0;
	virtual bool hit(const vec3& p1, const vec3& p2, float t_max, hit_record& rec) const = 0;
};



#endif // !HITABLEH

