#pragma once
#include <math.h>
#include <stdlib.h>
#include <iostream>



class vec3 {
public:
	vec3() {}
	vec3(float e0, float e1, float e2) {
		e[0] = e0;
		e[1] = e1;
		e[2] = e2;
	}
	//Position names
	inline float x() const { return e[0]; }
	inline float y() const { return e[1]; }
	inline float z() const { return e[2]; }
	//Color Names
	inline float r() const { return e[0]; }
	inline float g() const { return e[1]; }
	inline float b() const { return e[2]; }

	//Math Definitions
	inline const vec3& operator+() const { return *this; }
	inline vec3 operator-() const { return vec3(-e[0], -e[1], -e[2]); }
	inline float operator[](int i) const { return e[i]; }
	inline float& operator[](int i) { return e[i]; }

	//Itterative Maths 
	inline vec3& operator+=(const vec3& v2);
	inline vec3& operator-=(const vec3& v2);
	inline vec3& operator*=(const vec3& v2);
	inline vec3& operator/=(const vec3& v2);
	inline vec3& operator*=(const float t);
	inline vec3& operator/=(const float t);

	//Function definitions
	inline float length() const {
		return sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
	}
	inline float squared_length() const {
		return e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
	}
	inline void make_unit_vector();
	inline vec3& unit_vec();
	float e[3] = { 0.0f , 0.0f, 0.0f };
};



	inline std::istream& operator>>(std::istream& is, vec3& t) {
		is >> t.e[0] >> t.e[1] >> t.e[2];
		return is;
	}

	inline std::ostream& operator<<(std::ostream& os, vec3& t) {
		os << t.e[0] << t.e[1] << t.e[2];
		return os;
	}

	inline vec3 operator+(const vec3& v1, const vec3& v2) {
		return vec3(v1.e[0] + v2.e[0], v1.e[1] + v2.e[1], v1.e[2] + v2.e[2]);
	}

	inline vec3 operator-(const vec3& v1, const vec3& v2) {
		return vec3(v1.e[0] - v2.e[0], v1.e[1] - v2.e[1], v1.e[2] - v2.e[2]);
	}

	inline vec3 operator*(const vec3& v1, const vec3& v2) {
		return vec3(v1.e[0] * v2.e[0], v1.e[1] * v2.e[1], v1.e[2] * v2.e[2]);
	}

	inline vec3 operator/(const vec3& v1, const vec3& v2) {
		return vec3(v1.e[0] / v2.e[0], v1.e[1] / v2.e[1], v1.e[2] / v2.e[2]);
	}

	inline vec3 operator*(float t, const vec3& v1) {
		return vec3(t * v1.e[0], t * v1.e[1], t * v1.e[2]);
	}

	inline vec3 operator*(const vec3& v1, float t) {
		return vec3(t * v1.e[0], t * v1.e[1], t * v1.e[2]);
	}

	inline vec3 operator/(const vec3& v1, float t) {
		return vec3( v1.e[0]/t,  v1.e[1]/t, v1.e[2]/t);
	}


	inline float dot(const vec3& v1, const vec3& v2) {
		return v1.e[0] * v2.e[0] + v1.e[1] * v2.e[1] + v1.e[2] * v2.e[2];
		return v1.e[0] * v2.e[0] + v1.e[1] * v2.e[1] + v1.e[2] * v2.e[2];
	}
	inline vec3 cross(const vec3& v1, const vec3& v2) {
		return vec3((v1.e[1] * v2.e[2] - v1.e[2] * v2.e[1]),
			(-(v1.e[0] * v2.e[2] - v1.e[2] * v2.e[0])),
			(v1.e[0] * v2.e[1] - v1.e[1] * v2.e[0]));
	}

	inline vec3 unit_vector(vec3 v) {
		return v / v.length();
	}

	inline vec3& vec3::operator+=(const vec3& v) {
		e[0] += v.e[0];
		e[1] += v.e[1];
		e[2] += v.e[2];
		return *this;
	}

	inline vec3& vec3::operator-=(const vec3& v) {
		e[0] -= v.e[0];
		e[1] -= v.e[1];
		e[2] -= v.e[2];
		return *this;
	}

	inline vec3& vec3::operator*=(const vec3& v) {
		e[0] *= v.e[0];
		e[1] *= v.e[1];
		e[2] *= v.e[2];
		return *this;
	}

	inline vec3& vec3::operator/=(const vec3& v) {
		e[0] /= v.e[0];
		e[1] /= v.e[1];
		e[2] /= v.e[2];
		return *this;
	}

	inline vec3& vec3::operator*=(const float t) {
		e[0] *= t;
		e[1] *= t;
		e[2] *= t;
		return *this;
	}

	inline vec3& vec3::operator/=(const float t) {
		float k = 1.0f / t;
		e[0] *= k;
		e[1] *= k;
		e[2] *= k;
		return *this;
	}

	inline void vec3::make_unit_vector() {
		float k = 1.0f / sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
		e[0] *= k;
		e[1] *= k;
		e[2] *= k;
	}

	inline vec3& vec3::unit_vec() {
		float k = 1.0f / sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
		e[0] *= k;
		e[1] *= k;
		e[2] *= k;
		return *this;
	}

	float get_uniform_rand() {
		return (float(rand()) / (float(RAND_MAX) + 1.0f));
	}


	vec3 random_in_unit_sphere() {
		vec3 p;
		do {
			p = 2.0 * vec3(get_uniform_rand(), get_uniform_rand(), get_uniform_rand()) - vec3(1.0f, 1.0f, 1.0f);
		} while (p.squared_length() >= 1.0);
		return p;
	}

	vec3 reflect(const vec3& v, const vec3& n) {
		return v - 2 * dot(v, n) * n;
	}


	bool refract(const vec3& v, const vec3& n, float ni_over_nt, vec3& refracted) {
		vec3 uv = unit_vector(v);
		float dt = dot(uv, n);
		float disciminant = 1.0 - ni_over_nt * ni_over_nt * (1 - dt * dt);
		if (disciminant > 0) {
			refracted = ni_over_nt * (uv - n * dt) - n * sqrt(disciminant);
			return true;
		}
		return false;
	}


