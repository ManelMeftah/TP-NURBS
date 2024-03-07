#ifndef VEC3_H
#define VEC3_H

#include <cmath>
#include <iostream>

using std::sqrt;

class vec3 {
    public:
            double e[3];
    public:
        vec3() : e{0,0,0} {} //what does the : mean ?? is it return type?
        vec3(double e0, double e1, double e2) : e{e0,e1,e2} {}

        double x() const { return e[0]; }
        double y() const { return e[1]; }
        double z() const { return e[2]; }

        //-V
        vec3 operator-() const { return vec3(-e[0], -e[1], -e[2]); }
        // V[i]
        double operator[](int i) const { return e[i]; }
        // &V[i] ? TODO try this
        double& operator[](int i) { return e[i]; }

        // V1 += V2
        vec3& operator+=(const vec3 &v)
        {
            e[0] += v.e[0];
            e[1] += v.e[1];
            e[2] += v.e[2];
            return *this;
        }
        
        // V1 *= n
        vec3& operator*=(const double n)
        {
            e[0] *=n;
            e[1] *= n;
            e[2] *= n;
            return *this;
        }

        // V1 /= n
        vec3& operator/=(const double n)
        {
            return *this *= 1/n;
        }

        double length() const
        {
            return sqrt(length_squared());
        }

        double length_squared() const
        {
            return e[0]*e[0] + e[1]*e[1] + e[2]*e[2];
        }

};

// point coordinates and color vectors are type vec3
using point3 = vec3;
using color = vec3;

// utilities and vector operations

inline std::ostream& operator<<(std::ostream &out, const vec3 &v)
{
    return out << "V[" << v.e[0] << ", " <<  v.e[1] << ", " <<  v.e[2] << "]\n";
}

inline vec3 operator+(const vec3 &u, const vec3 &v)
{
    return vec3(u.e[0] + v.e[0], u.e[1] + v.e[1], u.e[2] + v.e[2]);
}

inline vec3 operator-(const vec3 &u, const vec3 &v)
{
    return vec3(u.e[0] - v.e[0], u.e[1] - v.e[1], u.e[2] - v.e[2]);
}

inline vec3 operator*(const vec3 &u, const vec3 &v)
{
    return vec3(u.e[0] * v.e[0], u.e[1] * v.e[1], u.e[2] * v.e[2]);
}

inline vec3 operator*( double n, const vec3 &u)
{
    return vec3(u.e[0] * n, u.e[1] * n, u.e[2] * n);
}

inline vec3 operator*(const vec3 &u, double n)
{
    return n*u;
}

inline vec3 operator/(const vec3 &u, double n)
{
    return u * (1/n);
}

// U.V
inline double dot(const vec3 &u, const vec3 &v)
{
    return u.e[0] * v.e[0]
         + u.e[1] * v.e[1]
         + u.e[2] * v.e[2];
}

// U x V
inline vec3 cross(const vec3 &u, const vec3 &v)
{
    return vec3(u.e[1] * v.e[2] - u.e[2] * v.e[1],
                u.e[2] * v.e[0] - u.e[0] * v.e[2],
                u.e[0] * v.e[1] - u.e[1] * v.e[0]);
}

inline vec3 vecteur_unitaire(vec3 v)
{
    return v/ v.length();
}


#endif
