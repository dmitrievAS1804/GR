#include "Geometry.h"


bool Plane::Intersect(const Ray &ray, float t_min, float t_max, SurfHit &surf) const
{
  surf.t = dot((point - ray.o), normal) / dot(ray.d, normal);

  if(surf.t > t_min && surf.t < t_max)
  {
    surf.hit      = true;
    surf.hitPoint = ray.o + surf.t * ray.d;
    surf.normal   = normal;
    surf.m_ptr    = m_ptr;
    return true;
  }

  return false;
}

bool Sphere::Intersect(const Ray& ray, float t_min, float t_max, SurfHit& surf) const 
{
    float3 k = ray.o - center;
    float a = dot(ray.d, ray.d);
    float b = dot(2 * k, ray.d);
    float c = dot(k, k) - r_sq;
    float det = b * b - 4 * a * c;
    if (det < 0) return false;
        surf.t = (-b - sqrt(det)) / 2 * a;
        if (surf.t > t_min && surf.t < t_max)
        {
            surf.hit = true;
            surf.hitPoint = ray.o + surf.t * ray.d;
            surf.normal = normalize(surf.hitPoint - center);
            surf.m_ptr = m_ptr;
            return true;
        }
        surf.t = (-b + sqrt(det)) / 2 * a;
        if (surf.t > t_min && surf.t < t_max)
        {
            surf.hit = true;
            surf.hitPoint = ray.o + surf.t * ray.d;
            surf.normal = normalize(surf.hitPoint - center);
            surf.m_ptr = m_ptr;
            return true;
        }
    return false;
}

bool Triangle::Intersect(const Ray& ray, float t_min, float t_max, SurfHit& surf) const
{

    surf.normal = cross((B - A), (C - A));
    if (dot(ray.d, surf.normal) == 0) return false;
    float3 t = ray.o - A;
    float3 e1 = B - A;
    float3 e2 = C - A;
    float det = dot(e1, cross(ray.d, e2));
    float inv_det = 1 / det;
    if (det<t_min && det>t_max) return false;
    float u = dot(t, cross(ray.d, e2))* inv_det;
    float v = dot(ray.d, cross(t, e1))* inv_det;
    if (u < 0 || u>1) return false;
    if (v < 0 || u + v>1) return false;
    surf.t = dot(e2, cross(t, e1)) * inv_det;
    if (surf.t > t_min && surf.t < t_max) {
        surf.hit = true;
        surf.hitPoint = float3(surf.t, u, v);
        surf.m_ptr = m_ptr;
        return true;
    }

    return false;
}


/////////////////////////////////////////
