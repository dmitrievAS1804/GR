#include "Material.h"

SurfHit& SurfHit::operator= (const SurfHit &rhs)
{
  if (this == &rhs)
    return (*this);

  hit = rhs.hit;
  hitPoint = rhs.hitPoint;
  normal = rhs.normal;
  m_ptr = rhs.m_ptr;
  t = rhs.t;

  return (*this);
}

//////////////////

float3 reflect(const float3& v, const float3& normal)
{
    return v - 2 * dot(v, normal) * normal;
}


bool IdealMirror::Scatter(const Ray &ray_in, const SurfHit &surf, float3 &attenuation, Ray &ray_out)
{
    float3 reflection_dir = reflect(normalize(ray_in.d), surf.normal); 

    ray_out = Ray(surf.hitPoint + surf.normal * 1e-4, reflection_dir); 

    attenuation = color;

    return (dot(ray_out.d, surf.normal) > 0);
}

bool diffuse::Scatter(const Ray& ray_in, const SurfHit& surf, float3& attenuation, Ray& ray_out)
{
    ray_out = Ray(surf.hitPoint + surf.normal * 1e-4, surf.normal);

    attenuation = color;

    return (dot(ray_out.d, surf.normal) > 0);
}
