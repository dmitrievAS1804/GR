#include <iostream>
#include <vector>
#include <random>

#include "LiteMath.h"
#include "Geometry.h"
#include "Camera.h"

using namespace HydraLiteMath;

void RenderScene(uint32_t w, uint32_t h, uint32_t num_samples, const std::vector<std::shared_ptr<GeoObject>> &scene, const Camera &cam, const std::string &filename)
{
  auto  background_color = float3(1.0f, 0.3f, 0.2f);
  auto  film = std::make_unique<Film>(w, h, num_samples);
  auto  tracer = std::make_unique<SimpleRT>(16, background_color);
  float invWidth  = 1.0f / float(w);
  float invHeight = 1.0f / float(h);

  for (int y = 0; y < h; ++y)
  {
    for (int x = 0; x < w; ++x)
    {
      float3 pixel_color = float3(0.0f, 0.0f, 0.0f);

      for (int s = 0; s < num_samples; ++s)
      {
        Ray ray = cam.genRay(w, h, x, h - y); //генерируем луч из камеры через текущий пиксель
        pixel_color += tracer->TraceRay(ray, scene, 0); //трассируем луч и получаем цвет
      }
      pixel_color /= film->num_samples;      // усредняем полученные цвета
      pixel_color *= cam.getExposureTime();  // умножаем на время экспонирования сенсора - выдержка виртуальной камеры
      film->SetPixelColor(x, y, pixel_color); //записываем цвет на виртуальный сенсор
    }
  }
  film->SaveImagePPM(filename); //сохраняем картинку
}

void create_scene()
{
  std::vector<std::shared_ptr<GeoObject>> myScene;
  float3        eye   (0.0f, 2.0f, 20.0f);
  float3        lookat(0.0f, 2.0f, 0.0f);
  float3        up    (0.0f, 1.0f, 0.0f);
  float         field_of_view = 90.0f;
  unsigned int  w = 1920;
  unsigned int  h =  1080;
  Camera        cam(eye, lookat, up, field_of_view, float(w) / float(h));

  auto  pl= std::make_shared<Plane>(float3(+0.0f, -1.0f, +0.0f), float3(0.0f, 1.0f, 0.0f), new diffuse(float3(0.5f, 0.5f, 0.5f)));
  myScene.push_back(pl);
  auto sphere = std::make_shared<Sphere>(float3(4.0f, 3.0f, 5.0f),2.0,  new IdealMirror(float3(0.50f, 1.00f, 0.30f)));
  myScene.push_back(sphere);
  auto triangle = std::make_shared<Triangle>(float3(2.0f,5.0f, 6.0f),float3(1.0f, 12.0f, 5.0f),float3(0.0f, 6.0f, 8.0f), new diffuse(float3(0.80f, 0.00f, 0.00f)));
  myScene.push_back(triangle);
  RenderScene(w, h, 1, myScene, cam,  "basic_scene");
}

int main()
{
  create_scene();

  return 0;
}
