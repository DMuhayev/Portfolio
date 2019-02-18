#include "Tracer.h"
#include "iostream"
#include <vector>

using namespace glm;

SRay CTracer::MakeRay(glm::uvec2 pixelPos)
{
  SRay tmp;

  tmp.m_start = m_camera.m_pos;

  double X, Y, Z, R, U;
  R = (((pixelPos.x+0.5)/m_camera.m_resolution.x)-0.5);
  U = (((pixelPos.y+0.5)/m_camera.m_resolution.y)-0.5);

  X = m_camera.m_forward.x + R*m_camera.m_right.x + U*m_camera.m_up.x;
  Y = m_camera.m_forward.y + R*m_camera.m_right.y + U*m_camera.m_up.y;
  Z = m_camera.m_forward.z + R*m_camera.m_right.z + U*m_camera.m_up.z;
  tmp.m_dir = vec3(X,Y,Z);
  return tmp;
}

double distance(dvec3 a, dvec3 b)
{
	return sqrt( (a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y) + (a.z-b.z)*(a.z-b.z) );
}

bool CTracer::RSIntersection(SRay ray, double &dfs)
{														//dfc - distance from start
	double l = ray.m_dir.x;
	double m = ray.m_dir.y;
	double n = ray.m_dir.z;
	
	double x = ray.m_start.x - m_pScene->centre.x;
	double y = ray.m_start.y - m_pScene->centre.y;
	double z = ray.m_start.z - m_pScene->centre.z;

	double A = l*l + m*m + n*n;
	double B = (2*x*l + 2*y*m + 2*z*n);
	double C = (x*x + y*y + z*z - m_pScene->BlackHoleRad*m_pScene->BlackHoleRad);

	double D = B*B - 4*A*C;

	if (D<0) return FALSE;

	double sqrtd = sqrt(D);

	double t1 = (-B + sqrtd)/(2*A);
	double t2 = (-B - sqrtd)/(2*A);
	double t;

	t = min(t1, t2);

	x = ray.m_start.x + t*l;
	y = ray.m_start.y + t*m;
	z = ray.m_start.z + t*n;

	dvec3 dot(x, y, z);

	dfs = distance(dot, ray.m_start);

	return TRUE; 
}


bool CTracer::RPIntersection(SRay ray, double &tResult, double &dfs) //Plane is Oxy
{															//dfs - distance from start

	if (ray.m_dir.z==0) return FALSE;
	double t = (m_pScene->centre.z - ray.m_start.z)/ray.m_dir.z;

	if (t<0) return FALSE;

	dvec3 dot = ray.m_start + t*ray.m_dir;

	tResult = t;
	double r = distance(dot, m_pScene->centre);
	dfs = distance(dot, ray.m_start);

	if ( (r < m_pScene->DiskMaxRad) && (r > m_pScene->DiskMinRad) ) return TRUE;

	return FALSE;
}

glm::vec3 CTracer::TraceRay(SRay ray)
{
	bool I1 = 0, I2 = 0;
	double Pi = 3.141592653589;
	vec3 color (0.4, 0.4, 0.4);
	double t;
	double r1 = 0, r2 = 0;
	double xOffset = m_pScene->centre.x - m_pScene->DiskMaxRad;
	double yOffset = m_pScene->centre.y - m_pScene->DiskMaxRad;

	const double c = 3e+8;
	const double dt = (m_pScene->BlackHoleRad/(13*c));
	const double G = 6.674e-11;
	const double M = 8.57e+36;

	for (;;)
	{
		I1 = RSIntersection(ray, r1);
		if (I1) {I1 = (I1) && (r1 <= dt*c);}
		I2 = RPIntersection(ray, t, r2);
		if (I2) {I2 = (I2) && (r2 <= dt*c);}

		if (!I1 && !I2 && (distance(m_camera.m_pos, ray.m_start)>(distance(m_camera.m_pos, m_pScene->centre)) + 2*m_pScene->DiskMaxRad) ) 
		{
			double fi = atan2(ray.m_dir.y, ray.m_dir.x) + Pi;
			double psi = atan(ray.m_dir.z/sqrt(ray.m_dir.y*ray.m_dir.y + ray.m_dir.x*ray.m_dir.x)) + Pi/2;

			int i = (m_BackColmatr.size()/Pi)*psi;
			int j = (m_BackColmatr[0].size()/(6*Pi))*fi;

			return vec3(m_BackColmatr[i][j*3]/255, m_BackColmatr[i][j*3 + 1]/255, m_BackColmatr[i][j*3 + 2]/255);//}
		}
		if (I1 && !I2) { return vec3(0, 0, 0);}
		if (I2 && !I1) 
		{
			double x = ray.m_start.x + t*ray.m_dir.x;
			double y = ray.m_start.y + t*ray.m_dir.y;

			int i = (x - xOffset)*m_DiskColmatr.size()/(2*m_pScene->DiskMaxRad);
			int j = (y - yOffset)*m_DiskColmatr.size()/(2*m_pScene->DiskMaxRad); 

			float r = (float) m_DiskColmatr[i][j*4]/255;
			float g = (float) m_DiskColmatr[i][j*4+1]/255;
			float b = (float) m_DiskColmatr[i][j*4+2]/255;
			float alpha = (float) m_DiskColmatr[i][j*4+3]/255;

			if (alpha!=0) return vec3(r, g, b);
		}
		if (I1 && I2)
		{
			if (r1<r2) { return vec3(0, 0, 0);}
			else if (r2<r1)
			{
				double x = ray.m_start.x + t*ray.m_dir.x;
				double y = ray.m_start.y + t*ray.m_dir.y;

				int i = (x - xOffset)*m_DiskColmatr.size()/(2*m_pScene->DiskMaxRad);
				int j = (y - yOffset)*m_DiskColmatr.size()/(2*m_pScene->DiskMaxRad);

				float r = (float) m_DiskColmatr[i][j*4]/255;
				float g = (float) m_DiskColmatr[i][j*4+1]/255;
				float b = (float) m_DiskColmatr[i][j*4+2]/255;
				float alpha = (float) m_DiskColmatr[i][j*4+3]/255;

				if (alpha!=0) return vec3(r, g, b);
			}
		}

		dvec3 a = (m_pScene->centre-ray.m_start)*G*M/(distance(m_pScene->centre, ray.m_start)*distance(m_pScene->centre, ray.m_start)*distance(m_pScene->centre, ray.m_start));

		ray.m_start += ray.m_dir*dt + a*(dt*dt/2);
		ray.m_dir += a*dt;
		ray.m_dir = c*normalize(ray.m_dir);

	}

	return color;
}


void CTracer::RenderImage(int xRes, int yRes, int xCamPos, int yCamPos, int zCamPos)
{
  // Reading input texture sample
  CImage* pImage = LoadImageFromFile("data/disk_32.png");

  if(pImage->GetBPP() == 32)
  {
    auto pData = (unsigned char*)pImage->GetBits();
    auto pCurrentLine = pData;
    int pitch = pImage->GetPitch();

    for(int i = 0; i < pImage->GetHeight(); i++) // Image lines
    {
	  std::vector<unsigned char> line;
      for(int j = 0; j < pImage->GetWidth(); j++) // Pixels in line
      {
        unsigned char b = pCurrentLine[j * 4];
        unsigned char g = pCurrentLine[j * 4 + 1];
        unsigned char r = pCurrentLine[j * 4 + 2];
        unsigned char alpha = pCurrentLine[j * 4 + 3];
		line.push_back(r);
		line.push_back(g);
		line.push_back(b);
		line.push_back(alpha);
      }
	  pCurrentLine += pitch;
	  m_DiskColmatr.push_back(line);
    }
  }

   pImage = LoadImageFromFile("data/stars.jpg");
 
    for(int i = 0; i < pImage->GetHeight(); i++) // Image lines
    {
	  std::vector<double> line;
      for(int j = 0; j < pImage->GetWidth(); j++) // Pixels in line
      {
		  double r = GetRValue(pImage->GetPixel(j, i));
		  double g = GetGValue(pImage->GetPixel(j, i));
		  double b = GetBValue(pImage->GetPixel(j, i));
		line.push_back(r);
		line.push_back(g);
		line.push_back(b);
	  }
	  m_BackColmatr.push_back(line);
    }
  
 // Rendering
 
  m_camera.m_resolution = uvec2(xRes, yRes);
  m_camera.m_pixels.resize(xRes * yRes);
  ////////////////////////Set! Camera! Motor!
  double Pi = 3.141592653589;

  m_camera.m_viewAngle = vec2(Pi/4*(xRes/yRes), Pi/4);
  m_camera.m_pos = vec3(xCamPos, yCamPos, zCamPos);
  m_camera.m_pos *= 1e+9;

  double X, Y, Z;

  double F = xRes/(2*tan(m_camera.m_viewAngle.x));

  m_camera.m_forward = F*normalize(m_pScene->centre - m_camera.m_pos)*1e+6;

  if ( (m_camera.m_forward.x == 0) && (m_camera.m_forward.y == 0) ) 
  {
	  m_camera.m_right = vec3(xRes, 0, 0);
	  m_camera.m_right = m_camera.m_right*1e+6;
	  m_camera.m_up = vec3(0, yRes, 0);
	  m_camera.m_up = m_camera.m_up*1e+6;
  }
  else 
  {
	  X = -m_camera.m_forward.y;
	  Y = m_camera.m_forward.x;
	  double Mod = std::sqrt(X*X+Y*Y);

	  m_camera.m_right = vec3((X/Mod)*xRes, (Y/Mod)*xRes, 0);
	  m_camera.m_right *= 1e+6;

	  X = -m_camera.m_forward.x*m_camera.m_forward.z;
	  Y = -m_camera.m_forward.y*m_camera.m_forward.z;
	  Z = m_camera.m_forward.x*m_camera.m_forward.x + m_camera.m_forward.y*m_camera.m_forward.y;
	  Mod = std::sqrt(X*X + Y*Y + Z*Z);

	  m_camera.m_up = vec3((X/Mod)*yRes, (Y/Mod)*yRes, (Z/Mod)*yRes);
	  m_camera.m_up *= 1e+6;
  }
  ////////////////////////////

  for(int i = 0; i < yRes; i++)
  {
    for(int j = 0; j < xRes; j++)
    {
      SRay ray = MakeRay(uvec2(j, i));
      m_camera.m_pixels[i * xRes + j] = TraceRay(ray);
    }
  }
}

void CTracer::SaveImageToFile(std::string fileName)
{
  CImage image;

  int width = m_camera.m_resolution.x;
  int height = m_camera.m_resolution.y;

  image.Create(width, height, 24);
    
	int pitch = image.GetPitch();
	unsigned char* imageBuffer = (unsigned char*)image.GetBits();

	if (pitch < 0)
	{
		imageBuffer += pitch * (height - 1);
		pitch =- pitch;
	}

	int i, j;
	int imageDisplacement = 0;
	int textureDisplacement = 0;

	for (i = 0; i < height; i++)
	{
    for (j = 0; j < width; j++)
    {
      vec3 color = m_camera.m_pixels[textureDisplacement + j];

      imageBuffer[imageDisplacement + j * 3] = clamp(color.b, 0.0f, 1.0f) * 255.0f;
      imageBuffer[imageDisplacement + j * 3 + 1] = clamp(color.g, 0.0f, 1.0f) * 255.0f;
      imageBuffer[imageDisplacement + j * 3 + 2] = clamp(color.r, 0.0f, 1.0f) * 255.0f;
    }

		imageDisplacement += pitch;
		textureDisplacement += width;
	}

  image.Save(fileName.c_str());
	image.Destroy();
}

CImage* CTracer::LoadImageFromFile(std::string fileName)
{
  CImage* pImage = new CImage;

  if(SUCCEEDED(pImage->Load(fileName.c_str())))
    return pImage;
  else
  {
    delete pImage;
    return NULL;
  }
}