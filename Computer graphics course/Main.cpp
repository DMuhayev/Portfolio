#include "Tracer.h"
#include "stdio.h"


void main(int argc, char** argv)
{
  int xRes = 1024;  // Default resolution
  int yRes = 768;
  int xCamPos = 0; //Default camera position
  int yCamPos = 0;
  int zCamPos = 0;
  float k = 0;

  if(argc == 2) // There is input file in parameters
  {
    FILE* file = fopen(argv[1], "r");
    if(file)
    {
      int xResFromFile = 0;
      int yResFromFile = 0;
	  int xCamPosFromFile = 0;
	  int yCamPosFromFile = 0;
	  int zCamPosFromFile = 0;
	  float kFromFile = 0;

      if(fscanf(file, "%d %d %d %d %d %f", &xResFromFile, &yResFromFile, &xCamPosFromFile, &yCamPosFromFile, &zCamPosFromFile, &kFromFile) == 6)
      {
        xRes = xResFromFile;
        yRes = yResFromFile;
		xCamPos = xCamPosFromFile;
		yCamPos = yCamPosFromFile;
		zCamPos = zCamPosFromFile;
		k = kFromFile;
      }
      else
        printf("Invalid config format! Using default parameters.\r\n");

      fclose(file);
    }
    else
      printf("Invalid config path! Using default parameters.\r\n");
  }
  else
    printf("No config! Using default parameters.\r\n");

  CTracer tracer;
  CScene scene(12e+10, 12e+10, 0, 12.7e+9, 12.7e+9, k*12.7e+9);

  tracer.m_pScene = &scene;
  tracer.RenderImage(xRes, yRes, xCamPos, yCamPos, zCamPos);
  tracer.SaveImageToFile("Result.png");
}