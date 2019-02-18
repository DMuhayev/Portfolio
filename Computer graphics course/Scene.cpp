#include "Scene.h"
#include "glm/glm.hpp"

CScene::CScene(double x, double y, double z, double bhr, double minr, double maxr)
{
	centre = glm::dvec3(x, y, z);
	BlackHoleRad = bhr;
	DiskMaxRad = maxr;
	DiskMinRad = minr;
}