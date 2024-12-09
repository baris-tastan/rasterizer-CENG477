#ifndef _SCENE_H_
#define _SCENE_H_
#include "Vec3.h"
#include "Vec4.h"
#include "Color.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Camera.h"
#include "Mesh.h"
#include "Matrix4.h"
#include "Line.h"
class Scene
{
public:
	Color backgroundColor;
	bool cullingEnabled;

	std::vector<std::vector<Color> > image;
	std::vector<Camera *> cameras;
	std::vector<Vec3 *> vertices;
	std::vector<Color *> colorsOfVertices;
	std::vector<Scaling *> scalings;
	std::vector<Rotation *> rotations;
	std::vector<Translation *> translations;
	std::vector<Mesh *> meshes;

	Scene(const char *xmlPath);

	void initializeImage(Camera *camera);
	int makeBetweenZeroAnd255(double value);
	void writeImageToPPMFile(Camera *camera);
	void convertPPMToPNG(std::string ppmFileName, int osType);
	void forwardRenderingPipeline(Camera *camera);
	Matrix4 translationMatrix(int transID);
    Matrix4 scalingMatrix(int scalingID);
    Matrix4 rotationMatrix(int rotationID);

	void triangleRasterization(Vec4 v0, Vec4 v1, Vec4 v2,std::vector<std::vector<double>> &depthBuffer, Camera* camera);
	void lineRasterization(Line line,std::vector<std::vector<double>> &depthBuffer, Camera* camera);
	bool Liang_Barsky(Line& line, double xMin, double xMax, double yMin, double yMax, double zMin, double zMax);
	bool isVisible(double den,double num,double& t_E,double& t_L);
};

#endif