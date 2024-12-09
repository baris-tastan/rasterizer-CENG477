#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>
#include <math.h>
#include "tinyxml2.h"
#include "Triangle.h"
#include "Helpers.h"
#include "Scene.h"
#include "Line.h"
#include <unordered_map>
#include <limits>
#define M_PI  3.14159265358979323846
using namespace tinyxml2;
using namespace std;

/*
	Parses XML file
*/
Scene::Scene(const char *xmlPath)
{
	const char *str;
	XMLDocument xmlDoc;
	XMLElement *xmlElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode *rootNode = xmlDoc.FirstChild();

	// read background color
	xmlElement = rootNode->FirstChildElement("BackgroundColor");
	str = xmlElement->GetText();
	sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

	// read culling
	xmlElement = rootNode->FirstChildElement("Culling");
	if (xmlElement != NULL)
	{
		str = xmlElement->GetText();

		if (strcmp(str, "enabled") == 0)
		{
			this->cullingEnabled = true;
		}
		else
		{
			this->cullingEnabled = false;
		}
	}

	// read cameras
	xmlElement = rootNode->FirstChildElement("Cameras");
	XMLElement *camElement = xmlElement->FirstChildElement("Camera");
	XMLElement *camFieldElement;
	while (camElement != NULL)
	{
		Camera *camera = new Camera();

		camElement->QueryIntAttribute("id", &camera->cameraId);

		// read projection type
		str = camElement->Attribute("type");

		if (strcmp(str, "orthographic") == 0)
		{
			camera->projectionType = ORTOGRAPHIC_PROJECTION;
		}
		else
		{
			camera->projectionType = PERSPECTIVE_PROJECTION;
		}

		camFieldElement = camElement->FirstChildElement("Position");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->position.x, &camera->position.y, &camera->position.z);

		camFieldElement = camElement->FirstChildElement("Gaze");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->gaze.x, &camera->gaze.y, &camera->gaze.z);

		camFieldElement = camElement->FirstChildElement("Up");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->v.x, &camera->v.y, &camera->v.z);

		camera->gaze = normalizeVec3(camera->gaze);
		camera->u = crossProductVec3(camera->gaze, camera->v);
		camera->u = normalizeVec3(camera->u);

		camera->w = inverseVec3(camera->gaze);
		camera->v = crossProductVec3(camera->u, camera->gaze);
		camera->v = normalizeVec3(camera->v);

		camFieldElement = camElement->FirstChildElement("ImagePlane");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
			   &camera->left, &camera->right, &camera->bottom, &camera->top,
			   &camera->near, &camera->far, &camera->horRes, &camera->verRes);

		camFieldElement = camElement->FirstChildElement("OutputName");
		str = camFieldElement->GetText();
		camera->outputFilename = string(str);

		this->cameras.push_back(camera);

		camElement = camElement->NextSiblingElement("Camera");
	}

	// read vertices
	xmlElement = rootNode->FirstChildElement("Vertices");
	XMLElement *vertexElement = xmlElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (vertexElement != NULL)
	{
		Vec3 *vertex = new Vec3();
		Color *color = new Color();

		vertex->colorId = vertexId;

		str = vertexElement->Attribute("position");
		sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

		str = vertexElement->Attribute("color");
		sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

		this->vertices.push_back(vertex);
		this->colorsOfVertices.push_back(color);

		vertexElement = vertexElement->NextSiblingElement("Vertex");

		vertexId++;
	}

	// read translations
	xmlElement = rootNode->FirstChildElement("Translations");
	XMLElement *translationElement = xmlElement->FirstChildElement("Translation");
	while (translationElement != NULL)
	{
		Translation *translation = new Translation();

		translationElement->QueryIntAttribute("id", &translation->translationId);

		str = translationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

		this->translations.push_back(translation);

		translationElement = translationElement->NextSiblingElement("Translation");
	}

	// read scalings
	xmlElement = rootNode->FirstChildElement("Scalings");
	XMLElement *scalingElement = xmlElement->FirstChildElement("Scaling");
	while (scalingElement != NULL)
	{
		Scaling *scaling = new Scaling();

		scalingElement->QueryIntAttribute("id", &scaling->scalingId);
		str = scalingElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		this->scalings.push_back(scaling);

		scalingElement = scalingElement->NextSiblingElement("Scaling");
	}

	// read rotations
	xmlElement = rootNode->FirstChildElement("Rotations");
	XMLElement *rotationElement = xmlElement->FirstChildElement("Rotation");
	while (rotationElement != NULL)
	{
		Rotation *rotation = new Rotation();

		rotationElement->QueryIntAttribute("id", &rotation->rotationId);
		str = rotationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

		this->rotations.push_back(rotation);

		rotationElement = rotationElement->NextSiblingElement("Rotation");
	}

	// read meshes
	xmlElement = rootNode->FirstChildElement("Meshes");

	XMLElement *meshElement = xmlElement->FirstChildElement("Mesh");
	while (meshElement != NULL)
	{
		Mesh *mesh = new Mesh();

		meshElement->QueryIntAttribute("id", &mesh->meshId);

		// read projection type
		str = meshElement->Attribute("type");

		if (strcmp(str, "wireframe") == 0)
		{
			mesh->type = WIREFRAME_MESH;
		}
		else
		{
			mesh->type = SOLID_MESH;
		}

		// read mesh transformations
		XMLElement *meshTransformationsElement = meshElement->FirstChildElement("Transformations");
		XMLElement *meshTransformationElement = meshTransformationsElement->FirstChildElement("Transformation");

		while (meshTransformationElement != NULL)
		{
			char transformationType;
			int transformationId;

			str = meshTransformationElement->GetText();
			sscanf(str, "%c %d", &transformationType, &transformationId);

			mesh->transformationTypes.push_back(transformationType);
			mesh->transformationIds.push_back(transformationId);

			meshTransformationElement = meshTransformationElement->NextSiblingElement("Transformation");
		}

		mesh->numberOfTransformations = mesh->transformationIds.size();

		// read mesh faces
		char *row;
		char *cloneStr;
		int v1, v2, v3;
		XMLElement *meshFacesElement = meshElement->FirstChildElement("Faces");
		str = meshFacesElement->GetText();
		cloneStr = strdup(str);

		row = strtok(cloneStr, "\n");
		while (row != NULL)
		{
			int result = sscanf(row, "%d %d %d", &v1, &v2, &v3);

			if (result != EOF)
			{
				mesh->triangles.push_back(Triangle(v1, v2, v3));
			}
			row = strtok(NULL, "\n");
		}
		mesh->numberOfTriangles = mesh->triangles.size();
		this->meshes.push_back(mesh);

		meshElement = meshElement->NextSiblingElement("Mesh");
	}
}

/*
	Initializes image with background color
*/
void Scene::initializeImage(Camera *camera)
{
	if (this->image.empty())
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			vector<Color> rowOfColors;

			for (int j = 0; j < camera->verRes; j++)
			{
				rowOfColors.push_back(this->backgroundColor);
			}

			this->image.push_back(rowOfColors);
		}
	}
	else
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			for (int j = 0; j < camera->verRes; j++)
			{
				this->image[i][j].r = this->backgroundColor.r;
				this->image[i][j].g = this->backgroundColor.g;
				this->image[i][j].b = this->backgroundColor.b;
			}
		}
	}
}

/*
	If given value is less than 0, converts value to 0.
	If given value is more than 255, converts value to 255.
	Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value)
{
	if (value >= 255.0)
		return 255;
	if (value <= 0.0)
		return 0;
	return (int)(value);
}


/*
	Writes contents of image (Color**) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera)
{
	ofstream fout;

	fout.open(camera->outputFilename.c_str());

	fout << "P3" << endl;
	fout << "# " << camera->outputFilename << endl;
	fout << camera->horRes << " " << camera->verRes << endl;
	fout << "255" << endl;

	for (int j = camera->verRes - 1; j >= 0; j--)
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
		}
		fout << endl;
	}
	fout.close();
}

/*
	Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
	os_type == 1 		-> Ubuntu
	os_type == 2 		-> Windows
	os_type == other	-> No conversion
*/
void Scene::convertPPMToPNG(string ppmFileName, int osType)
{
	string command;

	// call command on Ubuntu
	if (osType == 1)
	{
		command = "./magick " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// call command on Windows
	else if (osType == 2)
	{
		command = "magick " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// default action - don't do conversion
	else
	{
	}
}


/*
	Transformations, clipping, culling, rasterization are done here.
*/
Matrix4 Scene::translationMatrix(int transID) {
    double values[4][4] = {
        {1, 0, 0, (this->translations)[transID-1]->tx},
        {0, 1, 0, (this->translations)[transID-1]->ty},
        {0, 0, 1, (this->translations)[transID-1]->tz},
        {0, 0, 0, 1}
    };
    Matrix4 result(values);
    return result;
}

Matrix4 Scene::scalingMatrix(int scalingID) {
    double values[4][4] = {
        {this->scalings[scalingID-1]->sx, 0, 0, 0},
        {0, this->scalings[scalingID-1]->sy, 0, 0},
        {0, 0, this->scalings[scalingID-1]->sz, 0},
        {0, 0, 0, 1}
    };
    Matrix4 result(values);
    return result;
}

Matrix4 Scene::rotationMatrix(int rotationID) {
    Rotation* chosenRotation = rotations[rotationID-1];
    double cosT = cos((chosenRotation->angle * M_PI) / 180);
    double sinT = sin((chosenRotation->angle * M_PI) / 180);

    double values[4][4] = {
        {1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1}
    };
    values[0][0] = (1 - cosT) * chosenRotation->ux * chosenRotation->ux + cosT;
    values[0][1] = (1 - cosT) * chosenRotation->ux * chosenRotation->uy - (sinT * chosenRotation->uz);
    values[0][2] = (1 - cosT) * chosenRotation->ux * chosenRotation->uz + sinT * chosenRotation->uy;
    values[1][0] = (1 - cosT) * chosenRotation->ux * chosenRotation->uy + (sinT * chosenRotation->uz);
    values[1][1] = (1 - cosT) * chosenRotation->uy * chosenRotation->uy + cosT;
    values[1][2] = (1 - cosT) * chosenRotation->uy * chosenRotation->uz - (sinT * chosenRotation->ux);
    values[2][0] = (1 - cosT) * chosenRotation->ux * chosenRotation->uz - (sinT * chosenRotation->uy);
    values[2][1] = (1 - cosT) * chosenRotation->uy * chosenRotation->uz + (sinT * chosenRotation->ux);
    values[2][2] = (1 - cosT) * chosenRotation->uz * chosenRotation->uz + cosT;
    Matrix4 result(values);
    return result;
}

void orthProjPoints(unordered_map<int,Vec4>& vertexMap, Camera* camera){
	for (auto & p : vertexMap){
			p.second = multiplyMatrixWithVec4(camera->orthographicProjectionMatrix(), p.second);
		}
}

void persProjPoints(unordered_map<int,Vec4>& vertexMap, Camera* camera){
	for (auto & p : vertexMap){
			p.second = multiplyMatrixWithVec4(camera->perspectiveProjectionMatrix(), p.second);
		}
}



void Scene::forwardRenderingPipeline(Camera *camera) {
	std::vector<std::vector<double>> depthBuffer(camera->horRes, std::vector<double>(camera->verRes, std::numeric_limits<double>::max())); //verRes horREs yeri?

    for (Mesh* mesh : this->meshes) {

		unordered_map<int,Vec4> vertexMap;

		//find all vertices used in the given mesh, extend them into 4D
		for (Triangle triangle : mesh->triangles){

			for (int i=0; i<3; i++){
				Vec3 tDPoint = *((this->vertices)[triangle.vertexIds[i]-1]);
				Vec4 fourDPoint((tDPoint).x, (tDPoint).y, (tDPoint).z, 1, (tDPoint).colorId);
				vertexMap[triangle.vertexIds[i]] = fourDPoint;
			}
			
			
		}

		
		//these coordinates are such that the camera adjustments are done
		/*for (auto & p : vertexMap){
			p.second = multiplyMatrixWithVec4(multiplyMatrixWithMatrix(camera->rotateCameraAtOrigin(), camera->translateCameraToOrigin()), p.second);
		}*/
		for (auto & p : vertexMap){ // for every point apply the mesh transformations
			for (int i=0; i<mesh->numberOfTransformations; i++){
				if (mesh->transformationTypes[i] == 't'){
					p.second = multiplyMatrixWithVec4(translationMatrix(mesh->transformationIds[i]), p.second);
				}
				else if (mesh->transformationTypes[i] == 'r'){
					
					p.second = multiplyMatrixWithVec4(rotationMatrix(mesh->transformationIds[i]), p.second);
					
				}
				else if (mesh->transformationTypes[i] == 's'){
					
					p.second = multiplyMatrixWithVec4(scalingMatrix(mesh->transformationIds[i]), p.second);
					
				}
			
			}

		}
		for (auto & p : vertexMap){
			p.second = multiplyMatrixWithVec4(multiplyMatrixWithMatrix(camera->rotateCameraAtOrigin(), camera->translateCameraToOrigin()), p.second);
		}
		
		if (camera->projectionType == 0){
			orthProjPoints(vertexMap, camera);
		}

		else //(camera->projectionType == 1)
		{
			persProjPoints(vertexMap, camera);
		}


		if (this->cullingEnabled && (mesh->type==0)){ //wireframe with culling (clipping)
			for (Triangle triangle : mesh->triangles){

				Vec4 v2=vertexMap[triangle.vertexIds[2]];
				Vec3 v2_to3D(v2.x,v2.y,v2.z,v2.colorId);

				Vec4 v1=vertexMap[triangle.vertexIds[1]];
				Vec3 v1_to3D(v1.x,v1.y,v1.z,v1.colorId);

				Vec4 v0=vertexMap[triangle.vertexIds[0]];
				Vec3 v0_to3D(v0.x,v0.y,v0.z,v0.colorId);

				Vec3 normal = crossProductVec3(subtractVec3(v1_to3D,v0_to3D), subtractVec3(v2_to3D,v0_to3D));
				if(dotProductVec3(normal,v0_to3D)<0) continue;


				Color color_v0= *(this->colorsOfVertices[v0.colorId-1]);		
				Color color_v1= *(this->colorsOfVertices[v1.colorId-1]);
				Color color_v2= *(this->colorsOfVertices[v2.colorId-1]);

				Line line0(v0,v1,color_v0,color_v1);
				Line line1(v1,v2,color_v1,color_v2);
				Line line2(v2,v0,color_v2,color_v0);
				line0.v0.perspectiveDivide(line0.v0);
				line0.v1.perspectiveDivide(line0.v1);
				line1.v0.perspectiveDivide(line1.v0);
				line1.v1.perspectiveDivide(line1.v1);
				line2.v0.perspectiveDivide(line2.v0);
				line2.v1.perspectiveDivide(line2.v1);

				//Liang Barsky updates the vertices of lines
				bool line0_visible = Liang_Barsky(line0,-1,1,-1,1,-1,1);
				bool line1_visible = Liang_Barsky(line1,-1,1,-1,1,-1,1);
				bool line2_visible = Liang_Barsky(line2,-1,1,-1,1,-1,1);

				if(line0_visible) {
					line0.v0=multiplyMatrixWithVec4(camera->viewportMatrix(), line0.v0);
					line0.v1=multiplyMatrixWithVec4(camera->viewportMatrix(), line0.v1);

					lineRasterization(line0,depthBuffer,camera);
				}

				if(line1_visible) {
					line1.v0=multiplyMatrixWithVec4(camera->viewportMatrix(), line1.v0);
					line1.v1=multiplyMatrixWithVec4(camera->viewportMatrix(), line1.v1);
					lineRasterization(line1,depthBuffer,camera);
				}
				
				if(line2_visible) {
					line2.v0=multiplyMatrixWithVec4(camera->viewportMatrix(), line2.v0);
					line2.v1=multiplyMatrixWithVec4(camera->viewportMatrix(), line2.v1);
					lineRasterization(line2,depthBuffer,camera);
				}
			}
		}
		else if(mesh->type==0) {//wireframe without culling (clipping)
			//clip codu to do
			for (Triangle triangle : mesh->triangles){
			 	Vec4 v2=vertexMap[triangle.vertexIds[2]];

			 	Vec4 v1=vertexMap[triangle.vertexIds[1]];

			 	Vec4 v0=vertexMap[triangle.vertexIds[0]];

				Color color_v0= *(this->colorsOfVertices[v0.colorId-1]);		
				Color color_v1= *(this->colorsOfVertices[v1.colorId-1]);
				Color color_v2= *(this->colorsOfVertices[v2.colorId-1]);

				Line line0(v0,v1,color_v0,color_v1);
				Line line1(v1,v2,color_v1,color_v2);
				Line line2(v2,v0,color_v2,color_v0);
				line0.v0.perspectiveDivide(line0.v0);
				line0.v1.perspectiveDivide(line0.v1);
				line1.v0.perspectiveDivide(line1.v0);
				line1.v1.perspectiveDivide(line1.v1);
				line2.v0.perspectiveDivide(line2.v0);
				line2.v1.perspectiveDivide(line2.v1);

				//Liang Barsky updates the vertices of lines
				bool line0_visible = Liang_Barsky(line0,-1,1,-1,1,-1,1);
				bool line1_visible = Liang_Barsky(line1,-1,1,-1,1,-1,1);
				bool line2_visible = Liang_Barsky(line2,-1,1,-1,1,-1,1);

				if(line0_visible) {
					line0.v0=multiplyMatrixWithVec4(camera->viewportMatrix(), line0.v0);
					line0.v1=multiplyMatrixWithVec4(camera->viewportMatrix(), line0.v1);

					lineRasterization(line0,depthBuffer,camera);
				}

				if(line1_visible) {
					line1.v0=multiplyMatrixWithVec4(camera->viewportMatrix(), line1.v0);
					line1.v1=multiplyMatrixWithVec4(camera->viewportMatrix(), line1.v1);
					lineRasterization(line1,depthBuffer,camera);
				}
				
				if(line2_visible) {
					line2.v0=multiplyMatrixWithVec4(camera->viewportMatrix(), line2.v0);
					line2.v1=multiplyMatrixWithVec4(camera->viewportMatrix(), line2.v1);
					lineRasterization(line2,depthBuffer,camera);
				}
			}
		}

		else if(this->cullingEnabled  && (mesh->type==1)) {//solid with culling
			for (Triangle triangle : mesh->triangles){
				Vec4 v2=vertexMap[triangle.vertexIds[2]];
				Vec3 v2_to3D(v2.x,v2.y,v2.z,v2.colorId);

				Vec4 v1=vertexMap[triangle.vertexIds[1]];
				Vec3 v1_to3D(v1.x,v1.y,v1.z,v1.colorId);

				Vec4 v0=vertexMap[triangle.vertexIds[0]];
				Vec3 v0_to3D(v0.x,v0.y,v0.z,v0.colorId);

				Vec3 normal = crossProductVec3(subtractVec3(v1_to3D,v0_to3D), subtractVec3(v2_to3D,v0_to3D));
				if(dotProductVec3(normal,v0_to3D)<0) continue;

				//now v0,v1 and v2 are the viewport coordinates
			 	v0 = multiplyMatrixWithVec4(camera->viewportMatrix(), v0.perspectiveDivide(v0));
			 	v1 = multiplyMatrixWithVec4(camera->viewportMatrix(), v1.perspectiveDivide(v1));
			 	v2 = multiplyMatrixWithVec4(camera->viewportMatrix(), v2.perspectiveDivide(v2));

				triangleRasterization(v0,v1,v2,depthBuffer,camera);

				
			}
			
		}

		else{  
		
		 //solid without culling
			for (Triangle triangle : mesh->triangles){
				
			 	Vec4 v2=vertexMap[triangle.vertexIds[2]];

			 	Vec4 v1=vertexMap[triangle.vertexIds[1]];

			 	Vec4 v0=vertexMap[triangle.vertexIds[0]];

			 	//now v0,v1 and v2 are the viewport coordinates
			 	v0 = multiplyMatrixWithVec4(camera->viewportMatrix(), v0.perspectiveDivide(v0));
			 	v1 = multiplyMatrixWithVec4(camera->viewportMatrix(), v1.perspectiveDivide(v1));
			 	v2 = multiplyMatrixWithVec4(camera->viewportMatrix(), v2.perspectiveDivide(v2));

			 	triangleRasterization(v0,v1,v2,depthBuffer, camera);

		}
		
			
			
	}
}
}
bool Scene::isVisible(double den,double num,double& t_E,double& t_L){
	double t=0;
	if(den>0){
		t=num/den;
		if(t > t_L) return false;
		if(t > t_E) t_E = t;
	}
	else if(den<0){
		t=num/den;
		if(t < t_E) return false;
		if(t < t_L) t_L = t;
	}
	else if(num>0) {
		return false;
	}
	return true;
}
//line gönderdim, line coordinatları ve renkleri güncellendi. Triangleların tüm linelarını sırayla gönder, sonra 
//sana cliplenmiş linelar gelecek. O lineları viewportta çizme zamanı. -- midpoint algorithm ile.
//midpoint algo viewporta gönderilmiş lineı alsın ve image'a yazsın.
bool Scene::Liang_Barsky(Line& line, double xMin, double xMax, double yMin, double yMax, double zMin, double zMax){
	double t_E=0;
	double t_L=1;
	bool visible=false;
	Vec4 v1 = line.v1;
	Vec4 v0 = line.v0;
	if( isVisible(v1.x-v0.x, -1 - v0.x, t_E, t_L) && 
		isVisible(v0.x-v1.x, v0.x - 1, t_E, t_L) &&
		isVisible(v1.y-v0.y, -1 - v0.y, t_E, t_L) && 
		isVisible(v0.y-v1.y, v0.y - 1, t_E, t_L) &&
		isVisible(v1.z-v0.z, -1 - v0.z, t_E, t_L)&&
		isVisible(v0.z-v1.z, v0.z - 1, t_E, t_L)) 
	{
			Color dc = Color(line.c1.r - line.c0.r, line.c1.g - line.c0.g, line.c1.b - line.c0.b);
			visible=true;
			if(t_L < 1){
				line.v1.x = line.v0.x + (line.v1.x - line.v0.x)*t_L;
				line.v1.y = line.v0.y + (line.v1.y - line.v0.y)*t_L;
				line.v1.z = line.v0.z + (line.v1.z - line.v0.z)*t_L;
				line.c1.r = line.c0.r + (t_L * dc.r);
				line.c1.g = line.c0.g + (t_L * dc.g);
				line.c1.b = line.c0.b + (t_L * dc.b);
			}

			if(t_E > 0) {
				line.v0.x = line.v0.x + (line.v1.x - line.v0.x)*t_E;
				line.v0.y = line.v0.y + (line.v1.y - line.v0.y)*t_E;
				line.v0.z = line.v0.z + (line.v1.z - line.v0.z)*t_E;
				line.c0.r = line.c0.r + (t_E * dc.r);
				line.c0.g = line.c0.g + (t_E * dc.g);
				line.c0.b = line.c0.b + (t_E * dc.b);

			}

	}
	return visible;
		

}

void Scene::lineRasterization(Line line,std::vector<std::vector<double>> &depthBuffer, Camera* camera){
	int startX=(int) line.v0.x;
	int startY= (int)line.v0.y;
	int endX = (int)line.v1.x;
	int endY = (int)line.v1.y;
	double dx = line.v1.x-line.v0.x;
	double dy = line.v1.y-line.v0.y;
	Color color,dc; double depth,delta ,deltaDepth;
	int slope_sign=1;

	double slope; 
	if(startX==endX) slope=2*((line.v1.y-line.v0.y)); //random val instead of infinity
	else slope= (line.v1.y-line.v0.y) / (line.v1.x-line.v0.x);

	if((slope>=0 && slope<=1) || (slope<=0 && slope>-1)){ //horizontal line

		if(endX < startX) {
			startX=line.v1.x;
			startY= line.v1.y;
			endX = line.v0.x;
			endY = line.v0.y;	
			line.swapLine(line);
			dx= -1*dx;
			dy=-1*dy;
	
		}

		if(endY<startY) slope_sign = -1;


		int y=startY;
		color=line.c0;
		depth=line.v0.z;
		delta=2*(-dy) + slope_sign* dx;
		deltaDepth= (line.v1.z - line.v0.z) / dx;
		dc= (line.c1- line.c0) / dx;


		for(int x= startX; x<=endX ; x++){ //endxe eşit mi değil mi bak ??
			if(!((x<0) || (y<0) || x>camera->horRes-1 || y>camera->verRes-1)){
			if(depthBuffer[x][y]>=depth){
			depthBuffer[x][y]=depth;

			image[x][y].r= color.r;
			image[x][y].g= color.g;
			image[x][y].b= color.b; } }
			if(delta*slope_sign <0){
				y+=slope_sign;
				delta+= 2*((startY-endY) + slope_sign*dx);
				depth+= deltaDepth;
				color=color+ dc;
			}
			else{
				delta+=2*(startY-endY);
				depth+= deltaDepth;
				color=color+ dc;

			}
		}

	}
	else{
		if(endY < startY) {
			startX=(int)line.v1.x;
			startY= (int)line.v1.y;
			endX = (int)line.v0.x;
			endY = (int)line.v0.y;
			line.swapLine(line);
			dx= -1*dx;
			dy=-1*dy;
		}

		if(endX < startX) slope_sign = -1;

		int x=startX;
		color=line.c0;
		depth=line.v0.z;
		delta=2*(dx) + slope_sign*(-dy);
		deltaDepth= (line.v1.z - line.v0.z) / dy;
		dc=(line.c1- line.c0) / dy;

		for(int y=startY; y<=endY ; y++){ //endxe eşit mi değil mi bak ??
			if(!((x<0) || (y<0) || x>camera->horRes-1 || y>camera->verRes-1)){
			if(depthBuffer[x][y]>=depth){
			depthBuffer[x][y]=depth;

			image[x][y].r= color.r;
			image[x][y].g= color.g;
			image[x][y].b= color.b; } }
			if(delta*slope_sign > 0){
				x+=slope_sign;
				delta+= 2*(dx + slope_sign*(-dy));
				depth+= deltaDepth;
				color=color+ dc;
			}
			else{
				delta+=2*dx;
				depth+= deltaDepth;
				color=color+ dc;
			}
		}
	}




}

void Scene::triangleRasterization(Vec4 v0, Vec4 v1, Vec4 v2,std::vector<std::vector<double>> &depthBuffer, Camera* camera){
	//v0,v1 and v2 is already in the viewport. 
	//just calculate all the inside pixels and color them in the scene
	//first check the depth. if it s bigger then the depth buffer, then return. do nothing to image.


	//bounding square of the triangle
	double minX=min(min(v0.x,v1.x),v2.x);
	double minY=min(min(v0.y,v1.y),v2.y);
	double maxX=max(max(v0.x,v1.x),v2.x);
	double maxY=max(max(v0.y,v1.y),v2.y);
	double depth;
	if(minX < 0) minX = 0;

    if(minY < 0) minY = 0;

    if(maxX > camera->horRes -1 )   maxX = camera->horRes -1;

    if(maxY > camera->verRes -1 )   maxY = camera->verRes -1;
	//iterate over each pixel in the bounding square
	for(int i=minX; i<=maxX ;i++){
		for(int j=minY; j<=maxY; j++){

			double alpha = ((v1.y - v2.y)*i + (v2.x - v1.x)*j + v1.x*v2.y - v2.x*v1.y)
                           / ((v1.y - v2.y)*v0.x + (v2.x - v1.x)*v0.y + v1.x*v2.y - v2.x*v1.y);

            double beta = ((v2.y - v0.y)*i + (v0.x - v2.x)*j + v2.x*v0.y - v0.x*v2.y)
                          / ((v1.y - v2.y)*v0.x + (v2.x - v1.x)*v0.y + v1.x*v2.y - v2.x*v1.y);

            double gamma = ((v0.y-v1.y)*i + (v1.x-v0.x)*j + v0.x*v1.y - v0.y*v1.x) / (v2.x*(v0.y-v1.y) + v2.y*(v1.x-v0.x) + v0.x*v1.y - v0.y*v1.x);

			if(alpha >= 0 && beta >= 0 && gamma >= 0){ //inside triangle
				depth= v0.z*alpha + v1.z*beta + v2.z*gamma;
				
				if(depthBuffer[i][j]<depth) continue;

			else { 
					depthBuffer[i][j]=depth;


					Color color_v0=*(this->colorsOfVertices[v0.colorId-1]);		
					Color color_v1= *(this->colorsOfVertices[v1.colorId-1]);
					Color color_v2= *(this->colorsOfVertices[v2.colorId-1]);
					image[i][j].r = color_v0.r *alpha + color_v1.r*beta + color_v2.r*gamma;
					image[i][j].g = color_v0.g *alpha + color_v1.g*beta + color_v2.g*gamma;
					image[i][j].b = color_v0.b *alpha + color_v1.b*beta + color_v2.b*gamma;
				}
			}
		}
	}	
	
}
		

