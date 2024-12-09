#ifndef __CAMERA_H__
#define __CAMERA_H__
#define ORTOGRAPHIC_PROJECTION 0
#define PERSPECTIVE_PROJECTION 1
#include "Vec3.h"
#include "Matrix4.h"
#include "Helpers.h"



class Camera
{

public:
    int cameraId, projectionType, horRes, verRes; // projectionType=0 for orthographic, projectionType=1 for perspective
    Vec3 position, gaze, u, v, w;
    double left, right, bottom, top, near, far;
    std::string outputFilename;

    Camera();

    Camera(int cameraId,
           int projectionType,
           Vec3 position, Vec3 gaze,
           Vec3 u, Vec3 v, Vec3 w,
           double left, double right, double bottom, double top,
           double near, double far,
           int horRes, int verRes,
           std::string outputFilename);

    Camera(const Camera &other);

    Matrix4 translateCameraToOrigin();
    Matrix4 rotateCameraAtOrigin();
    Matrix4 orthographicProjectionMatrix();
    Matrix4 pers2orthProjectionMatrix();
    Matrix4 perspectiveProjectionMatrix();
    Matrix4 cameraTransformationMatrix();

    Matrix4 viewportMatrix();

    friend std::ostream &operator<<(std::ostream &os, const Camera &c);
};

#endif