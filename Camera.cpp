#include <iostream>
#include <iomanip>
#include "Camera.h"

Camera::Camera() {}



Camera::Camera(int cameraId,
               int projectionType,
               Vec3 position, Vec3 gaze,
               Vec3 u, Vec3 v, Vec3 w,
               double left, double right, double bottom, double top,
               double near, double far,
               int horRes, int verRes,
               std::string outputFilename)
{

    this->cameraId = cameraId;
    this->projectionType = projectionType;
    this->position = position;
    this->gaze = gaze;
    this->u = u;
    this->v = v;
    this->w = w;
    this->left = left;
    this->right = right;
    this->bottom = bottom;
    this->top = top;
    this->near = near;
    this->far = far;
    this->horRes = horRes;
    this->verRes = verRes;
    this->outputFilename = outputFilename;
}

Camera::Camera(const Camera &other)
{
    this->cameraId = other.cameraId;
    this->projectionType = other.projectionType;
    this->position = other.position;
    this->gaze = other.gaze;
    this->u = other.u;
    this->v = other.v;
    this->w = other.w;
    this->left = other.left;
    this->right = other.right;
    this->bottom = other.bottom;
    this->top = other.top;
    this->near = other.near;
    this->far = other.far;
    this->horRes = other.horRes;
    this->verRes = other.verRes;
    this->outputFilename = other.outputFilename;
}



std::ostream &operator<<(std::ostream &os, const Camera &c)
{
    const char *camType = c.projectionType ? "perspective" : "orthographic";

    os << std::fixed << std::setprecision(6) << "Camera " << c.cameraId << " (" << camType << ") => pos: " << c.position << " gaze: " << c.gaze << std::endl
       << "\tu: " << c.u << " v: " << c.v << " w: " << c.w << std::endl
       << std::fixed << std::setprecision(3) << "\tleft: " << c.left << " right: " << c.right << " bottom: " << c.bottom << " top: " << c.top << std::endl
       << "\tnear: " << c.near << " far: " << c.far << " resolutions: " << c.horRes << "x" << c.verRes << " fileName: " << c.outputFilename;

    return os;
}

Matrix4 Camera::translateCameraToOrigin() {
    Vec3 position = this->position;
    double values[4][4] = {
        {1, 0, 0, -position.x},
        {0, 1, 0, -position.y},
        {0, 0, 1, -position.z},
        {0, 0, 0, 1}
    };

    Matrix4 translationMatrix(values);
    return translationMatrix;
}

Matrix4 Camera::rotateCameraAtOrigin() {
    Vec3 u = this->u;
    Vec3 v = this->v;
    Vec3 w = this->w;

    double values[4][4] = {
        {u.x, u.y, u.z, 0},
        {v.x, v.y, v.z, 0},
        {w.x, w.y, w.z, 0},
        {0, 0, 0, 1}
    };

    Matrix4 rotationMatrix(values);
    return rotationMatrix;
}

Matrix4 Camera::orthographicProjectionMatrix() {
    double values[4][4] = {
        {2.0 / (right - left), 0, 0, -(right + left) / (right - left)},
        {0, 2.0 / (top - bottom), 0, -(top + bottom) / (top - bottom)},
        {0, 0, -2 / (far - near), -(far + near) / (far - near)},
        {0, 0, 0, 1}
    };
    Matrix4 orthProjMatrix(values);
    return orthProjMatrix;
}

Matrix4 Camera::pers2orthProjectionMatrix() {
    double values[4][4] = {
        {near, 0, 0, 0},
        {0, near, 0, 0},
        {0, 0, far + near, far * near},
        {0, 0, -1, 0}
    };
    Matrix4 p2oMatrix(values);
    return p2oMatrix;
}

Matrix4 Camera::perspectiveProjectionMatrix() {
    return multiplyMatrixWithMatrix(this->orthographicProjectionMatrix(), this->pers2orthProjectionMatrix());
}

Matrix4 Camera::cameraTransformationMatrix() {
    return multiplyMatrixWithMatrix(this->rotateCameraAtOrigin(), this->translateCameraToOrigin());
}

Matrix4 Camera::viewportMatrix() {
    double values[4][4] = {
        {horRes / 2.0, 0, 0, (horRes - 1) / 2.0},
        {0, verRes / 2.0, 0, (verRes - 1) / 2.0},
        {0, 0, far + near, far * near},
        {0, 0, 0, 1}
    };
    Matrix4 viewportMatrix(values);
    return viewportMatrix;
}

