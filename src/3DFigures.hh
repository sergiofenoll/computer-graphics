#ifndef PROJECT_3DFIGURES_H
#define PROJECT_3DFIGURES_H

#include <iostream>
#include <vector>
#include <cmath>
#include <stack>

#include "vector.hh"
#include "easy_image.hh"
#include "l_parser.hh"
#include "ini_configuration.hh"

const double PI = std::atan(1.0) * 4;

// typedef std::vector<Line2D> Lines2D;

class Face {
public:
    std::vector<int> point_indexes;

    Face(){};
    Face(const std::vector<int> &indexes) {
        Face::point_indexes = indexes;
    }
    Face &operator=(Face const &face) {
        point_indexes = face.point_indexes;
        return (*this);
    }
};

class Figure {
public:
    std::vector<Vector3D> points;
    std::vector<Face> faces;
    img::Color color;
    img::Color colorGr = img::Color();
    bool isGradient = false;

    Figure &operator=(Figure const &figure) {
        points = figure.points;
        faces = figure.faces;
        color = figure.color;
        colorGr = figure.colorGr;
        isGradient = figure.isGradient;
        return (*this);
    }
};

typedef std::vector<Figure> Figures3D;

Matrix scaleFigure(const double &scale);

Matrix rotateX(const double &angle);

Matrix rotateY(const double &angle);

Matrix rotateZ(const double &angle);

Matrix translate(const Vector3D &vector);

void applyTransformation(Figure &fig, const Matrix &m);

void applyTransformation(Figures3D &figs, const Matrix &m);

std::vector<double> toPolar(const Vector3D &point);

Matrix eyePointTrans(const Vector3D &eyepoint);

double X(const Vector3D &p1, const Vector3D &p2, const Vector3D &p3);

double Y(const Vector3D &p1, const Vector3D &p2, const Vector3D &p3);

double Z(const Vector3D &p1, const Vector3D &p2, const Vector3D &p3);

void triangulate(Figure& figure);

Figure createLineDrawing(const ini::Configuration &configuration, unsigned int i);

Figure createCube();

Figure createTetrahedron();

Figure createOctahedron();

Figure createIcosahedron();

Figure createDodecahedron();

Figure createSphere(double &radius, const unsigned int &n);

Figure createCone(const double &h, const unsigned int &n);

Figure createCylinder(const double &h, const unsigned int &n);

Figure createTorus(const double &r, const double &R, const int &m, const int &n);

Figure createSpecialSphere(const double &r, const double &R, const int &m, const int &n);

Figure create3DLSystem();

#endif
