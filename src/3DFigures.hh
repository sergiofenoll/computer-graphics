#ifndef PROJECT_3DFIGURES_H
#define PROJECT_3DFIGURES_H

#include <iostream>
#include <vector>
#include <cmath>

#include "vector.hh"
#include "easy_image.hh"

class Face {
public:
    std::vector<int> point_indexes;
};

class Figure {
public:
    std::vector<Vector3D> points;
    std::vector<Face> faces;
    img::Color color;
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

#endif
