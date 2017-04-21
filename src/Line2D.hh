//
// Created by sergio on 21/02/17.
//

#ifndef PROJECT_LINE2D_H
#define PROJECT_LINE2D_H

#include <iostream>
#include <vector>
#include <stack>
#include <limits>
#include "easy_image.hh"
#include "l_parser.hh"
#include "vector.hh"
#include "3DFigures.hh"

enum FigureType {
    Wires,
    ZBuff,
    Trian
};

class Point2D{
public:
    double x;
    double y;
    double z; // DO NOT USE FFS

    Point2D(){};
    Point2D(double x, double y){Point2D::x=x; Point2D::y=y;};
};

class Line2D{
public:
    Point2D p1;
    Point2D p2;
    img::Color color;
    img::Color colorGr;
    bool isGradient;

    Line2D(){};
    Line2D(Point2D p1, Point2D p2,
           img::Color c = img::Color(0, 0, 0),
           img::Color cGr = img::Color(0, 0, 0),
           bool isGradient = false)
        {
            Line2D::p1=p1;
            Line2D::p2=p2;
            Line2D::color=c;
            Line2D::colorGr=cGr;
            Line2D::isGradient = isGradient;
        };
};

typedef std::vector<Point2D> Points2D;
typedef std::vector<Line2D> Lines2D;

Point2D doProjection(const Vector3D& point, const double& d=1.0);

Lines2D doProjection(Figures3D& figures);

img::EasyImage draw2DLines(
        Lines2D& lines, const unsigned int& size, const img::Color& bc, const img::Color& bcGr=img::Color(), bool isGrB=false);

img::EasyImage draw3DLines(
        Figures3D& figures, FigureType& type, const unsigned int& size, const img::Color& bc, const img::Color& bcGr=img::Color(), bool isGrB=false);

img::EasyImage drawZBufLines(
        Lines2D &lines, const unsigned int& size, const img::Color& bc, const img::Color& bcGr=img::Color(), bool isGrB=false);

img::EasyImage drawTriangLines(
        Figures3D& figures, const unsigned int& size, const img::Color& bc, const img::Color& bcGr=img::Color(), bool isGrB=false);

void recursivePrintString(
        LParser::LSystem2D& l_system, std::string& print_string, unsigned int& cur_it, unsigned int& max_it,
        double& alph_angle, double& delt_angle, Lines2D& lines, std::stack<double>& brStack,
        img::Color& lc, img::Color& lcGr, bool isGrL, double& x, double& y);

void recursivePrintString(
        LParser::LSystem3D& l_system, std::string& print_string, unsigned int& cur_it, unsigned int& max_it,
        double& delt_angle, Figure& figure, std::stack<double>& doubleStack, std::stack<Vector3D>& vectorStack,
        double& x, double& y, double& z, Vector3D& H, Vector3D& L, Vector3D& U);

#endif //PROJECT_LINE2D_H
