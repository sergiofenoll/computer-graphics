//
// Created by sergio on 21/02/17.
//

#ifndef PROJECT_LINE2D_H
#define PROJECT_LINE2D_H

#include <iostream>
#include <vector>
#include "easy_image.hh"

class Point2D{
public:
    double x;
    double y;

    Point2D(){};
    Point2D(double x, double y){Point2D::x=x; Point2D::y=y;};
};

class Line2D{
public:
    Point2D p1;
    Point2D p2;
    img::Color color;

    Line2D(){};
    Line2D(Point2D p1, Point2D p2, img::Color c = img::Color(0, 0, 0)){Line2D::p1=p1; Line2D::p2=p2; Line2D::color=c;};
};

typedef std::vector<Line2D> Lines2D;

#endif //PROJECT_LINE2D_H
