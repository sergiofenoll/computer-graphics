//
// Created by sergio on 30/04/17.
//

#ifndef COMPGRAPHX_PROJECTION_HH
#define COMPGRAPHX_PROJECTION_HH


#include <vector>
#include "vector.hh"
#include "figures.hh"
#include "color.hh"

namespace prj {
    class Point{
    public: // TODO: Change to private
        double x;
        double y;
        double z; // DO NOT USE FFS
    public:
        Point();

        Point(const double& x, const double& y);

        Point(const Vector3D& point, const double& d = 1.0);

        void project(const Vector3D& point, const double& d = 1.0);
    };

    typedef std::vector<Point> Points;

    class Line{
    public: // TODO: Change to private
        Point p1;
        Point p2;
        col::Color color;
    public:
        Line();

        Line(Point p1, Point p2, col::Color color);
    };

    typedef std::vector<Line> Lines;

    Lines project_figures(const fig::Figures& figures);

    void recursivePrintString(
            LParser::LSystem2D& l_system, std::string& print_string,
            unsigned int& cur_it, unsigned int& max_it, double& alph_angle, double& delt_angle,
            Lines& lines, std::stack<double>& brStack, col::Color& lc, double& x, double& y);

    void recursivePrintString(
            LParser::LSystem3D& l_system, std::string& print_string,
            unsigned int& cur_it, unsigned int& max_it, double& delt_angle,
            fig::Figure& figure, std::stack<double>& doubleStack, std::stack<Vector3D>& vectorStack,
            double& x, double& y, double& z, Vector3D& H, Vector3D& L, Vector3D& U);
};


#endif //COMPGRAPHX_PROJECTION_HH
