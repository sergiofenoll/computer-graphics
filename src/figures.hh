#ifndef PROJECT_3DFIGURES_H
#define PROJECT_3DFIGURES_H

#include <iostream>
#include <cmath>
#include <stack>
#include <vector>
#include <algorithm>
#include <assert.h>
#include "easy_image.hh"
#include "l_parser.hh"
#include "vector.hh"
#include "color.hh"
#include "ini_configuration.hh"

const double PI = M_PI;

namespace fig {

    class Face {
    public: // TODO: Change to private after testing
        std::vector<unsigned int> point_indexes;
    public:
        Face();

        Face(const std::vector<unsigned int>& point_indexes);

        bool operator==(const Face& rhs) const;

        bool operator!=(const Face& rhs) const;

        Face& operator=(const Face& rhs);

        Face& operator=(const std::vector<unsigned int>& rhs);

        unsigned int& operator[](const std::size_t idx);

        unsigned int operator[](const std::size_t idx) const;
    };

    typedef std::vector<Face> Faces;

    class Figure {
    public: // TODO: Change to protected after testing
        Vectors3D points;
        Faces faces;
        col::Color ambient_reflection = {1, 1, 1};
        col::Color diffuse_reflection = {0, 0, 0};
        col::Color specular_reflection = {0, 0, 0};
        double reflection_coefficient;
    public:
        bool operator==(const Figure& rhs);

        bool operator!=(const Figure& rhs);

        Figure& operator=(const Figure& rhs);

        void set_point_at(std::size_t idx, const Vector3D& point);

        Vector3D get_point_at(std::size_t idx) const;

        void push_back_point(const Vector3D& point);

        void set_face_at(std::size_t idx, const Face& face);

        Face get_face_at(std::size_t idx) const;

        void push_back_face(const Face& face);

        void set_ambient(const col::Color& ambient);

        col::Color get_ambient() const;

        void set_diffuse(const col::Color& diffuse);

        col::Color get_diffuse() const;

        void set_specular(const col::Color& specular);

        col::Color get_specular() const;

        void set_reflection(const double& reflection);

        double get_reflection() const;

        void apply_transformation(const Matrix& m);

        void triangulate();
    };

    class LineDrawing : public Figure {
        LineDrawing(const ini::Configuration &configuration, unsigned int i);
    };

    class Cube : public Figure {
    public:
        Cube();

        void generate_menger(const unsigned int& max_it);
    };

    class Tetrahedron : public Figure {
    public:
        Tetrahedron();
    };

    class Octahedron : public Figure {
    public:
        Octahedron();
    };

    class Icosahedron : public Figure {
    public:
        Icosahedron();
    };

    class Dodecahedron : public Figure {
    public:
        Dodecahedron();
    };

    class Sphere : public Figure {
    public:
        Sphere(const double& radius, const unsigned int& n);
    };

    class Cone : public Figure {
    public:
        Cone(const double& height, const unsigned int &n);
    };

    class Cylinder : public Figure {
    public:
        Cylinder(const double& height, const unsigned int& n);
    };

    class Torus : public Figure {
    public:
        Torus(const double& r, const double& R, const unsigned int& m, const unsigned int& n);
    };

    class Buckyball : public Figure {
    public:
        Buckyball();
    };

    typedef std::vector<Figure*> Figures;

    void generate_fractal(
            Figure& figure, Figures& figures, const double& scale, unsigned int cur_it, const unsigned int& max_it);
}

//Figure createSpecialSphere(const double &r, const double &R, const unsigned int &m, const unsigned int &n);

fig::Figure createDodecahedron();

Matrix scaleFigure(const double& scale);

Matrix rotateX(const double& angle);

Matrix rotateY(const double& angle);

Matrix rotateZ(const double& angle);

Matrix translate(const Vector3D& vector);

void toPolar(const Vector3D& point, double& r, double& theta, double& phi);

Matrix eyePointTrans(const Vector3D &eyepoint);

// Figure createLineDrawing(const ini::Configuration &configuration, unsigned int i);

#endif
