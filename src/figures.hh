#ifndef PROJECT_3DFIGURES_H
#define PROJECT_3DFIGURES_H

#include <iostream>
#include <cmath>
#include <stack>
#include <vector>
#include <algorithm>
#include <assert.h>
#include "l_parser.hh"
#include "vector.hh"

const double PI = M_PI;

namespace fig {

    /**
     * @brief A class that stores RGB color values.
     *
     * The values are stored as doubles between 0 and 1.
     */
    class Color {
    private:
        /**
         * @brief The intensity of the red color component.
         */
        double red;

        /**
         * @brief The intensity of the green color component.
         */
        double green;

        /**
         * @brief The intensity of the blue color component.
         */
        double blue;
    public:
        /**
         * @brief Default constructor
         */
        Color();

        /**
         * @brief Construct a color with the given intensities.
         * @param red The red color component.
         * @param green The green color component.
         * @param blue The blue color component.
         */
        Color(const double& red, const double& green, const double& blue);

        /**
         * @brief Equality operator overload by comparing the red, green and blue intensities of this color and @param rhs.
         * @param rhs The color that is compared to this.
         * @return true if both colors are the same, false otherwise.
         */
        bool operator==(const Color& rhs) const;

        /**
         * @brief Inequality operator overload by comparing the red, green and blue intensities of this color and @param rhs.
         * @param rhs The color that is compared to this.
         * @return true if both colors are different, false otherwise.
         */
        bool operator!=(const Color& rhs) const;

        /**
         * @brief Assignment operator overload for a Color object.
         * @param rhs The Color that is copied.
         * @return A reference to this color.
         */
        Color& operator=(const Color& rhs);

        /**
         * @brief Assignment operator overload for a vector of doubles.
         * @param rhs The vector of size 3 containing doubles whose values will be used as red, green and blue intensities.
         * @return A reference to this color.
         */
        Color& operator=(const std::vector<double>& rhs);

        /**
         * @brief Addition operator overload by adding the red, green and blue intensities of this color and @param rhs.
         * @param rhs The color that will be added to this color.
         * @return A reference to this color.
         */
        Color& operator+(const Color& rhs);

        Color& operator+=(const Color& rhs);

        /**
         * @brief Mulitplication operator overload by mulitplying the red, green and blue intensities by @param rhs.
         * @param rhs The double that will be mulitplied to this color.
         * @return A reference to this color.
         */
        Color& operator*(const double& rhs);

        Color& operator*=(const double& rhs);

        void set_red_value(const double& red);

        double get_red_value() const;

        uint8_t red_int_value() const;

        void set_green_value(const double& green);

        double get_green_value() const;

        uint8_t green_int_value() const;

        void set_blue_value(const double& blue);

        double get_blue_value() const;

        uint8_t blue_int_value() const;
    };

    class Point2D{
    public:
        double x;
        double y;
        double z; // DO NOT USE FFS

        Point2D(){};
        Point2D(double x, double y){Point2D::x=x; Point2D::y=y;};
    };

    typedef std::vector<Point2D> Points2D;

    Point2D doProjection(const Vector3D& point, const double& d=1.0);

    class Line2D{
    public:
        Point2D p1;
        Point2D p2;
        Color color;
        Color colorGr;
        bool isGradient;

        Line2D(){};
        Line2D(Point2D p1, Point2D p2,
               Color c = Color(0, 0, 0),
               Color cGr = Color(0, 0, 0),
               bool isGradient = false)
        {
            Line2D::p1=p1;
            Line2D::p2=p2;
            Line2D::color=c;
            Line2D::colorGr=cGr;
            Line2D::isGradient = isGradient;
        };
    };

    typedef std::vector<Line2D> Lines2D;

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
        fig::Color ambient_reflection = {1, 1, 1};
        fig::Color diffuse_reflection = {0, 0, 0};
        fig::Color specular_reflection = {0, 0, 0};
        double reflection_coefficient;
        fig::Color colorGr = fig::Color();
        bool isGradient = false;
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

        void set_ambient(const fig::Color& ambient);

        fig::Color get_ambient() const;

        void set_diffuse(const fig::Color& diffuse);

        fig::Color get_diffuse() const;

        void set_specular(const fig::Color& specular);

        fig::Color get_specular() const;

        void set_reflection(const double& reflection);

        double get_reflection() const;

        void apply_transformation(const Matrix& m);


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

    typedef std::vector<Figure*> Figures3D;

    Lines2D doProjection(Figures3D& figures);

    void recursivePrintString(
            LParser::LSystem2D& l_system, std::string& print_string, unsigned int& cur_it, unsigned int& max_it,
            double& alph_angle, double& delt_angle, Lines2D& lines, std::stack<double>& brStack,
            Color& lc, Color& lcGr, bool isGrL, double& x, double& y);

    void recursivePrintString(
            LParser::LSystem3D& l_system, std::string& print_string, unsigned int& cur_it, unsigned int& max_it,
            double& delt_angle, Figure& figure, std::stack<double>& doubleStack, std::stack<Vector3D>& vectorStack,
            double& x, double& y, double& z, Vector3D& H, Vector3D& L, Vector3D& U);
}

//Figure createSpecialSphere(const double &r, const double &R, const unsigned int &m, const unsigned int &n);

Matrix scaleFigure(const double& scale);

Matrix rotateX(const double& angle);

Matrix rotateY(const double& angle);

Matrix rotateZ(const double& angle);

Matrix translate(const Vector3D& vector);

void toPolar(const Vector3D& point, double& r, double& theta, double& phi);

Matrix eyePointTrans(const Vector3D &eyepoint);

double X(const Vector3D &p1, const Vector3D &p2, const Vector3D &p3);

double Y(const Vector3D &p1, const Vector3D &p2, const Vector3D &p3);

double Z(const Vector3D &p1, const Vector3D &p2, const Vector3D &p3);

//void triangulate(Figure& figure);

//void recursiveGenerateFractal(
//        Figure& figure, Figures3D& figures, const double& scale, unsigned int cur_it, const unsigned int& max_it);

// Figure createLineDrawing(const ini::Configuration &configuration, unsigned int i);

#endif
