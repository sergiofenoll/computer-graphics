#include "figures.hh"

namespace fig {
    Point2D doProjection(const Vector3D& point, const double& d){
        Point2D p;
        p.x = (d * point.x) / -(point.z);
        p.y = (d * point.y) / -(point.z);
        p.z = point.z;
        return p;
    }

    Lines2D doProjection(Figures3D& figures){
        Lines2D lines;
        for (Figure* fig : figures){
            for (Face &face : fig->faces){
                Points2D points;
                for (int index : face.point_indexes){
                    Point2D p = doProjection(fig->points[index]);
                    points.push_back(p);
                }
                for (int i=0; i<points.size(); i++) {
                    Line2D line;
                    line.p1 = points[i];
                    line.p2 = points[(i + 1) % points.size()];
                    line.color = Color(fig->ambient_reflection.red_int_value(),
                                            fig->ambient_reflection.green_int_value(),
                                            fig->ambient_reflection.blue_int_value());
                    // line.colorGr = fig.colorGr; TODO: Uncomment this
                    line.isGradient = fig->isGradient;
                    lines.push_back(line);
                }
            }
        }
        return lines;
    }

    Color::Color() :
    red(0), green(0), blue(0) {}
    
    Color::Color(const double& red, const double& green, const double& blue) :
    red(red), green(green), blue(blue) {
        Color::red = Color::red <= 1.0 ? Color::red : 1.0;
        Color::green = Color::green <= 1.0 ? Color::green : 1.0;
        Color::blue = Color::blue <= 1.0 ? Color::blue : 1.0;
    }
    
    bool Color::operator==(const Color& rhs) const {
        return (*this).red == rhs.red &&
               (*this).green == rhs.green &&
               (*this).blue == rhs.blue;
    }
    
    bool Color::operator!=(const Color& rhs) const {
        return !((*this) == rhs);
    }
    
    Color& Color::operator=(const Color& rhs) {
        if ((*this) != rhs) {
            (*this).red = rhs.red;
            (*this).green = rhs.green;
            (*this).blue = rhs.blue;
        }
        return (*this);
    }
    
    Color& Color::operator=(const std::vector<double>& rhs) {
        assert(rhs.size() == 3);
        (*this).red = rhs[0] <= 1.0 ? rhs[0] : 1.0;
        (*this).green = rhs[1] <= 1.0 ? rhs[1] : 1.0;
        (*this).blue = rhs[2] <= 1.0 ? rhs[2] : 1.0;
        return (*this);
    }
    
    Color& Color::operator+(const Color& rhs) {
        (*this).red = (*this).red + rhs.red <= 1.0 ? (*this).red + rhs.red : 1.0;
        (*this).green = (*this).green + rhs.green <= 1.0 ? (*this).green + rhs.green : 1.0;
        (*this).blue = (*this).blue + rhs.blue <= 1.0 ? (*this).blue + rhs.blue : 1.0;
        return (*this);
    }
    
    Color& Color::operator+=(const Color& rhs) {
        return (*this) + rhs;
    }
    
    Color& Color::operator*(const double& rhs) {
        (*this).red = (*this).red * rhs <= 1.0 ? (*this).red * rhs: 1.0;
        (*this).green = (*this).green * rhs <= 1.0 ? (*this).green * rhs : 1.0;
        (*this).blue = (*this).blue * rhs <= 1.0 ? (*this).blue * rhs : 1.0;
        return (*this);
    }
    
    Color& Color::operator*=(const double& rhs) {
        return (*this) * rhs;
    }
    
    void Color::set_red_value(const double& red) {
        assert(red >= 0 and red <= 1);
        (*this).red = red;
    }
    
    double Color::get_red_value() const {
        return (*this).red;
    }
    
    uint8_t Color::red_int_value() const{
        return (uint8_t) ((*this).red * 255);
    }
    
    void Color::set_green_value(const double& green) {
        assert(green >= 0 and green <= 1);
        (*this).green = green;
    }
    
    double Color::get_green_value() const {
        return (*this).green;
    }
    
    uint8_t Color::green_int_value() const{
        return (uint8_t) ((*this).green * 255);
    }
    
    void Color::set_blue_value(const double& blue) {
        assert(blue >= 0 and blue <= 1);
        Color::blue = blue;
    }
    
    double Color::get_blue_value() const {
        return (*this).blue;
    }
    
    uint8_t Color::blue_int_value() const{
        return (uint8_t) ((*this).blue * 255);
    }
    
    Face::Face() {}
    
    Face::Face(const std::vector<unsigned int>& point_indexes) :
    point_indexes(point_indexes){}
    
    bool Face::operator==(const Face& rhs) const {
        return (*this).point_indexes == rhs.point_indexes;
    }
    
    bool Face::operator!=(const Face& rhs) const {
        return !((*this) == rhs);
    }
    
    Face& Face::operator=(const Face& rhs) {
        if ((*this) != rhs) {
            point_indexes = rhs.point_indexes;
        }
        return (*this);
    }
    
    Face& Face::operator=(const std::vector<unsigned int>& rhs) {
        Face::point_indexes = rhs;
        return (*this);
    }
    
    unsigned int& Face::operator[](const std::size_t idx) {
        return point_indexes[idx];
    }
    
    unsigned int Face::operator[](const std::size_t idx) const {
        return point_indexes[idx];
    }
    
    bool Figure::operator==(const Figure& rhs) {
        return ((*this).points == rhs.points) and
               ((*this).faces == rhs.faces);
    }
    
    bool Figure::operator!=(const Figure& rhs) {
        return !((*this) == rhs);
    }
    
    Figure& Figure::operator=(const Figure& rhs) {
        if ((*this) != rhs) {
            points = rhs.points;
            faces = rhs.faces;
            ambient_reflection = rhs.ambient_reflection;
            diffuse_reflection = rhs.diffuse_reflection;
            specular_reflection = rhs.specular_reflection;
            reflection_coefficient = rhs.reflection_coefficient;
            colorGr = rhs.colorGr;
            isGradient = rhs.isGradient;
        }
        return (*this);
    }
    
    void Figure::set_point_at(std::size_t idx, const Vector3D& point) {
        assert(point.is_point());
        (*this).points[idx] = point;
    }
    
    Vector3D Figure::get_point_at(std::size_t idx) const {
        return (*this).points[idx];
    }
    
    void Figure::push_back_point(const Vector3D& point) {
        points.push_back(point);
    }
    
    void Figure::set_face_at(std::size_t idx, const Face& face) {
        (*this).faces[idx] = face;
    }
    
    Face Figure::get_face_at(std::size_t idx) const {
        return (*this).faces[idx];
    }
    
    void Figure::push_back_face(const Face& face) {
        faces.push_back(face);
    }
    
    void Figure::set_ambient(const Color& ambient) {
        (*this).ambient_reflection = ambient;
    }
    
    Color Figure::get_ambient() const {
        return (*this).ambient_reflection;
    }
    
    void Figure::set_diffuse(const Color& diffuse) {
        (*this).diffuse_reflection = diffuse;
    }
    
    Color Figure::get_diffuse() const {
        return (*this).diffuse_reflection;
    }
    
    void Figure::set_specular(const Color& specular) {
        (*this).specular_reflection;
    }
    
    Color Figure::get_specular() const {
        return (*this).specular_reflection;
    }
    
    void Figure::set_reflection(const double& reflection) {
        (*this).reflection_coefficient = reflection;
    }
    
    double Figure::get_reflection() const {
        return (*this).reflection_coefficient;
    }
    
    void Figure::apply_transformation(const Matrix &m){
        for(Vector3D &point : (*this).points){
            point *= m;
        }
    }
    
    Cube::Cube() {
        Vector3D p0 = Vector3D::point(1, -1, -1);
        Vector3D p1 = Vector3D::point(-1, 1, -1);
        Vector3D p2 = Vector3D::point(1, 1, 1);
        Vector3D p3 = Vector3D::point(-1, -1, 1);
        Vector3D p4 = Vector3D::point(1, 1, -1);
        Vector3D p5 = Vector3D::point(-1, -1, -1);
        Vector3D p6 = Vector3D::point(1, -1, 1);
        Vector3D p7 = Vector3D::point(-1, 1, 1);
        points = {p0, p1, p2, p3, p4, p5, p6, p7};
    
        Face f0({0, 4, 2, 6});
        Face f1({4, 1, 7, 2});
        Face f2({1, 5, 3, 7});
        Face f3({5, 0, 6, 3});
        Face f4({6, 2, 7, 3});
        Face f5({0, 5, 1, 4});
        faces = {f0, f1, f2, f3, f4, f5};
    }
    
    void Cube::generate_menger(const unsigned int& max_it) {
        for (unsigned int i = 0; i < max_it; i++ ) {
            /* VISUAL REPRESENTATION OF FACE DIVISION
             *
             * Originial face:
             *
             * 3 ------------- 2
             * |               |
             * |               |
             * |               |
             * |               |
             * |               |
             * 0 ------------- 1
             *
             * Divided face (view from front):
             *
             * 15 -- 13 -- 11 -- 10
             * |  6  |  5  |  4  |
             * 14 -- 12 -- 9 --- 8
             * |  7  |     |  3  |
             * 3 --- 2 --- 5 --- 7
             * |  0  |  1  |  2  |
             * 0 --- 1 --- 4 --- 6
             *
             * Divided face (center):
             *
             *     19--18
             *   / |  / |
             * 12 -- 9  |
             * |  16-|-17
             * | /   |/
             * 2 --- 5
             */
            std::vector<Face> divided_faces;
            std::vector<Vector3D> divided_points;
            for (auto& face : faces) {
                Vector3D p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19;
                unsigned int i0, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14, i15, i16, i17, i18, i19;
    
                // Assign corners
                p0 = points[face[0]]; i0 = (unsigned int) divided_points.size();
                p6 = points[face[1]]; i6 = i0 + 6;
                p10 = points[face[2]]; i10 = i0 + 10;
                p15 = points[face[3]]; i15 = i0 + 15;
    
                // Create outer points first
                p1 = p0 + ((p6 - p0) / 3.0); i1 = i0 + 1;
                p4 = p0 + (2.0 * (p6 - p0) / 3.0); i4 = i0 + 4;
                p7 = p6 + ((p10 - p6) / 3.0); i7 = i0 + 7;
                p8 = p6 + (2.0 * (p10 - p6) / 3.0); i8 = i0 + 8;
                p13 = p15 + ((p10 - p15) / 3.0); i13 = i0 + 13;
                p11 = p15 + (2.0 * (p10 - p15) / 3.0); i11 = i0 + 11;
                p3 = p0 + ((p15 - p0) / 3.0); i3 = i0 + 3;
                p14 = p0 + (2.0 * (p15 - p0) / 3.0); i14 = i0 + 14;
    
                // Create inner points
                p2 = p3 + ((p7 - p3) / 3.0); i2 = i0 + 2;
                p5 = p3 + (2.0 * (p7 - p3) / 3.0); i5 = i0 + 5;
                p12 = p14 + ((p8 - p14) / 3.0); i12 = i0 + 12;
                p9 = p14 + (2.0 * (p8 - p14) / 3.0); i9 = i0 + 9;
    
                // Create center points
                p16 = p3.cross_equals(p1); i16 = i0 + 16;
                p17 = p7.cross_equals(p4); i17 = i0 + 17;
                p18 = p11.cross_equals(p8); i18 = i0 + 18;
                p19 = p14.cross_equals(p13); i19 = i0 + 19;
    
                divided_points.insert(
                        divided_points.end(),
                        {p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19});
                // Create faces
                Face f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11;
                f0 = {i0, i1, i2, i3};
                f1 = {i1, i4, i5, i2};
                f2 = {i4, i6, i7, i5};
                f3 = {i5, i7, i8, i9};
                f4 = {i9, i8, i10, i11};
                f5 = {i12, i9, i11, i13};
                f6 = {i14, i12, i13, i15};
                f7 = {i3, i2, i12, i14};
                f8 = {i2, i5, i17, i16};
                f9 = {i5, i17, i18, i9};
                f10 = {i9, i18, i19, i12};
                f11 = {i2, i16, i9, i12};
    
                divided_faces.insert(divided_faces.end(), {f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11});
            }
            faces = divided_faces;
            points = divided_points;
        }
    }
    
    Tetrahedron::Tetrahedron() {
        Vector3D p0 = Vector3D::point(1, -1, -1);
        Vector3D p1 = Vector3D::point(-1, 1, -1);
        Vector3D p2 = Vector3D::point(1, 1, 1);
        Vector3D p3 = Vector3D::point(-1, -1, 1);
        (*this).points = {p0, p1, p2, p3};
    
        Face f0({0, 1, 2});
        Face f1({1, 3, 2});
        Face f2({0, 1, 3});
        Face f3({0, 2, 3});
        (*this).faces = {f0, f1, f2, f3};
    }
    
    Octahedron::Octahedron() {
        Vector3D p0 = Vector3D::point(1, 0, 0);
        Vector3D p1 = Vector3D::point(0, 1, 0);
        Vector3D p2 = Vector3D::point(-1, 0, 0);
        Vector3D p3 = Vector3D::point(0, -1, 0);
        Vector3D p4 = Vector3D::point(0, 0, -1);
        Vector3D p5 = Vector3D::point(0, 0, 1);
        (*this).points = {p0, p1, p2, p3, p4, p5};
        Face f0({0, 1, 5});
        Face f1({1, 2, 5});
        Face f2({2, 3, 5});
        Face f3({3, 0, 5});
        Face f4({1, 0, 4});
        Face f5({2, 1, 4});
        Face f6({3, 2, 4});
        Face f7({0, 3, 4});
        (*this).faces = {f0, f1, f2, f3, f4, f5, f6, f7};
    }
    
    Icosahedron::Icosahedron() {
        Vector3D p0 = Vector3D::point(0, 0, std::sqrt(5) / 2);
        (*this).points.push_back(p0);
        for (int i=2; i<7; i++){
            Vector3D p = Vector3D::point(std::cos((i-2) * (2 * PI / 5)), std::sin((i-2) * (2 * PI / 5)), 0.5);
            (*this).points.push_back(p);
        }
        for (int i=7; i<12; i++){
            Vector3D p = Vector3D::point(
                    std::cos((PI / 5) + ((i-7) * (2 * PI / 5))),
                    std::sin((PI / 5) + ((i-7) * (2 * PI / 5))),
                    -0.5);
            (*this).points.push_back(p);
        }
        Vector3D p12 = Vector3D::point(0, 0, -std::sqrt(5) / 2);
        (*this).points.push_back(p12);
        Face f0({0, 1, 2});
        Face f1({0, 2, 3});
        Face f2({0, 3, 4});
        Face f3({0, 4, 5});
        Face f4({0, 5, 1});
        Face f5({1, 6, 2});
        Face f6({2, 6, 7});
        Face f7({2, 7, 3});
        Face f8({3, 7, 8});
        Face f9({3, 8, 4});
        Face f10({4, 8, 9});
        Face f11({4, 9, 5});
        Face f12({5, 9, 10});
        Face f13({5, 10, 1});
        Face f14({1, 10, 6});
        Face f15({11, 7, 6});
        Face f16({11, 8, 7});
        Face f17({11, 9, 8});
        Face f18({11, 10, 9});
        Face f19({11, 6, 10});
        (*this).faces = {f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17, f18, f19};
    }
    
    Dodecahedron::Dodecahedron() {
        Face f0({0, 1, 2, 3, 4});
        Face f1({0, 5, 6, 7, 1});
        Face f2({1, 7, 8, 9, 2});
        Face f3({2, 9, 10, 11, 3});
        Face f4({3, 11, 12, 13, 4});
        Face f5({4, 13, 14, 5, 0});
        Face f6({19, 18, 17, 16, 15});
        Face f7({19, 14, 13, 12, 18});
        Face f8({18, 12, 11, 10, 17});
        Face f9({17, 10, 9, 8, 16});
        Face f10({16, 8, 7, 6, 15});
        Face f11({15, 6, 7, 14, 19});
        (*this).faces = {f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11};
    }

    Sphere::Sphere(const double& radius, const unsigned int& n) {
        Icosahedron fig;
        // Reserve enough space for faces
        double szFaces = sizeof(fig.faces[0]) * fig.faces.size() * std::pow(3, n);
        fig.faces.reserve((unsigned long) std::round(szFaces));
        // Reserve enough space for points
        double szPoints = sizeof(fig.points[0]) * fig.points.size() * std::pow(3, n);
        fig.points.reserve((unsigned long) std::round(szPoints));
        // Repeat n times
        for (unsigned int j = 0; j < n; j++){
            unsigned int size = (unsigned int) fig.faces.size();
            // Divide every triangle (face) into four smaller triangles
            for (unsigned int i=0; i<size; i++){
                Face& face = fig.faces[i];
                /* Visual aid:
                 * Assume an upright triangle where the peak is the 1st point,
                 * the 2nd point is to the left and the 3rd point is to the right
                 */
                fig.points[face.point_indexes[0]].normalise();
                fig.points[face.point_indexes[1]].normalise();
                fig.points[face.point_indexes[2]].normalise();
                Vector3D& first = fig.points[face.point_indexes[0]];
                Vector3D& second = fig.points[face.point_indexes[1]];
                Vector3D& third = fig.points[face.point_indexes[2]];
                // Then p1 is the middle of the 1st and 2nd points
                Vector3D p1 = (first + second) / 2.0;
                p1.normalise();
                fig.points.push_back(p1);
                unsigned int p1Pos = (unsigned int) fig.points.size() - 1;
                // And p2 is the middle of the 2nd and the 3rd points
                Vector3D p2 = (second + third) / 2.0;
                p2.normalise();
                fig.points.push_back(p2);
                unsigned int p2Pos = p1Pos + 1;
                // And p3 is the middle of the 3rd and the 1st points
                Vector3D p3 = (first + third) / 2.0;
                p3.normalise();
                fig.points.push_back(p3);
                unsigned int p3Pos = p2Pos + 1;
                // Now, three new faces will be made
                // The triangle {first, p1, p3}
                Face face1({face.point_indexes[0], p1Pos, p3Pos});
                fig.faces.push_back(face1);
                // The triangle {p1, second, p2}
                Face face2({p1Pos, face.point_indexes[1], p2Pos});
                fig.faces.push_back(face2);
                // The triangle {p3, p2, third}
                Face face3({p3Pos, p2Pos, face.point_indexes[2]});
                fig.faces.push_back(face3);
                // The middle triangle {p1, p2, p3}
                Face face4({p1Pos, p2Pos, p3Pos});
                face = face4;
            }
        }
        return fig;
    }

    Cone::Cone(const double& height, const unsigned int& n) {
        for (unsigned int i=0; i<n; i++){
            Vector3D p = Vector3D::point(std::cos((2 * i * PI) / n),
                                         std::sin((2 * i * PI) / n),
                                         0);
            (*this).points.push_back(p);
        }
        Vector3D pN = Vector3D::point(0, 0, height);
        (*this).points.push_back(pN);
        for (unsigned int i = 0; i < n; i++){
            Face f({i, (i + 1) % n, (unsigned int) ((*this).points.size() - 1)});
            (*this).faces.push_back(f);
        }
        Face fN;
        for (unsigned int i = (unsigned int) (*this).points.size()-1; i>0; i--){
            fN.point_indexes.push_back(i);
        }
    }

    Cylinder::Cylinder(const double& height, const unsigned int& n) {
        // Base points
        for (unsigned int i=0; i<n; i++) {
            Vector3D p = Vector3D::point(std::cos((2 * i * PI) / n),
                                         std::sin((2 * i * PI) / n),
                                         0);
            (*this).points.push_back(p);
        }
        // Base face
        Face fBase;
        for (int i = (unsigned int) fig.points.size() - 1; i >= 0; i--){
            fBase.point_indexes.push_back((unsigned int) i);
        }
        (*this).faces.push_back(fBase);
        // Top points
        for (unsigned int i = n; i < 2*n; i++){
            Vector3D p = Vector3D::point(std::cos((2 * i * PI) / n),
                                         std::sin((2 * i * PI) / n),
                                         height);
            (*this).points.push_back(p);
        }
        // Side faces
        for (unsigned int i=0; i < n; i++){
            Face f;
            f.point_indexes = {i, (i + 1) % n, ((i + 1) % n) + n, i + n};
            (*this).faces.push_back(f);
        }
        // Top face
        Face fTop;
//     for (unsigned int i = (unsigned int) fig.points.size() - 1; i >= n; i--){
//         fTop.point_indexes.push_back(i);
//     }
        for (unsigned int i = n; i < (*this).points.size(); i++){
            fTop.point_indexes.push_back(i);
        }
        (*this).faces.push_back(fTop);
    }

    Torus::Torus(const double& r, const double& R, const unsigned int& m, const unsigned int& n) {
        for (int i=0; i<n; i++){
            double u = (2*i*PI) / n;
            for (int j=0; j<m; j++){
                double v = (2*j*PI) / m;
                Vector3D p1 = Vector3D::point(((R + std::cos(v))*std::cos(u)),
                                              ((R + std::cos(v))*std::sin(u)),
                                              r*std::sin(v));
                (*this).points.push_back(p1);
            }
        }
        unsigned int size = (unsigned int) (*this).points.size();
        for (int k=0; k<size; k++){
            unsigned int p1 = k % size;
            unsigned int p2 = (k + m) % size;
            unsigned int p3 = ((k + m + 1)) % size;
            unsigned int p4 = ((k + 1)) % size;
            if (k % m == m - 1) {
                p3 = ((k + m + 1) - m) % size;
                p4 = ((k + 1) - m) % size;
            }
            Face f({p1, p2, p3, p4});
            (*this).faces.push_back(f);
        }
    }

    Buckyball::Buckyball() {
        Icosahedron fig;
        /* VISUAL REPRESENTATION OF FACE DIVISION
         *
         * Original face:
         *
         *      2
         *    /  \
         *   /    \
         *  /      \
         * 0 ------ 1
         *
         * Divided face:
         *
         *        C
         *      / 2\
         *     4 -- 3
         *    /      \
         *   5   3   2
         *  /0\     /1\
         * A - 0---1 - B
         */
        std::vector<Vector3D> divided_points;
        std::vector<Face> divided_faces;
        // Bucky ball has twelve pentagons, so twelve empty lists
        for (auto& face : fig.faces) {
            Vector3D pA, pB, pC, p0, p1, p2, p3, p4, p5;
            unsigned int i0, i1, i2, i3, i4, i5;

            // Assign corners
            pA = fig.points[face.point_indexes[0]];
            pB = fig.points[face.point_indexes[1]];
            pC = fig.points[face.point_indexes[2]];

            // Create outer points first
            p0 = pA + ((pB - pA) / 3.0); i0 = (unsigned int) divided_points.size();
            p1 = pA + (2.0 * (pB - pA) / 3.0); i1 = i0 + 1;
            p2 = pB + ((pC - pB) / 3.0); i2 = i0 + 2;
            p3 = pB + (2.0 * (pC - pB) / 3.0); i3 = i0 + 3;
            p4 = pA + (2.0 * (pC - pA) / 3.0); i4 = i0 + 4;
            p5 = pA + ((pC - pA) / 3.0); i5 = i0 + 5;
            divided_points.insert(divided_points.end(), {p0, p1, p2, p3, p4, p5});

            // Create faces
            Face f({i0, i1, i2, i3, i4, i5});
            divided_faces.push_back(f);
        }
        (*this).points = divided_points;
        (*this).faces = divided_faces;
//        std::vector<unsigned int> p0 = {0, 2, 6, 10,14,18,22,  30,36,42,48,53};
//        std::vector<unsigned int> p1 = {5, 1, 4, 8, 12,16,21,  28,34,40,46,58};
//        std::vector<unsigned int> p2 = {9, 19,3, 7, 11,15,49,  26,32,38,44,56};
//        std::vector<unsigned int> p3 = {13,47,23,29,35,41,52,  25,31,37,43,54};
//        std::vector<unsigned int> p4 = {17,20,27,33,39,45,24,  51,55,57,59,50};
//        for(int i = 0; i < 12; i++){
//            Face newFace = Face();
//            newFace.point_indexes = {p0[i],p1[i],p2[i],p3[i],p4[i]};
//            fig.faces.push_back(newFace);
//        }
    }

    void recursivePrintString(
            LParser::LSystem2D& l_system, std::string& print_string, unsigned int& cur_it, unsigned int& max_it,
            double& alph_angle, double& delt_angle, Lines2D& lines, std::stack<double>& brStack,
            Color& lc, Color& lcGr, bool isGrL, double& x, double& y){
        for (char &c : print_string){
            // Iterate over characters of string
            if (cur_it < max_it) {
                // Recursion ahead!
                if (c=='+') alph_angle += delt_angle;
                else if (c=='-') alph_angle -= delt_angle;
                else if (c=='(') {
                    brStack.push(x);
                    brStack.push(y);
                    brStack.push(alph_angle);
                } else if (c==')') {
                    alph_angle = brStack.top();
                    brStack.pop();
                    y = brStack.top();
                    brStack.pop();
                    x = brStack.top();
                    brStack.pop();
                } else {
                    std::string repl = l_system.get_replacement(c);
                    recursivePrintString(l_system, repl, ++cur_it, max_it,
                                         alph_angle, delt_angle,
                                         lines, brStack, lc, lcGr, isGrL, x, y);
                }
            }
            else {
                // You have reached max iterations
                if (c=='+') alph_angle += delt_angle;
                else if (c=='-') alph_angle -= delt_angle;
                else if (c=='(') {
                    brStack.push(x);
                    brStack.push(y);
                    brStack.push(alph_angle);
                } else if (c==')') {
                    alph_angle = brStack.top();
                    brStack.pop();
                    y = brStack.top();
                    brStack.pop();
                    x = brStack.top();
                    brStack.pop();
                } else if (l_system.draw(c)) {
                    Point2D p1 (x, y);
                    x += std::cos(alph_angle);
                    y += std::sin(alph_angle);
                    Point2D p2 (x, y);
                    Line2D line(p1, p2, lc, lcGr, isGrL);
                    lines.push_back(line);
                } else {
                    x += std::cos(alph_angle);
                    y += std::sin(alph_angle);
                }
            }
        }
        --cur_it;
    }

    void recursivePrintString(
            LParser::LSystem3D& l_system, std::string& print_string, unsigned int& cur_it, unsigned int& max_it,
            double& delt_angle, Figure& figure, std::stack<double>& doubleStack, std::stack<Vector3D>& vectorStack,
            double& x, double& y, double& z, Vector3D& H, Vector3D& L, Vector3D& U){
        for (char &c : print_string){
            // Iterate over characters of string
            if (cur_it < max_it) {
                // Recursion ahead!
                if (c=='+') {
                    Vector3D H_New; H_New =
                                            H * std::cos(delt_angle) + L * std::sin(delt_angle);
                    Vector3D L_New; L_New =
                                            -H * std::sin(delt_angle) + L * std::cos(delt_angle);
                    H = H_New;
                    L = L_New;
                } else if (c=='-') {
                    Vector3D H_New; H_New =
                                            H * std::cos(-delt_angle) + L * std::sin(-delt_angle);
                    Vector3D L_New; L_New =
                                            -H * std::sin(-delt_angle) + L * std::cos(-delt_angle);
                    H = H_New;
                    L = L_New;
                } else if (c == '^') {
                    Vector3D H_New; H_New =
                                            H * std::cos(delt_angle) + U * std::sin(delt_angle);
                    Vector3D U_New; U_New =
                                            -H * std::sin(delt_angle) + U * std::cos(delt_angle);
                    H = H_New;
                    U = U_New;
                } else if (c == '&') {
                    Vector3D H_New; H_New =
                                            H * std::cos(-delt_angle) + U * std::sin(-delt_angle);
                    Vector3D U_New; U_New =
                                            -H * std::sin(-delt_angle) + U * std::cos(-delt_angle);
                    H = H_New;
                    U = U_New;
                } else if (c == '\\') {
                    Vector3D L_New; L_New =
                                            L * std::cos(delt_angle) - U * std::sin(delt_angle);
                    Vector3D U_New; U_New =
                                            L * std::sin(delt_angle) + U * std::cos(delt_angle);
                    L = L_New;
                    U = U_New;
                } else if (c == '/') {
                    Vector3D L_New; L_New =
                                            L * std::cos(-delt_angle) - U * std::sin(-delt_angle);
                    Vector3D U_New; U_New =
                                            L * std::sin(-delt_angle) + U * std::cos(-delt_angle);
                    L = L_New;
                    U = U_New;
                } else if (c == '|') {
                    H = -H;
                    L = -L;
                } else if (c=='(') {
                    doubleStack.push(x);
                    doubleStack.push(y);
                    doubleStack.push(z);
                    doubleStack.push(delt_angle);
                    vectorStack.push(H);
                    vectorStack.push(L);
                    vectorStack.push(U);
                } else if (c==')') {
                    delt_angle = doubleStack.top();
                    doubleStack.pop();
                    z = doubleStack.top();
                    doubleStack.pop();
                    y = doubleStack.top();
                    doubleStack.pop();
                    x = doubleStack.top();
                    doubleStack.pop();
                    U = vectorStack.top();
                    vectorStack.pop();
                    L = vectorStack.top();
                    vectorStack.pop();
                    H = vectorStack.top();
                    vectorStack.pop();
                } else {
                    std::string repl = l_system.get_replacement(c);
                    recursivePrintString(l_system, repl, ++cur_it, max_it,
                                         delt_angle,
                                         figure, doubleStack, vectorStack,
                                         x, y, z, H, L, U);
                }
            }
            else {
                // You have reached max iterations
                if (c=='+') {
                    Vector3D H_New; H_New =
                                            H * std::cos(delt_angle) + L * std::sin(delt_angle);
                    Vector3D L_New; L_New =
                                            -H * std::sin(delt_angle) + L * std::cos(delt_angle);
                    H = H_New;
                    L = L_New;
                } else if (c=='-') {
                    Vector3D H_New; H_New =
                                            H * std::cos(-delt_angle) + L * std::sin(-delt_angle);
                    Vector3D L_New; L_New =
                                            -H * std::sin(-delt_angle) + L * std::cos(-delt_angle);
                    H = H_New;
                    L = L_New;
                } else if (c == '^') {
                    Vector3D H_New; H_New =
                                            H * std::cos(delt_angle) + U * std::sin(delt_angle);
                    Vector3D U_New; U_New =
                                            -H * std::sin(delt_angle) + U * std::cos(delt_angle);
                    H = H_New;
                    U = U_New;
                } else if (c == '&') {
                    Vector3D H_New; H_New =
                                            H * std::cos(-delt_angle) + U * std::sin(-delt_angle);
                    Vector3D U_New; U_New =
                                            -H * std::sin(-delt_angle) + U * std::cos(-delt_angle);
                    H = H_New;
                    U = U_New;
                } else if (c == '\\') {
                    Vector3D L_New; L_New =
                                            L * std::cos(delt_angle) - U * std::sin(delt_angle);
                    Vector3D U_New; U_New =
                                            L * std::sin(delt_angle) + U * std::cos(delt_angle);
                    L = L_New;
                    U = U_New;
                } else if (c == '/') {
                    Vector3D L_New; L_New =
                                            L * std::cos(-delt_angle) - U * std::sin(-delt_angle);
                    Vector3D U_New; U_New =
                                            L * std::sin(-delt_angle) + U * std::cos(-delt_angle);
                    L = L_New;
                    U = U_New;
                } else if (c == '|') {
                    H = -H;
                    L = -L;
                } else if (c=='(') {
                    doubleStack.push(x);
                    doubleStack.push(y);
                    doubleStack.push(z);
                    doubleStack.push(delt_angle);
                    vectorStack.push(H);
                    vectorStack.push(L);
                    vectorStack.push(U);
                } else if (c==')') {
                    delt_angle = doubleStack.top();
                    doubleStack.pop();
                    z = doubleStack.top();
                    doubleStack.pop();
                    y = doubleStack.top();
                    doubleStack.pop();
                    x = doubleStack.top();
                    doubleStack.pop();
                    U = vectorStack.top();
                    vectorStack.pop();
                    L = vectorStack.top();
                    vectorStack.pop();
                    H = vectorStack.top();
                    vectorStack.pop();
                } else if (l_system.draw(c)) {
                    Vector3D p1 = Vector3D::point(x, y, z);
                    figure.points.push_back(p1);
                    unsigned int p1Ind = figure.points.size() - 1;
                    x += H.x;
                    y += H.y;
                    z += H.z;
                    Vector3D p2 = Vector3D::point(x, y, z);
                    figure.points.push_back(p2);
                    unsigned int p2Ind = p1Ind + 1;
                    Face f({p1Ind, p2Ind});
                    figure.faces.push_back(f);
                } else {
                    x += H.x;
                    y += H.y;
                    z += H.z;
                }
            }
        }
        --cur_it;
    }
}

Matrix scaleFigure(const double& scale){
    Matrix m;
    m(1, 1) = scale;
    m(2, 2) = scale;
    m(3, 3) = scale;
    return m;
}

Matrix rotateX(const double& angle){
    Matrix m;
    m(2, 2) = std::cos(angle);
    m(2, 3) = std::sin(angle);
    m(3, 2) = -std::sin(angle);
    m(3, 3) = std::cos(angle);
    return m;
}

Matrix rotateY(const double& angle){
    Matrix m;
    m(1, 1) = std::cos(angle);
    m(1, 3) = -std::sin(angle);
    m(3, 1) = std::sin(angle);
    m(3, 3) = std::cos(angle);
    return m;
}

Matrix rotateZ(const double& angle){
    Matrix m;
    m(1, 1) = std::cos(angle);
    m(1, 2) = std::sin(angle);
    m(2, 1) = -std::sin(angle);
    m(2, 2) = std::cos(angle);
    return m;
}

Matrix translate(const Vector3D &vector){
    Matrix m;
    m(4, 1) = vector.x;
    m(4, 2) = vector.y;
    m(4, 3) = vector.z;
    return m;
}

void toPolar(const Vector3D &point, double &r, double &theta, double &phi){
    // Original coordinates
    double x = point.x;
    double y = point.y;
    double z = point.z;
    // Distance to eye
    r = std::sqrt(x*x + y*y + z*z);
    // Angles
    theta = std::atan2(y, x);
    phi = std::acos(z/r);
}

Matrix eyePointTrans(const Vector3D &eyepoint){
    double r;
    double theta;
    double phi;
    toPolar(eyepoint, r, theta, phi);
    Matrix m;
    m(1, 1) = -std::sin(theta);
    m(1, 2) = -std::cos(theta) * std::cos(phi);
    m(1, 3) = std::cos(theta) * std::sin(phi);
    m(2, 1) = std::cos(theta);
    m(2, 2) = -std::sin(theta) * std::cos(phi);
    m(2, 3) = std::sin(theta) * std::sin(phi);
    m(3, 2) = std::sin(phi);
    m(3, 3) = std::cos(phi);
    m(4, 3) = -r;
    return m;
}

/* Figure createLineDrawing(const ini::Configuration &configuration, unsigned int i){
    Figure fig;
    ini::Section figure = configuration["Figure" + std::to_string(i)];
    int nrPoints = figure["nrPoints"].as_int_or_die();
    int nrLines = figure["nrLines"].as_int_or_die();

    std::vector<Vector3D> points;
    std::vector<Face> faces;

    for (int j=0; j<nrPoints; j++){
        std::vector<double> p = figure["point" + std::to_string(j)].as_double_tuple_or_die();
        Vector3D point = Vector3D::point(p[0], p[1], p[2]);
        points.push_back(point);
    }

    for (int j=0; j<nrLines; j++){
        Face f;
        std::vector<int> p = figure["line" + std::to_string(j)].as_int_tuple_or_die();
        for (auto& k : p){
            f.point_indexes.push_back((unsigned int) k);
        }
        faces.push_back(f);
    }
    fig.points = points;
    fig.faces = faces;
    return fig;
} */

double X(const Vector3D &p1, const Vector3D &p2, const Vector3D &p3){
    return (p1.x + p2.x + p3.x) / 3;
}

double Y(const Vector3D &p1, const Vector3D &p2, const Vector3D &p3){
    return (p1.y + p2.y + p3.y) / 3;
}

double Z(const Vector3D &p1, const Vector3D &p2, const Vector3D &p3){
    return (p1.z + p2.z + p3.z) / 3;
}

fig::Figure createDodecahedron(){
    fig::Icosahedron fig;
    std::vector<Vector3D> t = fig.points;
    Vector3D p1 = Vector3D::point(X(t[0], t[1], t[2]), Y(t[0], t[1], t[2]), Z(t[0], t[1], t[2]));
    Vector3D p2 = Vector3D::point(X(t[0], t[2], t[3]), Y(t[0], t[2], t[3]), Z(t[0], t[2], t[3]));
    Vector3D p3 = Vector3D::point(X(t[0], t[3], t[4]), Y(t[0], t[3], t[4]), Z(t[0], t[3], t[4]));
    Vector3D p4 = Vector3D::point(X(t[0], t[4], t[5]), Y(t[0], t[4], t[5]), Z(t[0], t[4], t[5]));
    Vector3D p5 = Vector3D::point(X(t[0], t[5], t[1]), Y(t[0], t[5], t[1]), Z(t[0], t[5], t[1]));
    Vector3D p6 = Vector3D::point(X(t[1], t[6], t[2]), Y(t[1], t[6], t[2]), Z(t[1], t[6], t[2]));
    Vector3D p7 = Vector3D::point(X(t[2], t[6], t[7]), Y(t[2], t[6], t[7]), Z(t[2], t[6], t[7]));
    Vector3D p8 = Vector3D::point(X(t[2], t[7], t[3]), Y(t[2], t[7], t[3]), Z(t[2], t[7], t[3]));
    Vector3D p9 = Vector3D::point(X(t[3], t[7], t[8]), Y(t[3], t[7], t[8]), Z(t[3], t[7], t[8]));
    Vector3D p10 = Vector3D::point(X(t[3], t[8], t[4]), Y(t[3], t[8], t[4]), Z(t[3], t[8], t[4]));
    Vector3D p11 = Vector3D::point(X(t[4], t[8], t[9]), Y(t[4], t[8], t[9]), Z(t[4], t[8], t[9]));
    Vector3D p12 = Vector3D::point(X(t[4], t[9], t[5]), Y(t[4], t[9], t[5]), Z(t[4], t[9], t[5]));
    Vector3D p13 = Vector3D::point(X(t[5], t[9], t[10]), Y(t[5], t[9], t[10]), Z(t[5], t[9], t[10]));
    Vector3D p14 = Vector3D::point(X(t[5], t[10], t[1]), Y(t[5], t[10], t[1]), Z(t[5], t[10], t[1]));
    Vector3D p15 = Vector3D::point(X(t[1], t[10], t[6]), Y(t[1], t[10], t[6]), Z(t[1], t[10], t[6]));
    Vector3D p16 = Vector3D::point(X(t[11], t[7], t[6]), Y(t[11], t[7], t[6]), Z(t[11], t[7], t[6]));
    Vector3D p17 = Vector3D::point(X(t[11], t[8], t[7]), Y(t[11], t[8], t[7]), Z(t[11], t[8], t[7]));
    Vector3D p18 = Vector3D::point(X(t[11], t[9], t[8]), Y(t[11], t[9], t[8]), Z(t[11], t[9], t[8]));
    Vector3D p19 = Vector3D::point(X(t[11], t[10], t[9]), Y(t[11], t[10], t[9]), Z(t[11], t[10], t[9]));
    Vector3D p20 = Vector3D::point(X(t[11], t[6], t[10]), Y(t[11], t[6], t[10]), Z(t[11], t[6], t[10]));
    fig.points = {p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20};
    fig::Face f1({0, 1, 2, 3, 4});
    fig::Face f2({0, 5, 6, 7, 1});
    fig::Face f3({1, 7, 8, 9, 2});
    fig::Face f4({2, 9, 10, 11, 3});
    fig::Face f5({3, 11, 12, 13, 4});
    fig::Face f6({4, 13, 14, 5, 0});
    fig::Face f7({19, 18, 17, 16, 15});
    fig::Face f8({19, 14, 13, 12, 18});
    fig::Face f9({18, 12, 11, 10, 17});
    fig::Face f10({17, 10, 9, 8, 16});
    fig::Face f11({16, 8, 7, 6, 15});
    fig::Face f12({15, 6, 7, 14, 19});
    fig.faces = {f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12};
    return fig;
}

/* Figure createSpecialSphere(const double &r, const double &R, const unsigned int &m, const unsigned int &n) {
    Figure fig;
    for (int i=0; i<n; i++){
        double u = (2*i*PI) / n;
        for (int j=0; j<m; j++){
            double v = (2*j*PI) / m;
            Vector3D p1 = Vector3D::point((R + (std::cos(v)*std::cos(u))),
                                          (R + (std::cos(v)*std::sin(u))),
                                          r*std::sin(v));
            fig.points.push_back(p1);
        }
    }
    unsigned int size = (unsigned int) fig.points.size();
    for (int k=0; k<n*m; k++){
        unsigned int p1 = k % size;
        unsigned int p2 = (k + m) % size;
        unsigned int p3;
        unsigned int p4;
        if (k % m == m - 1) {
            p3 = ((k + m + 1) - m) % size;
            p4 = ((k + 1) - m) % size;
        } else {
            p3 = ((k + m + 1)) % size;
            p4 = ((k + 1)) % size;
        }
        Face f({p1, p2, p3, p4});
        fig.faces.push_back(f);
    }
    return fig;
} */

//void triangulate(Figure& figure) {
//    std::vector<Face> triang_faces;
//    for (auto& face : figure.faces) {
//        unsigned int i = 1; // Start at second index
//        unsigned int last_index = (unsigned int) face.point_indexes.size() - 2; // We only need to get to the second last element
//        for (i; i<=last_index; i++) {
//            Face f({face.point_indexes[0], face.point_indexes[i], face.point_indexes[i+1]});
//            triang_faces.push_back(f);
//        }
//    }
//    figure.faces = triang_faces;
//}

//void recursiveGenerateFractal(
//        Figure& figure, Figures3D& figures, const double& scale, unsigned int cur_it, const unsigned int & max_it) {
//    if (cur_it < max_it) {
//        // Recursion ahead!
//        for (unsigned int i = 0; i < figure.points.size(); i++) {
//            Figure fractal = figure;
//            // Rescale fractal
//            Matrix sclM = scaleFigure(1.0 / scale);
//            fractal.apply_transformation(sclM);
//            // Move fractal
//            Vector3D move; move = figure.points[i] - fractal.points[i];
//            Matrix trsM = translate(move);
//            fractal.apply_transformation(trsM);
//            recursiveGenerateFractal(fractal, figures, scale, cur_it + 1, max_it);
//        }
//    } else if (cur_it == max_it) {
//        // We have reached max depth!
//        for (unsigned int i = 0; i < figure.points.size(); i++) {
//            Figure fractal = figure;
//            // Rescale fractal
//            Matrix sclM = scaleFigure(1.0 / scale);
//            fractal.apply_transformation(sclM);
//            // Move fractal
//            Vector3D move; move = figure.points[i] - fractal.points[i];
//            Matrix trsM = translate(move);
//            fractal.apply_transformation(trsM);
//            // Add the fractal!
//            figures.push_back(&fractal);
//        }
//    }
//}

//void generateMengerSponge(Figure& figure, const unsigned int& max_it) {
//    for (unsigned int i = 0; i < max_it; i++ ) {
//        /* VISUAL REPRESENTATION OF FACE DIVISION
//         *
//         * Originial face:
//         *
//         * 3 ------------- 2
//         * |               |
//         * |               |
//         * |               |
//         * |               |
//         * |               |
//         * 0 ------------- 1
//         *
//         * Divided face (view from front):
//         *
//         * 15 -- 13 -- 11 -- 10
//         * |  6  |  5  |  4  |
//         * 14 -- 12 -- 9 --- 8
//         * |  7  |     |  3  |
//         * 3 --- 2 --- 5 --- 7
//         * |  0  |  1  |  2  |
//         * 0 --- 1 --- 4 --- 6
//         *
//         * Divided face (center):
//         *
//         *     19--18
//         *   / |  / |
//         * 12 -- 9  |
//         * |  16-|-17
//         * | /   |/
//         * 2 --- 5
//         */
//        std::vector<Face> divided_faces;
//        std::vector<Vector3D> divided_points;
//        for (auto& face : figure.faces) {
//            Vector3D p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19;
//            unsigned int i0, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14, i15, i16, i17, i18, i19;
//
//            // Assign corners
//            p0 = figure.points[face.point_indexes[0]]; i0 = (unsigned int) divided_points.size();
//            p6 = figure.points[face.point_indexes[1]]; i6 = i0 + 6;
//            p10 = figure.points[face.point_indexes[2]]; i10 = i0 + 10;
//            p15 = figure.points[face.point_indexes[3]]; i15 = i0 + 15;
//
//            // Create outer points first
//            p1 = p0 + ((p6 - p0) / 3.0); i1 = i0 + 1;
//            p4 = p0 + (2.0 * (p6 - p0) / 3.0); i4 = i0 + 4;
//            p7 = p6 + ((p10 - p6) / 3.0); i7 = i0 + 7;
//            p8 = p6 + (2.0 * (p10 - p6) / 3.0); i8 = i0 + 8;
//            p13 = p15 + ((p10 - p15) / 3.0); i13 = i0 + 13;
//            p11 = p15 + (2.0 * (p10 - p15) / 3.0); i11 = i0 + 11;
//            p3 = p0 + ((p15 - p0) / 3.0); i3 = i0 + 3;
//            p14 = p0 + (2.0 * (p15 - p0) / 3.0); i14 = i0 + 14;
//
//            // Create inner points
//            p2 = p3 + ((p7 - p3) / 3.0); i2 = i0 + 2;
//            p5 = p3 + (2.0 * (p7 - p3) / 3.0); i5 = i0 + 5;
//            p12 = p14 + ((p8 - p14) / 3.0); i12 = i0 + 12;
//            p9 = p14 + (2.0 * (p8 - p14) / 3.0); i9 = i0 + 9;
//
//            // Create center points
////            p16 = (p3 - p2).cross_equals((p1 - p2)); i16 = i0 + 16;
////            p17 = (p7 - p5).cross_equals((p4 - p5)); i17 = i0 + 17;
////            p18 = (p11 - p9).cross_equals((p8 - p9)); i18 = i0 + 18;
////            p19 = (p14 - p12).cross_equals((p13 - p12)); i19 = i0 + 19;
//
//            p16 = p3.cross_equals(p1); i16 = i0 + 16;
//            p17 = p7.cross_equals(p4); i17 = i0 + 17;
//            p18 = p11.cross_equals(p8); i18 = i0 + 18;
//            p19 = p14.cross_equals(p13); i19 = i0 + 19;
//
//            divided_points.insert(
//                    divided_points.end(),
//                    {p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19});
//            // Create faces
//            Face f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11;
//            f0 = {i0, i1, i2, i3};
//            f1 = {i1, i4, i5, i2};
//            f2 = {i4, i6, i7, i5};
//            f3 = {i5, i7, i8, i9};
//            f4 = {i9, i8, i10, i11};
//            f5 = {i12, i9, i11, i13};
//            f6 = {i14, i12, i13, i15};
//            f7 = {i3, i2, i12, i14};
//            f8 = {i2, i5, i17, i16};
//            f9 = {i5, i17, i18, i9};
//            f10 = {i9, i18, i19, i12};
//            f11 = {i2, i16, i9, i12};
//
//            divided_faces.insert(divided_faces.end(), {f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11});
//        }
//        figure.faces = divided_faces;
//        figure.points = divided_points;
//    }
//}
