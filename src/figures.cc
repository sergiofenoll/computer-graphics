#include "figures.hh"

namespace fig {
    
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
    
    void Figure::set_ambient(const col::Color& ambient) {
        (*this).ambient_reflection = ambient;
    }
    
    col::Color Figure::get_ambient() const {
        return (*this).ambient_reflection;
    }
    
    void Figure::set_diffuse(const col::Color& diffuse) {
        (*this).diffuse_reflection = diffuse;
    }
    
    col::Color Figure::get_diffuse() const {
        return (*this).diffuse_reflection;
    }
    
    void Figure::set_specular(const col::Color& specular) {
        (*this).specular_reflection = specular;
    }
    
    col::Color Figure::get_specular() const {
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

    void Figure::set_colors(const ini::Section& section) {
        section["color"].as_fig_color_if_exists(ambient_reflection);
        section["ambientReflection"].as_fig_color_if_exists(ambient_reflection);
        section["diffuseReflection"].as_fig_color_if_exists(diffuse_reflection);
        section["specularReflection"].as_fig_color_if_exists(specular_reflection);
        section["reflectionCoefficient"].as_double_if_exists(reflection_coefficient);
    }

    void Figure::set_colors(Figure* other_figure) {
        ambient_reflection = other_figure->ambient_reflection;
        diffuse_reflection = other_figure->diffuse_reflection;
        specular_reflection = other_figure->specular_reflection;
        reflection_coefficient = other_figure->reflection_coefficient;
    }

    void Figure::set_opacity(double& opacity) {
        assert(opacity >= 0.0);
        assert(opacity <= 1.0);
        (*this).opacity = opacity;
    }

    double Figure::get_opacity() const {
        return opacity;
    }

    bool Figure::is_opaque() {
        return opacity != 1;
    }

    void Figure::triangulate() {
        std::vector<Face> triang_faces;
        for (auto& face : (*this).faces) {
            unsigned int i = 1; // Start at second index
            unsigned int last_index = (unsigned int) face.point_indexes.size() - 2; // We only need to get to the second last element
            for (i; i<=last_index; i++) {
                Face f({face.point_indexes[0], face.point_indexes[i], face.point_indexes[i+1]});
                triang_faces.push_back(f);
            }
        }
        (*this).faces = triang_faces;
    }

    LineDrawing::LineDrawing(const ini::Configuration &configuration, unsigned int i){
        ini::Section figure = configuration["Figure" + std::to_string(i)];
        int nrPoints = figure["nrPoints"].as_int_or_die();
        int nrLines = figure["nrLines"].as_int_or_die();

        Vectors3D points;
        Faces faces;

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
        (*this).points = points;
        (*this).faces = faces;
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

    Tetrahedron::Tetrahedron() {
        Vector3D p0 = Vector3D::point(1, -1, -1);
        Vector3D p1 = Vector3D::point(-1, 1, -1);
        Vector3D p2 = Vector3D::point(1, 1, 1);
        Vector3D p3 = Vector3D::point(-1, -1, 1);
        (*this).points = {p0, p1, p2, p3};
    
        Face f0({0, 1, 2});
        Face f1({1, 3, 2});
        Face f2({0, 3, 1});
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
        Vector3D p0 = Vector3D::point(0.436339, 0.317019, 0.706011);
        Vector3D p1 = Vector3D::point(-0.166667, 0.512947, 0.706011);
        Vector3D p2 = Vector3D::point(-0.539345, 7.40149e-17, 0.706011);
        Vector3D p3 = Vector3D::point(-0.166667, -0.512947, 0.706011);
        Vector3D p4 = Vector3D::point(0.436339, -0.317019, 0.706011);
        Vector3D p5 = Vector3D::point(0.706011, 0.512947, 0.166667);
        Vector3D p6 = Vector3D::point(0.269672, 0.829966, -0.166667);
        Vector3D p7 = Vector3D::point(-0.269672, 0.829966, 0.166667);
        Vector3D p8 = Vector3D::point(-0.706011, 0.512947, -0.166667);
        Vector3D p9 = Vector3D::point(-0.872678, 1.11022e-16, 0.166667);
        Vector3D p10 = Vector3D::point(-0.706011, -0.512947, -0.166667);
        Vector3D p11 = Vector3D::point(-0.269672, -0.829966, 0.166667);
        Vector3D p12 = Vector3D::point(0.269672, -0.829966, -0.166667);
        Vector3D p13 = Vector3D::point(0.706011, -0.512947, 0.166667);
        Vector3D p14 = Vector3D::point(0.872678, -7.40149e-17, -0.166667);
        Vector3D p15 = Vector3D::point(0.166667, 0.512947, -0.706011);
        Vector3D p16 = Vector3D::point(-0.436339, 0.317019, -0.706011);
        Vector3D p17 = Vector3D::point(-0.436339, -0.317019, -0.706011);
        Vector3D p18 = Vector3D::point(0.166667, -0.512947, -0.706011);
        Vector3D p19 = Vector3D::point(0.539345, -7.40149e-17, -0.706011);
        (*this).points = {p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19};
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
        Face f11({15, 6, 5, 14, 19});
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
                Vector3D& first = fig.points[face.point_indexes[0]];
                Vector3D& second = fig.points[face.point_indexes[1]];
                Vector3D& third = fig.points[face.point_indexes[2]];
                // Then p1 is the middle of the 1st and 2nd points
                Vector3D p1 = (first + second) / 2.0;
                fig.points.push_back(p1);
                unsigned int p1Pos = (unsigned int) fig.points.size() - 1;
                // And p2 is the middle of the 2nd and the 3rd points
                Vector3D p2 = (second + third) / 2.0;
                fig.points.push_back(p2);
                unsigned int p2Pos = p1Pos + 1;
                // And p3 is the middle of the 3rd and the 1st points
                Vector3D p3 = (first + third) / 2.0;
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
        for (auto& point : fig.points) {
            point.normalise();
            point *= radius;
        }
        (*this).points = fig.points;
        (*this).faces = fig.faces;
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

    Cylinder::Cylinder(const double& height, const unsigned int& n, bool top_and_bottom) {
        // Base points
        for (unsigned int i=0; i<n; i++) {
            Vector3D p = Vector3D::point(std::cos((2 * i * PI) / n),
                                         std::sin((2 * i * PI) / n),
                                         0);
            (*this).points.push_back(p);
        }
        // Base face
        if (top_and_bottom) {
            Face fBase;
            for (int i = (unsigned int) (*this).points.size() - 1; i >= 0; i--){
                fBase.point_indexes.push_back((unsigned int) i);
            }
            (*this).faces.push_back(fBase);
        }
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
        if (top_and_bottom) {
            Face fTop;
            for (unsigned int i = n; i < (*this).points.size(); i++){
                fTop.point_indexes.push_back(i);
            }
            (*this).faces.push_back(fTop);
        }
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
        // FIXME: shit's broken yo
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
            if (std::find(divided_points.begin(), divided_points.end(), p0) == divided_points.end()) divided_points.push_back(p0);
            if (std::find(divided_points.begin(), divided_points.end(), p1) == divided_points.end()) divided_points.push_back(p1);
            if (std::find(divided_points.begin(), divided_points.end(), p2) == divided_points.end()) divided_points.push_back(p2);
            if (std::find(divided_points.begin(), divided_points.end(), p3) == divided_points.end()) divided_points.push_back(p3);
            if (std::find(divided_points.begin(), divided_points.end(), p4) == divided_points.end()) divided_points.push_back(p4);
            if (std::find(divided_points.begin(), divided_points.end(), p5) == divided_points.end()) divided_points.push_back(p5);

            // Create faces
            Face f({i0, i1, i2, i3, i4, i5});
            divided_faces.push_back(f);
        }
        int count = 0;
//        for (auto& point : divided_points) {
//            std::cout << "Vector3D p" << count << " = Vector3D::point" << point << std::endl;
//            count++;
//        }
        Vector3D p0 = Vector3D::point(0.333333,        0, 0.912023);
        Vector3D p1 = Vector3D::point(0.666667,        0, 0.706011);
        Vector3D p2 = Vector3D::point(0.769672, 0.317019,      0.5);
        Vector3D p3 = Vector3D::point(0.539345, 0.634038,      0.5);
        Vector3D p4 = Vector3D::point(0.206011, 0.634038, 0.706011);
        Vector3D p5 = Vector3D::point(0.103006, 0.317019, 0.912023);
        Vector3D p6 = Vector3D::point(-0.063661, 0.829966,      0.5);
        Vector3D p7 = Vector3D::point(-0.436339, 0.708876,      0.5);
        Vector3D p8 = Vector3D::point(-0.539345, 0.391857, 0.706011);
        Vector3D p9 = Vector3D::point(-0.269672, 0.195928, 0.912023);
        Vector3D p10 = Vector3D::point(-0.809017, 0.195928,      0.5);
        Vector3D p11 = Vector3D::point(-0.809017, -0.195928,      0.5);
        Vector3D p12 = Vector3D::point(-0.539345, -0.391857, 0.706011);
        Vector3D p13 = Vector3D::point(-0.269672, -0.195928, 0.912023);
        Vector3D p14 = Vector3D::point(-0.436339, -0.708876,      0.5);
        Vector3D p15 = Vector3D::point(-0.063661, -0.829966,      0.5);
        Vector3D p16 = Vector3D::point(0.206011, -0.634038, 0.706011);
        Vector3D p17 = Vector3D::point(0.103006, -0.317019, 0.912023);
        Vector3D p18 = Vector3D::point(0.539345, -0.634038,      0.5);
        Vector3D p19 = Vector3D::point(0.769672, -0.317019,      0.5);
        Vector3D p20 = Vector3D::point(0.936339, 0.195928, 0.166667);
        Vector3D p21 = Vector3D::point(0.872678, 0.391857, -0.166667);
        Vector3D p22 = Vector3D::point( 0.64235, 0.708876, -0.166667);
        Vector3D p23 = Vector3D::point(0.475684, 0.829966, 0.166667);
        Vector3D p24 = Vector3D::point(0.475684, 0.829966, 0.166667);
        Vector3D p25 = Vector3D::point( 0.64235, 0.708876, -0.166667);
        Vector3D p26 = Vector3D::point(0.436339, 0.708876,     -0.5);
        Vector3D p27 = Vector3D::point(0.063661, 0.829966,     -0.5);
        Vector3D p28 = Vector3D::point(-0.103006, 0.951057, -0.166667);
        Vector3D p29 = Vector3D::point(0.103006, 0.951057, 0.166667);
        Vector3D p30 = Vector3D::point(-0.475684, 0.829966, -0.166667);
        Vector3D p31 = Vector3D::point(-0.64235, 0.708876, 0.166667);
        Vector3D p32 = Vector3D::point(-0.64235, 0.708876, 0.166667);
        Vector3D p33 = Vector3D::point(-0.475684, 0.829966, -0.166667);
        Vector3D p34 = Vector3D::point(-0.539345, 0.634038,     -0.5);
        Vector3D p35 = Vector3D::point(-0.769672, 0.317019,     -0.5);
        Vector3D p36 = Vector3D::point(-0.936339, 0.195928, -0.166667);
        Vector3D p37 = Vector3D::point(-0.872678, 0.391857, 0.166667);
        Vector3D p38 = Vector3D::point(-0.936339, -0.195928, -0.1666);
        Vector3D p39 = Vector3D::point(-0.872678, -0.391857, 0.166667);
        Vector3D p40 = Vector3D::point(-0.872678, -0.391857, 0.166667);
        Vector3D p41 = Vector3D::point(-0.936339, -0.195928, -0.1666);
        Vector3D p42 = Vector3D::point(-0.769672, -0.317019,     -0.);
        Vector3D p43 = Vector3D::point(-0.539345, -0.634038,     -0.5);
        Vector3D p44 = Vector3D::point(-0.475684, -0.829966, -0.16666);
        Vector3D p45 = Vector3D::point(-0.64235, -0.708876, 0.166667);
        Vector3D p46 = Vector3D::point(-0.103006, -0.951057, -0.166667);
        Vector3D p47 = Vector3D::point(0.103006, -0.951057, 0.166667);
        Vector3D p48 = Vector3D::point(0.103006, -0.951057, 0.166667);
        Vector3D p49 = Vector3D::point(-0.103006, -0.951057, -0.166667);
        Vector3D p50 = Vector3D::point(0.063661, -0.829966,     -0.5);
        Vector3D p51 = Vector3D::point(0.436339, -0.708876,     -0.5);
        Vector3D p52 = Vector3D::point( 0.64235, -0.708876, -0.166667);
        Vector3D p53 = Vector3D::point(0.475684, -0.829966, 0.166667);
        Vector3D p54 = Vector3D::point(0.872678, -0.391857, -0.166667);
        Vector3D p55 = Vector3D::point(0.936339, -0.195928, 0.166667);
        Vector3D p56 = Vector3D::point(0.936339, -0.195928, 0.166667);
        Vector3D p57 = Vector3D::point(0.872678, -0.391857, -0.166667);
        Vector3D p58 = Vector3D::point(0.809017, -0.195928,     -0.5);
        Vector3D p59 = Vector3D::point(0.809017, 0.195928,     -0.5);
        Vector3D p60 = Vector3D::point(-0.103006, 0.317019, -0.912023);
        Vector3D p61 = Vector3D::point(-0.206011, 0.634038, -0.706011);
        Vector3D p62 = Vector3D::point(0.063661, 0.829966,     -0.5);
        Vector3D p63 = Vector3D::point(0.436339, 0.708876,     -0.5);
        Vector3D p64 = Vector3D::point(0.539345, 0.391857, -0.706011);
        Vector3D p65 = Vector3D::point(0.269672, 0.195928, -0.912023);
        Vector3D p66 = Vector3D::point(-0.333333, 4.08216e-17, -0.912023);
        Vector3D p67 = Vector3D::point(-0.666667, 8.16431e-17, -0.706011);
        Vector3D p68 = Vector3D::point(-0.769672, 0.317019,     -0.5);
        Vector3D p69 = Vector3D::point(-0.539345, 0.634038,     -0.5);
        Vector3D p70 = Vector3D::point(-0.103006, -0.317019, -0.912023);
        Vector3D p71 = Vector3D::point(-0.206011, -0.634038, -0.706011);
        Vector3D p72 = Vector3D::point(-0.539345, -0.634038,     -0.5);
        Vector3D p73 = Vector3D::point(0.269672, -0.195928, -0.912023);
        Vector3D p74 = Vector3D::point(0.539345, -0.391857, -0.706011);
        Vector3D p75 = Vector3D::point(0.436339, -0.708876,     -0.5);
        Vector3D p76 = Vector3D::point(0.063661, -0.829966,     -0.5);
        Vector3D p77 = Vector3D::point(0.809017, 0.195928,     -0.5);
        Vector3D p78 = Vector3D::point(0.809017, -0.195928,     -0.5);
        (*this).points = {p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20,
                          p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, p37, p38, p39,
                          p40, p41, p42, p43, p44, p45, p46, p47, p48, p49, p50, p51, p52, p53, p54, p55, p56, p57, p58,
                          p59, p60, p61, p62, p63, p64, p65, p66, p67, p68, p69, p70, p71, p72, p73, p74, p75, p76, p77,
                          p78};
        Face f0({0, 1, 2, 3, 4, 5});
        Face f1({5, 4, 6, 7, 8, 9});
        Face f2({9, 8, 10, 11, 12, 13});
        Face f3({13, 12, 14, 15, 16, 17});
        Face f4({17, 16, 18, 19, 1, 0});
        (*this).faces = {f0, f1, f2, f3, f4};
    }

    void generate_fractal(
            Figure& figure, Figures& figures, const double& scale, unsigned int cur_it, const unsigned int & max_it) {
        if (cur_it < max_it) {
            // Recursion ahead!
            for (unsigned int i = 0; i < figure.points.size(); i++) {
                Figure* fractal = new Figure(figure);
                // Rescale fractal
                Matrix sclM = scaleFigure(1.0 / scale);
                fractal->apply_transformation(sclM);
                // Move fractal
                Vector3D move; move = figure.points[i] - fractal->points[i];
                Matrix trsM = translate(move);
                fractal->apply_transformation(trsM);
                generate_fractal((*fractal), figures, scale, cur_it + 1, max_it);
            }
        } else if (cur_it == max_it) {
            // We have reached max depth!
            for (unsigned int i = 0; i < figure.points.size(); i++) {
                Figure* fractal = new Figure(figure);
                // Rescale fractal
                Matrix sclM = scaleFigure(1.0 / scale);
                fractal->apply_transformation(sclM);
                // Move fractal
                Vector3D move; move = figure.points[i] - fractal->points[i];
                Matrix trsM = translate(move);
                fractal->apply_transformation(trsM);
                // Add the fractal!
                figures.push_back(fractal);
            }
        } else {
            // Zero iterations
            figures.push_back(&figure);
        }
    }

    void generate_thick(
            Figure& figure, Figures& figures, const double& r, const unsigned int& n, const unsigned int& m) {
        for (auto& point : figure.points) {
            Figure* sphere = new Sphere(r, m);
            Matrix trans = translate(point);
            sphere->apply_transformation(trans);
            sphere->set_colors(&figure);
            sphere->set_opacity(figure.opacity);
            figures.push_back(sphere);
        }
        for (auto& line : figure.faces) {
            int size = line.point_indexes.size();
            for (int i = 0; i < size; i++) {
                Vector3D p1 = figure.points[line[i]];
                Vector3D p2 = figure.points[line[(i + 1) % size]];
                Vector3D p = (p2 - p1);
                double h = p.length() / r;
                Figure* cylinder = new Cylinder(h, n, false);
                Matrix scale = scaleFigure(r);
                cylinder->apply_transformation(scale);
                Vector3D pR = Vector3D::point(p.x, p.y, p.z);
                double placeholder, theta, phi;
                toPolar(pR, placeholder, theta, phi);
                Matrix V = rotateY(phi) * rotateZ(theta) * translate(p1);
                cylinder->apply_transformation(V);
                cylinder->set_colors(&figure);
                cylinder->set_opacity(figure.opacity);
                figures.push_back(cylinder);
            }
        }
    }

    void generate_menger(
            Figure& figure, Figures& figures, const unsigned int& max_iter) {
        Figures cubes = {&figure};
        for (unsigned int it = 0; it < max_iter; it++) {
            Figures temp_cubes = cubes;
            cubes = {};
            for (auto fig : temp_cubes) {
                Vectors3D inner_points = {
                        (fig->points[0] + fig->points[4]) / 2,
                        (fig->points[6] + fig->points[2]) / 2,
                        (fig->points[6] + fig->points[0]) / 2,
                        (fig->points[4] + fig->points[2]) / 2,
                        (fig->points[4] + fig->points[1]) / 2,
                        (fig->points[7] + fig->points[1]) / 2,
                        (fig->points[7] + fig->points[2]) / 2,
                        (fig->points[3] + fig->points[7]) / 2,
                        (fig->points[3] + fig->points[6]) / 2,
                        (fig->points[5] + fig->points[3]) / 2,
                        (fig->points[5] + fig->points[0]) / 2,
                        (fig->points[5] + fig->points[1]) / 2,
                };
                // Only generate one layer of cubes
                if (it == max_iter - 1) generate_fractal((*fig), figures, 3, 1, 1);
                else generate_fractal((*fig), cubes, 3, 1, 1);
                for (unsigned int i = 0; i < inner_points.size(); i++) {
                    Figure* fractal = new Figure(*fig);
                    // Rescale fractal
                    Matrix sclM = scaleFigure(1.0 / 3.0);
                    fractal->apply_transformation(sclM);
                    Vectors3D fractal_inner_points = {
                            (fractal->points[0] + fractal->points[4]) / 2,
                            (fractal->points[6] + fractal->points[2]) / 2,
                            (fractal->points[6] + fractal->points[0]) / 2,
                            (fractal->points[4] + fractal->points[2]) / 2,
                            (fractal->points[4] + fractal->points[1]) / 2,
                            (fractal->points[7] + fractal->points[1]) / 2,
                            (fractal->points[7] + fractal->points[2]) / 2,
                            (fractal->points[3] + fractal->points[7]) / 2,
                            (fractal->points[3] + fractal->points[6]) / 2,
                            (fractal->points[5] + fractal->points[3]) / 2,
                            (fractal->points[5] + fractal->points[0]) / 2,
                            (fractal->points[5] + fractal->points[1]) / 2,
                    };
                    // Move fractal
                    Vector3D move;
                    move = inner_points[i] - fractal_inner_points[i];
                    Matrix trsM = translate(move);
                    fractal->apply_transformation(trsM);
                    // Add the fractal!
                    if (it == max_iter - 1) figures.push_back(fractal);
                    else cubes.push_back(fractal);
                }
            }
        }
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
    // m = rotateZ(-theta - (PI / 2.0)) * rotateX(-phi) * translate(Vector3D::vector(0, 0, -r));

    /*
     * This code causes segfaults, burning homes and XK-Class apocalyptic events.
     * The above line was supposed to fix that, but it didn't.
     * It is left there for posterity.
     * Thanks Andrei for the (almost but not) bugfix!
     */

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
