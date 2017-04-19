#include "3DFigures.hh"

Matrix scaleFigure(const double &scale){
    Matrix m;
    m(1, 1) = scale;
    m(2, 2) = scale;
    m(3, 3) = scale;
    return m;
}

Matrix rotateX(const double &angle){
    Matrix m;
    m(2, 2) = std::cos(angle);
    m(2, 3) = std::sin(angle);
    m(3, 2) = -std::sin(angle);
    m(3, 3) = std::cos(angle);
    return m;
}

Matrix rotateY(const double &angle){
    Matrix m;
    m(1, 1) = std::cos(angle);
    m(1, 3) = -std::sin(angle);
    m(3, 1) = std::sin(angle);
    m(3, 3) = std::cos(angle);
    return m;
}

Matrix rotateZ(const double &angle){
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

void applyTransformation(Figure &fig, const Matrix &m){
    for(Vector3D &point : fig.points){
        point *= m;
    }
}

void applyTransformation(Figures3D &figs, const Matrix &m){
    for(Figure &fig : figs){
        applyTransformation(fig, m);
    }
}

void toPolar(const Vector3D &point,
                            double &r,
                            double &theta,
                            double &phi){
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

Figure createLineDrawing(const ini::Configuration &configuration, unsigned int i){
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
        std::vector<int> p = figure["line" + std::to_string(j)].as_int_tuple_or_die();
        Face f;
        for (int k : p){
            f.point_indexes.push_back(k);
        }
        faces.push_back(f);
    }
    fig.points = points;
    fig.faces = faces;
    return fig;
}

Figure createCube(){
    Figure fig;
    // Cube points
    Vector3D p1 = Vector3D::point(1, -1, -1);
    Vector3D p2 = Vector3D::point(-1, 1, -1);
    Vector3D p3 = Vector3D::point(1, 1, 1);
    Vector3D p4 = Vector3D::point(-1, -1, 1);
    Vector3D p5 = Vector3D::point(1, 1, -1);
    Vector3D p6 = Vector3D::point(-1, -1, -1);
    Vector3D p7 = Vector3D::point(1, -1, 1);
    Vector3D p8 = Vector3D::point(-1, 1, 1);
    fig.points = {p1, p2, p3, p4, p5, p6, p7, p8};
    // Cube faces
    Face f1({0, 4, 2, 6});
    Face f2({4, 1, 7, 2});
    Face f3({1, 5, 3, 7});
    Face f4({5, 0, 6, 3});
    Face f5({6, 2, 7, 3});
    Face f6({0, 5, 1, 4});
    fig.faces = {f1, f2, f3, f4, f5, f6};
    return fig;
}

Figure createTetrahedron(){
    Figure fig;
    // Tetrahedron points
    Vector3D p1 = Vector3D::point(1, -1, -1);
    Vector3D p2 = Vector3D::point(-1, 1, -1);
    Vector3D p3 = Vector3D::point(1, 1, 1);
    Vector3D p4 = Vector3D::point(-1, -1, 1);
    fig.points = {p1, p2, p3, p4};
    // Tetrahedron faces
    Face f1({1, 2, 3});
    Face f2({2, 4, 3});
    Face f3({1, 2, 4});
    Face f4({1, 3, 4});
    fig.faces = {f1, f2, f3, f4};
    // Lazy fix
    for (Face &f : fig.faces){
        for (int &i : f.point_indexes){
            i--;
        }
    }
    return fig;
}

Figure createOctahedron(){
    Figure fig;
    // Octahedron points
    Vector3D p1 = Vector3D::point(1, 0, 0);
    Vector3D p2 = Vector3D::point(0, 1, 0);
    Vector3D p3 = Vector3D::point(-1, 0, 0);
    Vector3D p4 = Vector3D::point(0, -1, 0);
    Vector3D p5 = Vector3D::point(0, 0, -1);
    Vector3D p6 = Vector3D::point(0, 0, 1);
    fig.points = {p1, p2, p3, p4, p5, p6};
    // Octahedron faces
    Face f1({1, 2, 6});
    Face f2({2, 3, 6});
    Face f3({3, 4, 6});
    Face f4({4, 1, 6});
    Face f5({2, 1, 5});
    Face f6({3, 2, 5});
    Face f7({4, 3, 5});
    Face f8({1, 4, 5});
    fig.faces = {f1, f2, f3, f4, f5, f6, f7, f8};
    // Laxy fix
    for (Face &f : fig.faces){
        for (int &i : f.point_indexes){
            i--;
        }
    }
    return fig;
}

Figure createIcosahedron(){
    Figure fig;
    // Icosahedron points
    Vector3D p1 = Vector3D::point(0, 0, std::sqrt(5) / 2);
    fig.points.push_back(p1);
    for (int i=2; i<7; i++){
        Vector3D p = Vector3D::point(std::cos((i-2) * (2 * PI / 5)), std::sin((i-2) * (2 * PI / 5)), 0.5);
        fig.points.push_back(p);
    }
    for (int i=7; i<12; i++){
        Vector3D p = Vector3D::point(std::cos((PI / 5) + ((i-7) * (2 * PI / 5))), std::sin((PI / 5) + ((i-7) * (2 * PI / 5))), -0.5);
        fig.points.push_back(p);
    }
    Vector3D p12 = Vector3D::point(0, 0, -std::sqrt(5) / 2);
    fig.points.push_back(p12);
    // Icosahedron faces
    Face f1({1, 2, 3});
    Face f2({1, 3, 4});
    Face f3({1, 4, 5});
    Face f4({1, 5, 6});
    Face f5({1, 6, 2});
    Face f6({2, 7, 3});
    Face f7({3, 7, 8});
    Face f8({3, 8, 4});
    Face f9({4, 8, 9});
    Face f10({4, 9, 5});
    Face f11({5, 9, 10});
    Face f12({5, 10, 6});
    Face f13({6, 10, 11});
    Face f14({6, 11, 2});
    Face f15({2, 11, 7});
    Face f16({12, 8, 7});
    Face f17({12, 9, 8});
    Face f18({12, 10, 9});
    Face f19({12, 11, 10});
    Face f20({12, 7, 11});
    fig.faces = {f1, f2, f3, f4, f5, f6, f7, f8, f9, f10,
                 f11, f12, f13, f14, f15, f16, f17, f18, f19, f20};
    // Lazy fix
    for (Face &f : fig.faces){
        for (int &i : f.point_indexes){
            i--;
        }
    }
    return fig;
}

double X(const Vector3D &p1, const Vector3D &p2, const Vector3D &p3){
    return (p1.x + p2.x + p3.x) / 3;
}

double Y(const Vector3D &p1, const Vector3D &p2, const Vector3D &p3){
    return (p1.y + p2.y + p3.y) / 3;
}

double Z(const Vector3D &p1, const Vector3D &p2, const Vector3D &p3){
    return (p1.z + p2.z + p3.z) / 3;
}

Figure createDodecahedron(){
    Figure fig;
    std::vector<Vector3D> t = createIcosahedron().points;
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
    Vector3D p15 = Vector3D::point(X(t[1], t[10], t[1]), Y(t[1], t[10], t[6]), Z(t[1], t[10], t[6]));
    Vector3D p16 = Vector3D::point(X(t[11], t[7], t[6]), Y(t[11], t[7], t[6]), Z(t[11], t[7], t[6]));
    Vector3D p17 = Vector3D::point(X(t[11], t[8], t[7]), Y(t[11], t[8], t[7]), Z(t[11], t[8], t[7]));
    Vector3D p18 = Vector3D::point(X(t[11], t[9], t[8]), Y(t[11], t[9], t[8]), Z(t[11], t[9], t[8]));
    Vector3D p19 = Vector3D::point(X(t[11], t[10], t[9]), Y(t[11], t[10], t[9]), Z(t[11], t[10], t[9]));
    Vector3D p20 = Vector3D::point(X(t[11], t[6], t[10]), Y(t[11], t[6], t[10]), Z(t[11], t[6], t[10]));
    fig.points = {p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20};
    Face f1({1, 2, 3, 4, 5});
    Face f2({1, 6, 7, 8, 2});
    Face f3({2, 8, 9, 10, 3});
    Face f4({3, 10, 11, 12, 4});
    Face f5({4, 12, 13, 14, 5});
    Face f6({5, 14, 15, 6, 1});
    Face f7({20, 19, 18, 17, 16});
    Face f8({20, 15, 14, 13, 19});
    Face f9({19, 13, 12, 11, 18});
    Face f10({18, 11, 10, 9, 17});
    Face f11({17, 9, 8, 7, 16});
    Face f12({16, 7, 6, 15, 20});
    fig.faces = {f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12};
    // Lazy fix
    for (Face &f : fig.faces){
        for (int &i : f.point_indexes){
            i--;
        }
    }
    return fig;
}

Figure createSphere(double &radius, const unsigned int &n){
    Figure fig = createIcosahedron();
    // Reserve enough space for faces
    int szFaces = sizeof(fig.faces[0]) * fig.faces.size() * std::pow(3, n);
    fig.faces.reserve(szFaces);
    // Reserve enough space for points
    int szPoints = sizeof(fig.points[0]) * fig.points.size() * std::pow(3, n);
    fig.points.reserve(szPoints);
    // Repeat n times
    for (unsigned int j=0; j<n; j++){
        unsigned int size = fig.faces.size();
        // Divide every triangle (face) into four smaller triangles
        for (unsigned int i=0; i<size; i++){
            Face &face = fig.faces[i];
            /* Visual aid:
             * Assume an upright triangle where the peak is the 1st point,
             * the 2nd point is to the left and the 3rd point is to the right
             */
             fig.points[face.point_indexes[0]].normalise();
             fig.points[face.point_indexes[1]].normalise();
             fig.points[face.point_indexes[2]].normalise();
            Vector3D &first = fig.points[face.point_indexes[0]];
            Vector3D &second = fig.points[face.point_indexes[1]];
            Vector3D &third = fig.points[face.point_indexes[2]];
            // Then p1 is the middle of the 1st and 2nd points
            Vector3D p1 = Vector3D::point((first.x + second.x) / 2,
                                          (first.y + second.y) / 2,
                                          (first.z + second.z) / 2);
            p1.normalise();
            fig.points.push_back(p1);
            int p1Pos = fig.points.size() - 1;
            // And p2 is the middle of the 2nd and the 3rd points
            Vector3D p2 = Vector3D::point((second.x + third.x) / 2,
                                          (second.y + third.y) / 2,
                                          (second.z + third.z) / 2);
            p2.normalise();
            fig.points.push_back(p2);
            int p2Pos = p1Pos + 1;
            // And p3 is the middle of the 3rd and the 1st points
            Vector3D p3 = Vector3D::point((third.x + first.x) / 2,
                                          (third.y + first.y) / 2,
                                          (third.z + first.z) / 2);
            p3.normalise();
            fig.points.push_back(p3);
            int p3Pos = p2Pos + 1;
            // Now, three new faces will be made
            // The triangle {first, p1, p3}
            Face face1({face.point_indexes[0], p1Pos, p3Pos});
            fig.faces.push_back(face1);
            // The triangle {second, p1, p2}
            Face face2({face.point_indexes[1], p1Pos, p2Pos});
            fig.faces.push_back(face2);
            // The triangle {third, p2, p3}
            Face face3({face.point_indexes[2], p2Pos, p3Pos});
            fig.faces.push_back(face3);
            // The middle triangle {p1, p2, p3}
            Face face4({p1Pos, p2Pos, p3Pos});
            face = face4;
        }
    }
    return fig;
}

Figure createCone(const double &h, const unsigned int &n){
    Figure fig;
    for (unsigned int i=0; i<n; i++){
        Vector3D p = Vector3D::point(std::cos((2 * i * PI) / n),
                                     std::sin((2 * i * PI) / n),
                                     0);
        fig.points.push_back(p);
    }
    Vector3D pN = Vector3D::point(0, 0, h);
    fig.points.push_back(pN);
    for (unsigned int i=0; i<n; i++){
        Face f({i, (i + 1) % n, fig.points.size() - 1});
        fig.faces.push_back(f);
    }
    Face fN;
    for (unsigned int i=fig.points.size()-1; i>0; i--){
        fN.point_indexes.push_back(i);
    }
    return fig;
}

Figure createCylinder(const double &h, const unsigned int &n){
    Figure fig;
    // Base points
    for (unsigned int i=0; i<n; i++){
        Vector3D p = Vector3D::point(std::cos((2 * i * PI) / n),
                                     std::sin((2 * i * PI) / n),
                                     0);
        fig.points.push_back(p);
    }
    // Base face
    Face fBase;
    for (int i=fig.points.size()-1; i>=0; i--){
        fBase.point_indexes.push_back(i);
    }
    fig.faces.push_back(fBase);
    // Top points
    for (unsigned int i=n; i<2*n; i++){
        Vector3D p = Vector3D::point(std::cos((2 * i * PI) / n),
                                     std::sin((2 * i * PI) / n),
                                     h);
        fig.points.push_back(p);
    }
    // Side faces
    for (unsigned int i=0; i<n; i++){
        Face f;
        f.point_indexes = {i, (i + 1) % n, ((i + 1) % n) + n, i + n};
        fig.faces.push_back(f);
    }
    // Top face
    Face fTop;
    for (unsigned int i=fig.points.size()-1; i>=n; i--){
        fTop.point_indexes.push_back(i);
    }
    fig.faces.push_back(fTop);
    return fig;
}

Figure createTorus(const double &r, const double &R, const int &m, const int &n){
    Figure fig;
    for (int i=0; i<n; i++){
        double u = (2*i*PI) / n;
        for (int j=0; j<m; j++){
            double v = (2*j*PI) / m;
            Vector3D p1 = Vector3D::point(((R + std::cos(v))*std::cos(u)),
                                          ((R + std::cos(v))*std::sin(u)),
                                          r*std::sin(v));
            fig.points.push_back(p1);
        }
    }
    int size = fig.points.size();
    for (int k=0; k<size; k++){
        int p1 = k % size;
        int p2 = (k + m) % size;
        int p3 = ((k + m + 1)) % size;
        int p4 = ((k + 1)) % size;
        if (k % m == m - 1) {
            p3 = ((k + m + 1) - m) % size;
            p4 = ((k + 1) - m) % size;
        }
        Face f({p1, p2, p3, p4});
        fig.faces.push_back(f);
    }
    return fig;
}

Figure createSpecialSphere(const double &r, const double &R, const int &m, const int &n) {
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
    int size = fig.points.size();
    for (int k=0; k<n*m; k++){
        int p1 = k % size;
        int p2 = (k + m) % size;
        int p3;
        int p4;
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
}

Figure create3DLSystem() {
    return Figure();
}

void triangulate(Figure& figure) {
    std::vector<Face> triang_faces;
    for (auto& face : figure.faces) {
        int i = 1; // Start at second index
        int last_index = face.point_indexes.size() - 2; // We only need to get to the second last element
        for (i; i<=last_index; i++) {
            Face f({face.point_indexes[0], face.point_indexes[i], face.point_indexes[i+1]});
            triang_faces.push_back(f);
        }
    }
    figure.faces = triang_faces;
}
