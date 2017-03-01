#include "3DFigures.hh"

Matrix scaleFigure(const double &scale){
    Matrix m;
    m(1, 1) = scale;
    m(2, 2) = scale;
    m(3, 3) = scale;
    return m;
};

Matrix rotateX(const double &angle){
    Matrix m;
    m(2, 2) = std::cos(angle);
    m(2, 3) = std::sin(angle);
    m(3, 2) = -std::sin(angle);
    m(3, 3) = std::cos(angle);
    return m;
};

Matrix rotateY(const double &angle){
    Matrix m;
    m(1, 1) = std::cos(angle);
    m(1, 3) = -std::sin(angle);
    m(3, 1) = std::sin(angle);
    m(3, 3) = std::cos(angle);
    return m;
};

Matrix rotateZ(const double &angle){
    Matrix m;
    m(1, 1) = std::cos(angle);
    m(1, 2) = std::sin(angle);
    m(2, 1) = -std::sin(angle);
    m(2, 2) = std::cos(angle);
    return m;
};

Matrix translate(const Vector3D &vector){
    Matrix m;
    m(3, 1) = vector.x;
    m(3, 2) = vector.y;
    m(3, 3) = vector.z;
    return m;
};

void applyTransformation(Figure &fig, const Matrix &m){
    for(Vector3D &point : fig.points){
        point *= m;
    }
};

void applyTransformation(Figures3D &figs, const Matrix &m){
    for(Figure &fig : figs){
        applyTransformation(fig, m);
    }
};

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
};

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
};
