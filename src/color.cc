//
// Created by sergio on 30/04/17.
//

#include "color.hh"
#include "figures.hh"

namespace col {

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

    Color& Color::operator+=(const Color& rhs) {
        (*this).red = (*this).red + rhs.red <= 1.0 ? (*this).red + rhs.red : 1.0;
        (*this).green = (*this).green + rhs.green <= 1.0 ? (*this).green + rhs.green : 1.0;
        (*this).blue = (*this).blue + rhs.blue <= 1.0 ? (*this).blue + rhs.blue : 1.0;
        return (*this);
    }

    Color operator+(Color lhs, const Color& rhs) {
        return lhs += rhs;
    }

    Color& Color::operator*=(const Color& rhs) {
        (*this).red = (*this).red * rhs.red <= 1.0 ? (*this).red * rhs.red : 1.0;
        (*this).green = (*this).green * rhs.green <= 1.0 ? (*this).green * rhs.green : 1.0;
        (*this).blue = (*this).blue * rhs.blue <= 1.0 ? (*this).blue * rhs.blue : 1.0;
        return (*this);
    }

    Color operator*(Color lhs, const Color& rhs) {
        return lhs *= rhs;
    }

    Color& Color::operator*=(const double& rhs) {
        (*this).red = (*this).red * rhs <= 1.0 ? (*this).red * rhs: 1.0;
        (*this).green = (*this).green * rhs <= 1.0 ? (*this).green * rhs : 1.0;
        (*this).blue = (*this).blue * rhs <= 1.0 ? (*this).blue * rhs : 1.0;
        return (*this);
    }

    Color operator*(Color lhs, const double& rhs) {
        return lhs *= rhs;
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

    ZBuffer::ZBuffer() {}

    ZBuffer::ZBuffer(const unsigned int& width, const unsigned int& height) {
        std::vector<double> h(height, inf);
        std::vector<std::vector<double>> w(width, h);
        z_buf = w;
    }

    double& ZBuffer::operator()(unsigned int x, unsigned int y) {
        return z_buf[x][y];
    }

    OpacityMatrix::OpacityMatrix() {}

    OpacityMatrix::OpacityMatrix(const unsigned int& width, const unsigned int& height) {
        std::vector<double> h(height, 1);
        std::vector<std::vector<double>> w(width, h);
        matrix = w;
    }

    double& OpacityMatrix::operator()(unsigned int x, unsigned int y) {
        return matrix[x][y];
    }

    Light::Light() {
        (*this).ambient_light = Color(0, 0, 0);
        (*this).diffuse_light = Color(0, 0, 0);
        (*this).specular_light = Color(0, 0, 0);
    }

    Color& Light::ambient() {
        return (*this).ambient_light;
    }

    Color& Light::diffuse() {
        return (*this).diffuse_light;
    }

    Color& Light::specular() {
        return (*this).specular_light;
    }

    void Light::set_mask_size(int) {}

    int Light::get_mask_size() {
        return 0;
    }

    void Light::set_direction(Vector3D& direction) {}

    Vector3D Light::get_direction() {
        return Vector3D();
    }

    void Light::set_location(Vector3D& location) {}

    Vector3D Light::get_location() {
        return Vector3D();
    }

    void Light::set_shadow_mask(ZBuffer& shadow_mask) {}

    ZBuffer Light::get_shadow_mask() {
        return ZBuffer();
    }

    void Light::set_shadow_mask_at(unsigned int x, unsigned int y, double val) {}

    double Light::get_shadow_mask_at(unsigned int x, unsigned int y) { ;
        return 0;
    }

    void Light::set_eye(Matrix& eye) {}

    Matrix Light::get_eye() {
        return Matrix();
    }

    void Light::set_values(double& d, double& dx, double& dy) {}

    void Light::get_values(double& d, double& dx, double& dy) {}

    bool Light::is_ambient() {
        return !((*this).ambient_light == Color(0, 0, 0));
    }

    bool Light::is_diffuse_inf() {
        return !((*this).diffuse_light == Color(0, 0, 0));
    }

    bool Light::is_diffuse_pnt() {
        return !((*this).diffuse_light == Color(0, 0, 0));
    }

    bool Light::is_specular() {
        return !((*this).specular_light == Color(0, 0, 0));
    }

    bool Light::is_shadow() {
        return false;
    }

    void InfLight::set_direction(Vector3D& direction) {
        (*this).direction_vector = direction;
    }

    Vector3D InfLight::get_direction() {
        return (*this).direction_vector;
    }

    bool InfLight::is_diffuse_inf() {
        return true;
    }

    bool InfLight::is_diffuse_pnt() {
        return false;
    }

    void PntLight::set_location(Vector3D& location) {
        (*this).location_vector = location;
    }

    Vector3D PntLight::get_location() {
        return (*this).location_vector;
    }

    void PntLight::set_mask_size(int mask_size) {
        (*this).mask_size = mask_size;
    }

    int PntLight::get_mask_size() {
        return mask_size;
    }

    void PntLight::set_shadow_mask(ZBuffer& shadow_mask) {
        (*this).shadow_mask = shadow_mask;
    }

    ZBuffer PntLight::get_shadow_mask() {
        return shadow_mask;
    }

    void PntLight::set_shadow_mask_at(unsigned int x, unsigned int y, double val) {
        shadow_mask(x, y) = val;
    }

    double PntLight::get_shadow_mask_at(unsigned int x, unsigned int y) {
        return shadow_mask(x, y);
    }

    void PntLight::set_eye(Matrix& eye) {
        (*this).eye = eye;
    }

    Matrix PntLight::get_eye() {
        return eye;
    }

    void PntLight::set_values(double& d, double& dx, double& dy) {
        (*this).d = d;
        (*this).dx = dx;
        (*this).dy = dy;
    }

    void PntLight::get_values(double& d, double& dx, double& dy) {
        d = (*this).d;
        dx = (*this).dx;
        dy = (*this).dy;
    }

    bool PntLight::is_diffuse_inf() {
        return false;
    }

    bool PntLight::is_diffuse_pnt() {
        return true;
    }

    bool PntLight::is_shadow() {
        return mask_size != 0;
    }
}