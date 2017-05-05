//
// Created by sergio on 30/04/17.
//

#ifndef COMPGRAPHX_DRAW_HH
#define COMPGRAPHX_DRAW_HH


#include "easy_image.hh"
#include "color.hh"
#include "figures.hh"
#include "projection.hh"

namespace drw {

    enum DrawType {
        Wires,
        ZBuff,
        Trian,
    };

    img::EasyImage draw(
            fig::Figures& figures, DrawType& type, col::Lights& lights,
            const unsigned int& size, col::Color& bc);

    img::EasyImage draw(prj::Lines& lines, const unsigned int& size, col::Color& bc);

    void draw_lines(
            img::EasyImage& image,
            unsigned int& x0, unsigned int& y0,
            unsigned int& x1, unsigned int& y1,
            col::Color color);

    void draw_zbuf_lines(
            img::EasyImage& image, col::ZBuffer& z_buffer,
            unsigned int& x0, unsigned int& y0, double& z0,
            unsigned int& x1, unsigned int& y1, double& z1,
            col::Color color);

    void draw_zbuf_triangles(
            img::EasyImage& image, col::ZBuffer& z_buffer,
            Vector3D& A, Vector3D& B, Vector3D& C,
            const double& d, const double& dx, const double& dy,
            col::Color ambientReflection,
            col::Color diffuseReflection,
            col::Color specularReflection,
            double reflectionCoeff, col::Lights& lights);
};


#endif //COMPGRAPHX_DRAW_HH
