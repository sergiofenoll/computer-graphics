//
// Created by sergio on 30/04/17.
//

#ifndef COMPGRAPHX_DRAW_HH
#define COMPGRAPHX_DRAW_HH

#include "easy_image.hh"
#include "color.hh"
#include "figures.hh"
#include "projection.hh"

const int TRANSPARENT =  1;
const int OPAQUE = 0;

namespace drw {

    enum DrawType {
        Wires,
        ZBuff,
        Triangles,
        Opaque
    };

    img::EasyImage draw(
            fig::Figures& figures, DrawType& type, col::Lights& lights,
            Matrix& eyeMatrix, const unsigned int& size, col::Color& bc);

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
            Matrix& eyeMatrix, Vector3D& A, Vector3D& B, Vector3D& C,
            const double& d, const double& dx, const double& dy,
            fig::Figure* fig, col::Lights& lights);

    void draw_transparent(
            img::EasyImage& image, col::ZBuffer& z_buffer,
            col::OpacityMatrix& opacity_matrix,
            Matrix& eyeMatrix, Vector3D& A, Vector3D& B, Vector3D& C,
            const double& d, const double& dx, const double& dy,
            fig::Figure* fig, col::Lights& lights);

    void apply_face_light(col::Color& color_values,
                          std::set<col::Light*>& pnt_spec_lights,
                          col::Lights& lights,
                          fig::Figure* fig,
                          Vector3D& n);

    void apply_point_light(col::Color& color_values,
                           col::Color& final_color,
                           std::set<col::Light*>& pnt_spec_lights,
                           fig::Figure* fig, Matrix& eyeMatrix,
                           double z_inv, double z_inv_accurate,
                           unsigned int xI, unsigned int yI,
                           double dx, double dy, double d,
                           Vector3D& n);

    void setup_shadow_mask(Vector3D& A, Vector3D& B, Vector3D& C, col::Light* light);
};


#endif //COMPGRAPHX_DRAW_HH
