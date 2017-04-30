//
// Created by sergio on 30/04/17.
//

#ifndef COMPGRAPHX_DRAW_HH
#define COMPGRAPHX_DRAW_HH


#include "easy_image.hh"
#include "figures.hh"

namespace draw {
    enum FigureType {
        Wires,
        ZBuff,
        Trian,
    };

    img::EasyImage draw2DLines(
            fig::Lines2D& lines, const unsigned int& size, const fig::Color& bc, const fig::Color& bcGr, bool isGrB=false);

    img::EasyImage draw_lines(
            fig::Figures3D& figures, FigureType& type, const unsigned int& size, Lights& lights, const fig::Color& bc, const fig::Color& bcGr, bool isGrB=false);
};


#endif //COMPGRAPHX_DRAW_HH
