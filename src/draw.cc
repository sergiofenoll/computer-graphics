//
// Created by sergio on 30/04/17.
//

#include "draw.hh"

namespace draw {
    img::EasyImage draw2DLines(
            fig::Lines2D& lines, const unsigned int& size, const fig::Color& bc, const fig::Color& bcGr, bool isGrB) {
        double maxX = lines[0].p1.x;
        double maxY = lines[0].p1.y;
        double minX = lines[0].p1.x;
        double minY = lines[0].p1.y;
        // Iterate over all lines looking for xmax and ymax
        for (Line2D &l : lines) {
            if (l.p1.x > maxX) maxX = l.p1.x;
            if (l.p1.x < minX) minX = l.p1.x;
            if (l.p2.x > maxX) maxX = l.p2.x;
            if (l.p2.x < minX) minX = l.p2.x;
            if (l.p1.y > maxY) maxY = l.p1.y;
            if (l.p1.y < minY) minY = l.p1.y;
            if (l.p2.y > maxY) maxY = l.p2.y;
            if (l.p2.y < minY) minY = l.p2.y;
        }
        // Ranges
        double xRange = maxX - minX;
        double yRange = maxY - minY;
        // Size of image (width, height)
        double maxXY = std::max(xRange, yRange);
        double imageX = size * (xRange / maxXY);
        double imageY = size * (yRange / maxXY);
        // Scaling factor d
        double d = 0.95 * (imageX / xRange);
        // Centers
        double DCx = d * ((minX + maxX) / 2.0);
        double DCy = d * ((minY + maxY) / 2.0);
        double dx = (imageX / 2.0) - DCx;
        double dy = (imageY / 2.0) - DCy;
        // Multiply every coordinate with d,
        // add dx and dy
        // and draw lines
        img::EasyImage image((unsigned int) std::round(imageX), (unsigned int) std::round(imageY), bc, bcGr, isGrB);
        for(Line2D &l : lines){
            // Update the coordinates
            l.p1.x *= d;
            l.p1.x += dx;
            l.p2.x *= d;
            l.p2.x += dx;
            l.p1.y *= d;
            l.p1.y += dy;
            l.p2.y *= d;
            l.p2.y += dy;
            // Start drawing the lines
            unsigned int x0 = (unsigned int) std::round(l.p1.x);
            unsigned int y0 = (unsigned int) std::round(l.p1.y);
            unsigned int x1 = (unsigned int) std::round(l.p2.x);
            unsigned int y1 = (unsigned int) std::round(l.p2.y);
            image.draw_line(x0, y0, x1, y1, l.color, l.colorGr, l.isGradient);
        }
        return image;
    }

    img::EasyImage draw3DLines(
            fig::Figures3D& figures, FigureType& type, const unsigned int& size, Lights& lights, const fig::Color& bc, const fig::Color& bcGr, bool isGrB) {
        fig::Lines2D lines = doProjection(figures);
        double maxX = lines[0].p1.x;
        double maxY = lines[0].p1.y;
        double minX = lines[0].p1.x;
        double minY = lines[0].p1.y;
        // Iterate over all lines looking for xmax and ymax
        for (fig::Line2D &l : lines) {
            if (l.p1.x > maxX) maxX = l.p1.x;
            if (l.p1.x < minX) minX = l.p1.x;
            if (l.p2.x > maxX) maxX = l.p2.x;
            if (l.p2.x < minX) minX = l.p2.x;
            if (l.p1.y > maxY) maxY = l.p1.y;
            if (l.p1.y < minY) minY = l.p1.y;
            if (l.p2.y > maxY) maxY = l.p2.y;
            if (l.p2.y < minY) minY = l.p2.y;
        }
        // Ranges
        double xRange = maxX - minX;
        double yRange = maxY - minY;
        // Size of image (width, height)
        double maxXY = std::max(xRange, yRange);
        double imageX = size * (xRange / maxXY);
        double imageY = size * (yRange / maxXY);
        // Scaling factor d
        double d = 0.95 * (imageX / xRange);
        // Centers
        double DCx = d * ((minX + maxX) / 2.0);
        double DCy = d * ((minY + maxY) / 2.0);
        double dx = (imageX / 2.0) - DCx;
        double dy = (imageY / 2.0) - DCy;
        // Multiply every coordinate with d,
        // add dx and dy
        // and draw lines
        img::EasyImage image((unsigned int) imageX, (unsigned int) imageY, bc, bcGr, isGrB);
        ZBuffer z_buffer;
        switch (type) {
            case Wires :
                for(fig::Line2D &l : lines){
                    // Update the coordinates
                    l.p1.x *= d;
                    l.p1.x += dx;
                    l.p2.x *= d;
                    l.p2.x += dx;
                    l.p1.y *= d;
                    l.p1.y += dy;
                    l.p2.y *= d;
                    l.p2.y += dy;
                    // Start drawing the lines
                    unsigned int x0 = (unsigned int) std::round(l.p1.x);
                    unsigned int y0 = (unsigned int) std::round(l.p1.y);
                    unsigned int x1 = (unsigned int) std::round(l.p2.x);
                    unsigned int y1 = (unsigned int) std::round(l.p2.y);
                    image.draw_line(x0, y0, x1, y1, l.color, l.colorGr, l.isGradient);
                }
                break;
            case ZBuff :
                z_buffer = ZBuffer((unsigned int) imageX, (unsigned int) imageY);
                for (fig::Line2D &l : lines) {
                    // Update the coordinates
                    l.p1.x *= d;
                    l.p1.x += dx;
                    l.p2.x *= d;
                    l.p2.x += dx;
                    l.p1.y *= d;
                    l.p1.y += dy;
                    l.p2.y *= d;
                    l.p2.y += dy;
                    // Start drawing the lines
                    unsigned int x0 = (unsigned int) std::round(l.p1.x);
                    unsigned int y0 = (unsigned int) std::round(l.p1.y);
                    unsigned int x1 = (unsigned int) std::round(l.p2.x);
                    unsigned int y1 = (unsigned int) std::round(l.p2.y);
                    image.draw_zbuf_line(z_buffer, x0, y0, l.p1.z, x1, y1, l.p2.z, l.color, l.colorGr, l.isGradient);
                }
                break;
            case Trian :
                z_buffer = ZBuffer((unsigned int) imageX, (unsigned int) imageY);
                // Draw images
                for(auto& figure : figures){
                    triangulate((*figure));
                    for (auto& face : figure->faces) {
                        Vector3D A = figure->points[face[0]];
                        Vector3D B = figure->points[face[1]];
                        Vector3D C = figure->points[face[2]];
                        image.draw_zbuf_triang(z_buffer, A, B, C, d, dx, dy, figure->ambient_reflection,
                                               figure->diffuse_reflection, figure->specular_reflection,
                                               figure->reflection_coefficient, lights, figure->colorGr, figure->isGradient);
                    }
                }
                break;
        }
        return image;
    }
}
