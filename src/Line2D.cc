//
// Created by sergio on 21/02/17.
//

#include "Line2D.hh"

Point2D doProjection(const Vector3D& point, const double& d){
    Point2D p;
    p.x = (d * point.x) / -(point.z);
    p.y = (d * point.y) / -(point.z);
    p.z = point.z;
    return p;
}

Lines2D doProjection(Figures3D& figures){
    Lines2D lines;
    for (Figure &fig : figures){
        for (Face &face : fig.faces){
            Points2D points;
            for (int index : face.point_indexes){
                Point2D p = doProjection(fig.points[index]);
                points.push_back(p);
            }
            for (int i=0; i<points.size(); i++) {
                Line2D line;
                line.p1 = points[i];
                line.p2 = points[(i + 1) % points.size()];
                line.color = fig.color;
                line.colorGr = fig.colorGr;
                line.isGradient = fig.isGradient;
                lines.push_back(line);
            }
        }
    }
    return lines;
}

img::EasyImage draw2DLines(
        Lines2D& lines, const unsigned int& size, const img::Color& bc, const img::Color& bcGr, bool isGrB) {
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
        Figures3D& figures, FigureType& type, const unsigned int& size, const img::Color& bc, const img::Color& bcGr, bool isGrB) {
    Lines2D lines = doProjection(figures);
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
    img::EasyImage image((unsigned int) imageX, (unsigned int) imageY, bc, bcGr, isGrB);
    ZBuffer z_buffer;
    switch (type) {
        case Wires :
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
            break;
        case ZBuff :
            z_buffer = ZBuffer((unsigned int) imageX, (unsigned int) imageY);
            for (Line2D &l : lines) {
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
                triangulate(figure);
                for (auto& face : figure.faces) {
                    Vector3D A = figure.points[face.point_indexes[0]];
                    Vector3D B = figure.points[face.point_indexes[1]];
                    Vector3D C = figure.points[face.point_indexes[2]];
                    image.draw_zbuf_triang(z_buffer, A, B, C, d, dx, dy, figure.color, figure.colorGr, figure.isGradient);
                }
            }
            break;
    }
    return image;
}

void recursivePrintString(
        LParser::LSystem2D& l_system, std::string& print_string, unsigned int& cur_it, unsigned int& max_it,
        double& alph_angle, double& delt_angle, Lines2D& lines, std::stack<double>& brStack,
        img::Color& lc, img::Color& lcGr, bool isGrL, double& x, double& y){
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
                int p1Ind = figure.points.size() - 1;
                x += H.x;
                y += H.y;
                z += H.z;
                Vector3D p2 = Vector3D::point(x, y, z);
                figure.points.push_back(p2);
                int p2Ind = p1Ind + 1;
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
