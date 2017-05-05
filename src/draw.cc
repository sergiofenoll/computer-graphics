//
// Created by sergio on 30/04/17.
//

#include "draw.hh"

namespace drw {

    img::EasyImage draw(
            fig::Figures& figures, DrawType& type, col::Lights& lights,
            const unsigned int& size, col::Color& bc) {
        prj::Lines lines = prj::project_figures(figures);
        double maxX = lines[0].p1.x;
        double maxY = lines[0].p1.y;
        double minX = lines[0].p1.x;
        double minY = lines[0].p1.y;
        for (prj::Line &l : lines) {
            maxX = std::max(std::max(l.p1.x, l.p2.x), maxX);
            minX = std::min(std::min(l.p1.x, l.p2.x), minX);
            maxY = std::max(std::max(l.p1.y, l.p2.y), maxY);
            minY = std::min(std::min(l.p1.y, l.p2.y), minY);
        }
        double xRange = maxX - minX;
        double yRange = maxY - minY;
        double maxXY = std::max(xRange, yRange);
        double imageX = size * (xRange / maxXY);
        double imageY = size * (yRange / maxXY);
        double d = 0.95 * (imageX / xRange);
        double DCx = d * ((minX + maxX) / 2.0);
        double DCy = d * ((minY + maxY) / 2.0);
        double dx = (imageX / 2.0) - DCx;
        double dy = (imageY / 2.0) - DCy;
        img::EasyImage image((unsigned int) std::round(imageX),
                             (unsigned int) std::round(imageY),
                             img::Color(bc.red_int_value(), bc.green_int_value(), bc.blue_int_value()));
        col::ZBuffer* z_buffer;
        switch (type) {
            case Wires :
                for(prj::Line &l : lines){
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
                    draw_lines(image, x0, y0, x1, y1, l.color);
                }
                break;
            case ZBuff :
                z_buffer = new col::ZBuffer((unsigned int) std::round(imageX), (unsigned int) std::round(imageY));
                for (prj::Line &l : lines) {
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
                    draw_zbuf_lines(image, (*z_buffer), x0, y0, l.p1.z, x1, y1, l.p2.z, l.color);
                }
                delete z_buffer;
                break;
            case Trian :
                z_buffer = new col::ZBuffer((unsigned int) std::round(imageX), (unsigned int) std::round(imageY));
                // Draw images
                for (auto& figure : figures) {
                    figure->triangulate();
                    for (auto& face : figure->faces) {
                        Vector3D A = figure->get_point_at(face[0]);
                        Vector3D B = figure->get_point_at(face[1]);
                        Vector3D C = figure->get_point_at(face[2]);
                        draw_zbuf_triangles(image, (*z_buffer),
                                            A, B, C,
                                            d, dx, dy,
                                            figure->get_ambient(),
                                            figure->get_diffuse(),
                                            figure->get_specular(),
                                            figure->get_reflection(), lights);
                    }
                }
                delete z_buffer;
                break;
        }
    return image;
    }

    img::EasyImage draw(prj::Lines& lines, const unsigned int& size, col::Color& bc) {
        double maxX = lines[0].p1.x;
        double maxY = lines[0].p1.y;
        double minX = lines[0].p1.x;
        double minY = lines[0].p1.y;
        for (prj::Line &l : lines) {
            if (l.p1.x > maxX) maxX = l.p1.x;
            if (l.p1.x < minX) minX = l.p1.x;
            if (l.p2.x > maxX) maxX = l.p2.x;
            if (l.p2.x < minX) minX = l.p2.x;
            if (l.p1.y > maxY) maxY = l.p1.y;
            if (l.p1.y < minY) minY = l.p1.y;
            if (l.p2.y > maxY) maxY = l.p2.y;
            if (l.p2.y < minY) minY = l.p2.y;
        }
        double xRange = maxX - minX;
        double yRange = maxY - minY;
        double maxXY = std::max(xRange, yRange);
        double imageX = size * (xRange / maxXY);
        double imageY = size * (yRange / maxXY);
        double d = 0.95 * (imageX / xRange);
        double DCx = d * ((minX + maxX) / 2.0);
        double DCy = d * ((minY + maxY) / 2.0);
        double dx = (imageX / 2.0) - DCx;
        double dy = (imageY / 2.0) - DCy;
        img::EasyImage image((unsigned int) imageX,
                             (unsigned int) imageY,
                             img::Color(bc.red_int_value(), bc.green_int_value(), bc.blue_int_value()));
        for(prj::Line &l : lines){
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
            draw_lines(image, x0, y0, x1, y1, l.color);
        }
        return image;
    }

    void draw_lines(
            img::EasyImage& image,
            unsigned int& x0, unsigned int& y0,
            unsigned int& x1, unsigned int& y1,
            col::Color color) {
        assert(x0 < image.get_width() && y0 < image.get_height());
        assert(x1 < image.get_width() && y1 < image.get_height());
        img::Color img_color(color.red_int_value(), color.green_int_value(), color.blue_int_value());
        if (x0 == x1) {
            //special case for x0 == x1
            for (unsigned int i = std::min(y0, y1); i <= std::max(y0, y1); i++) {
                image(x0, i) = img_color;
            }
        }
        else if (y0 == y1) {
            //special case for y0 == y1
            for (unsigned int i = std::min(x0, x1); i <= std::max(x0, x1); i++) {
                image(i, y0) = img_color;
            }
        }
        else {
            if (x0 > x1) {
                //flip points if x1>x0: we want x0 to have the lowest value
                std::swap(x0, x1);
                std::swap(y0, y1);
            }
            double m = ((double) y1 - (double) y0) / ((double) x1 - (double) x0);
            if (-1.0 <= m && m <= 1.0) {
                for (unsigned int i = 0; i <= (x1 - x0); i++) {
                    image(x0 + i, (unsigned int) round(y0 + m * i)) = img_color;
                }
            }
            else if (m > 1.0) {
                for (unsigned int i = 0; i <= (y1 - y0); i++) {
                    image((unsigned int) round(x0 + (i / m)), y0 + i) = img_color;
                }
            }
            else if (m < -1.0) {
                for (unsigned int i = 0; i <= (y0 - y1); i++) {
                    image((unsigned int) round(x0 - (i / m)), y0 - i) = img_color;
                }
            }
        }
    }

    void draw_zbuf_lines(
            img::EasyImage& image, col::ZBuffer& z_buffer,
            unsigned int& x0, unsigned int& y0, double& z0,
            unsigned int& x1, unsigned int& y1, double& z1,
            col::Color color) {
        assert(x0 < image.get_width() && y0 < image.get_height());
        assert(x1 < image.get_width() && y1 < image.get_height());
        img::Color img_color(color.red_int_value(), color.green_int_value(), color.blue_int_value());
        if (x0 == x1)
        {
            int a = std::max(y0, y1) - std::min(y0, y1);
            int k = a;
            //special case for x0 == x1
            for (unsigned int i = std::min(y0, y1); i <= std::max(y0, y1); i++)
            {
                double z_inv = (((double) k / (double) a) / z0) + ((1 - ((double) k / (double) a)) / z1);
                if (z_inv < z_buffer(x0, i)) {
                    z_buffer(x0, i) = z_inv;
                    image(x0, i) = img_color;
                }
                k--;
            }
        }
        else if (y0 == y1)
        {
            int a = std::max(x0, x1) - std::min(x0, x1);
            int k = a;
            //special case for y0 == y1
            for (unsigned int i = std::min(x0, x1); i <= std::max(x0, x1); i++)
            {
                double z_inv = (((double) k / (double) a) / z0) + ((1 - ((double) k / (double) a)) / z1);
                if (z_inv < z_buffer(i, y0)) {
                    z_buffer(i, y0) = z_inv;
                    image(i, y0) = img_color;
                }
                k--;
            }
        }
        else
        {
            if (x0 > x1)
            {
                //flip points if x1>x0: we want x0 to have the lowest value
                std::swap(x0, x1);
                std::swap(y0, y1);
                std::swap(z0, z1);
            }
            double m = ((double) y1 - (double) y0) / ((double) x1 - (double) x0);
            if (-1.0 <= m && m <= 1.0)
            {
                int a = x1 - x0;
                int k = a;
                for (unsigned int i = 0; i <= (x1 - x0); i++)
                {
                    double z_inv = (((double) k / (double) a) / z0) + ((1 - ((double) k / (double) a)) / z1);
                    if (z_inv < z_buffer(x0 + i, (unsigned int) round(y0 + m * i))) {
                        z_buffer(x0 + i, (unsigned int) round(y0 + m * i)) = z_inv;
                        image(x0 + i, (unsigned int) round(y0 + m * i)) = img_color;
                    }
                    k--;
                }
            }
            else if (m > 1.0)
            {
                int a = y1 - y0;
                int k = a;
                for (unsigned int i = 0; i <= (y1 - y0); i++)
                {
                    double z_inv = (((double) k / (double) a) / z0) + ((1 - ((double) k / (double) a)) / z1);
                    if (z_inv < z_buffer((unsigned int) round(x0 + (i / m)), y0 + i)) {
                        z_buffer((unsigned int) round(x0 + (i / m)), y0 + i) = z_inv;
                        image((unsigned int) round(x0 + (i / m)), y0 + i) = img_color;
                    }
                    k--;
                }
            }
            else if (m < -1.0)
            {
                int a = y0 - y1;
                int k = a;
                for (unsigned int i = 0; i <= (y0 - y1); i++)
                {
                    double z_inv = (((double) k / (double) a) / z0) + ((1 - ((double) k / (double) a)) / z1);
                    if (z_inv < z_buffer((unsigned int) round(x0 - (i / m)), y0 - i)) {
                        z_buffer((unsigned int) round(x0 - (i / m)), y0 - i) = z_inv;
                        image((unsigned int) round(x0 - (i / m)), y0 - i) = img_color;
                    }
                    k--;
                }
            }
        }
    }

    void draw_zbuf_triangles(
            img::EasyImage& image, col::ZBuffer& z_buffer,
            Vector3D& A, Vector3D& B, Vector3D& C,
            const double& d, const double& dx, const double& dy,
            col::Color ambientReflection,
            col::Color diffuseReflection,
            col::Color specularReflection,
            double reflectionCoeff, col::Lights& lights) {

        double xA = -((d * A.x) / (A.z)) + dx;
        double yA = -((d * A.y) / (A.z)) + dy;
        double xB = -((d * B.x) / (B.z)) + dx;
        double yB = -((d * B.y) / (B.z)) + dy;
        double xC = -((d * C.x) / (C.z)) + dx;
        double yC = -((d * C.y) / (C.z)) + dy;

        double xG = (xA + xB + xC) / 3.0;
        double yG = (yA + yB + yC) / 3.0;
        double zG_inverted = (1.0 / (3.0 * A.z)) + (1.0 / (3.0 * B.z)) + (1.0 / (3.0 * C.z));

        int yMin = (int) std::round(std::min(std::min(yA, yB), yC) + 0.5);
        int yMax = (int) std::round(std::max(std::max(yA, yB), yC) - 0.5);

        Vector3D u = B - A;
        Vector3D v = C - A;
        Vector3D w; w = u.cross_equals(v);
        double k = (w.x * A.x) + (w.y * A.y) + (w.z * A.z);
        double dzdx = -(w.x / (d * k));
        double dzdy = -(w.y / (d * k));

        // Color
        col::Color color_values = {0, 0, 0};
        Vector3D n = w;
        n.normalise();
        std::vector<col::Light*> pnt_spc_lights;
        for (auto& light : lights) {
            if (light->is_ambient()) {
                color_values += (ambientReflection * light->ambient());
            }
            if (light->is_diffuse_inf()) {
                double cos_alpha = 0;
                Vector3D l = light->get_direction();
                l.normalise();
                l = -l;
                double val = n.dot(l);
                if (val > 0) cos_alpha = val;
                color_values += (diffuseReflection * light->diffuse() * cos_alpha);
            }
            if (light->is_diffuse_pnt()) {
                pnt_spc_lights.push_back(light);
            }
            if (light->is_specular() and
                std::find(pnt_spc_lights.begin(), pnt_spc_lights.end(), light) == pnt_spc_lights.end()) {
                pnt_spc_lights.push_back(light);
            }
        }
        img::Color color(color_values.red_int_value(), color_values.green_int_value(), color_values.blue_int_value());

        for (unsigned int yI = (unsigned int) yMin; yI <= yMax; yI++) {
            double xL_AB = inf;
            double xL_AC = inf;
            double xL_BC = inf;
            double xR_AB = -inf;
            double xR_AC = -inf;
            double xR_BC = -inf;

            // PQ == AB
            if ((((yI - yA) * (yI - yB)) <= 0) and yA != yB) {
                // double xI = xB + ((xA - xB) * ((yI - yB) / (yA - yB)));
                double xI = xA + ((xB - xA) * ((yI - yA) / (yB - yA)));
                xL_AB = xI;
                xR_AB = xI;
            }

            // PQ == AC
            if ((((yI - yA) * (yI - yC)) <= 0) and yA != yC) {
                // double xI = xC + ((xA - xC) * ((yI - yC) / (yA - yC)));
                double xI = xA + ((xC - xA) * ((yI - yA) / (yC - yA)));
                xL_AC = xI;
                xR_AC = xI;
            }

            // PQ == BC
            if ((((yI - yB) * (yI - yC)) <= 0) and yB != yC) {
                // double xI = xC + ((xB - xC) * ((yI - yC) / (yB - yC)));
                double xI = xB + ((xC - xB) * ((yI - yB) / (yC - yB)));
                xL_BC = xI;
                xR_BC = xI;
            }
            int xL = (int) std::round(std::min(std::min(xL_AB, xL_AC), xL_BC) + 0.5);
            int xR = (int) std::round(std::max(std::max(xR_AB, xR_AC), xR_BC) - 0.5);

            // if (xL < 0) continue;
            double z_inv_part = 1.0001 * zG_inverted;

            for (unsigned int xI = (unsigned int) xL; xI <= xR; xI++) {
                double z_inv = z_inv_part + (((double) xI - xG) * dzdx) + (((double) yI - yG) * dzdy);
                if (z_inv < z_buffer(xI, yI)) {
                    z_buffer(xI, yI) = z_inv;
                    double cos_alpha = 0;
                    Vector3D l = Vector3D::vector(0, 0, 0);
                    double zEye = 1.0 / z_inv;
                    double xEye = (xI - dx) * (-zEye) / d;
                    double yEye = (yI - dy) * (-zEye) / d;
                    Vector3D P = Vector3D::point(xEye, yEye, zEye);
                    col::Color pnt_spec_color = color_values;
                      for (auto light : pnt_spc_lights) {
                        if (light->is_diffuse_inf()) {
                            l = Vector3D::normalise(-light->get_direction());
                        }
                        else {
                            l = Vector3D::normalise(light->get_location() - P);
                        }
                        cos_alpha = n.dot(l);
                        cos_alpha =  cos_alpha < 0 ? 0 : cos_alpha;
                        pnt_spec_color += (diffuseReflection * light->diffuse()) * cos_alpha;
                        if (light->is_specular() && cos_alpha > 0) {
                            Vector3D r = Vector3D::normalise((2 * cos_alpha * n) + l);
                            Vector3D cam = Vector3D::normalise(-P);
                            double cos_beta = std::pow(r.dot(cam), reflectionCoeff);
                            cos_beta = cos_beta < 0 ? 0 : cos_beta;
                            pnt_spec_color += (specularReflection * light->specular()) * cos_beta;
                        }
                    }
                    color = img::Color(pnt_spec_color.red_int_value(),
                                       pnt_spec_color.green_int_value(),
                                       pnt_spec_color.blue_int_value());
                    image(xI, yI) = color;
                }
            }
        }

    }
}
