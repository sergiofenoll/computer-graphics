//
// Created by sergio on 30/04/17.
//

#include "draw.hh"

namespace drw {

    img::EasyImage draw(
            fig::Figures& figures, DrawType& type, col::Lights& lights,
            Matrix& eyeMatrix, const unsigned int& size, col::Color& bc) {
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
                             img::Color(bc.red_int_value(), 
                                        bc.green_int_value(), 
                                        bc.blue_int_value()));
        col::ZBuffer* z_buffer;
        for (auto light : lights[col::Shadow]) {
            Vector3D loc = light->get_location() * eyeMatrix;
            light->set_location(loc);
        }
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
                z_buffer = new col::ZBuffer((unsigned int) imageX, 
                                            (unsigned int) imageY);
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
                    draw_zbuf_lines(image, (*z_buffer),
                                    x0, y0, l.p1.z,
                                    x1, y1, l.p2.z,
                                    l.color);
                }
                delete z_buffer;
                break;
            case Triangles :
                z_buffer = new col::ZBuffer((unsigned int) imageX, 
                                            (unsigned int) imageY);
                for (auto& figure : figures) {
                    figure->triangulate();
                    for (auto& face : figure->faces) {
                        Vector3D A = figure->get_point_at(face[0]);
                        Vector3D B = figure->get_point_at(face[1]);
                        Vector3D C = figure->get_point_at(face[2]);
                        for (auto light : lights[col::Shadow]) {
                            Vector3D A_light = A * Matrix::inv(eyeMatrix) * light->get_eye();
                            Vector3D B_light = B * Matrix::inv(eyeMatrix) * light->get_eye();
                            Vector3D C_light = C * Matrix::inv(eyeMatrix) * light->get_eye();
                            setup_shadow_mask(A_light, B_light, C_light, light);
                        }
                        draw_zbuf_triangles(image, (*z_buffer),
                                            eyeMatrix, A, B, C,
                                            d, dx, dy, figure, lights);
                    }
                }
                delete z_buffer;
                break;
            case Opaque :
                z_buffer = new col::ZBuffer((unsigned int) imageX,
                                            (unsigned int) imageY);
                col::OpacityMatrix opacity_matrix((unsigned int) imageX,
                                                  (unsigned int) imageY);
                for (auto& figure : figures) {
                    figure->triangulate();
                    for (auto& face : figure->faces) {
                        Vector3D A = figure->get_point_at(face[0]);
                        Vector3D B = figure->get_point_at(face[1]);
                        Vector3D C = figure->get_point_at(face[2]);
                        draw_transparent(image, (*z_buffer), opacity_matrix,
                                         eyeMatrix, A, B, C,
                                         d, dx, dy, figure, lights);
                    }
                }
                delete z_buffer;
                break;
        }
    return image;
    }

    img::EasyImage draw(prj::Lines& lines,
            const unsigned int& size,
            col::Color& bc) {
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
                             img::Color(bc.red_int_value(),
                                        bc.green_int_value(),
                                        bc.blue_int_value()));
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
        img::Color img_color(color.red_int_value(),
                             color.green_int_value(),
                             color.blue_int_value());
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
        img::Color img_color(color.red_int_value(),
                             color.green_int_value(),
                             color.blue_int_value());
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
            Matrix& eyeMatrix, Vector3D& A, Vector3D& B, Vector3D& C,
            const double& d, const double& dx, const double& dy,
            fig::Figure* fig, col::Lights& lights) {
        double dzdx, dzdy;
        Vector3D n;
        {
            Vector3D u = B - A;
            Vector3D v = C - A;
            Vector3D w = Vector3D::cross(u, v);
            if (w.z < 0) return; // https://www.tutorialspoint.com/computer_graphics/visible_surface_detection.htm
            double k = Vector3D::dot(w, A);
            dzdx = -(w.x / (d *  k));
            dzdy = -(w.y / (d * k));
            n = Vector3D::normalise(w);
        }

        double xA = -((d * A.x) / (A.z)) + dx;
        double yA = -((d * A.y) / (A.z)) + dy;
        double xB = -((d * B.x) / (B.z)) + dx;
        double yB = -((d * B.y) / (B.z)) + dy;
        double xC = -((d * C.x) / (C.z)) + dx;
        double yC = -((d * C.y) / (C.z)) + dy;

        double xG = (xA + xB + xC) / 3.0;
        double yG = (yA + yB + yC) / 3.0;
        double zG_inverted =
                (1.0 / (3.0 * A.z)) + (1.0 / (3.0 * B.z)) + (1.0 / (3.0 * C.z));

        int yMin = (int) std::round(std::min(std::min(yA, yB), yC) + 0.5);
        int yMax = (int) std::round(std::max(std::max(yA, yB), yC) - 0.5);

        // Color
        col::Color color_values = {0, 0, 0};
        std::set<col::Light*> pnt_spec_lights;
        apply_face_light(color_values, pnt_spec_lights, lights, fig, n);

        for (unsigned int yI = (unsigned int) yMin; yI <= yMax; yI++) {
            int xL, xR;
            {
                double xL_AB = inf;
                double xL_AC = inf;
                double xL_BC = inf;
                double xR_AB = -inf;
                double xR_AC = -inf;
                double xR_BC = -inf;
                // PQ == AB
                if ((((yI - yA) * (yI - yB)) <= 0) and yA != yB) {
                    double xI = xA + ((xB - xA) * ((yI - yA) / (yB - yA)));
                    xL_AB = xI;
                    xR_AB = xI;
                }
                // PQ == AC
                if ((((yI - yA) * (yI - yC)) <= 0) and yA != yC) {
                    double xI = xA + ((xC - xA) * ((yI - yA) / (yC - yA)));
                    xL_AC = xI;
                    xR_AC = xI;
                }
                // PQ == BC
                if ((((yI - yB) * (yI - yC)) <= 0) and yB != yC) {
                    double xI = xB + ((xC - xB) * ((yI - yB) / (yC - yB)));
                    xL_BC = xI;
                    xR_BC = xI;
                }
                xL = (int) std::round(
                        std::min(std::min(xL_AB, xL_AC), xL_BC) + 0.5);
                xR = (int) std::round(
                        std::max(std::max(xR_AB, xR_AC), xR_BC) - 0.5);
            }

            double z_inv_constant = 1.0001 * zG_inverted;
            double z_inv_y_term = ((double) yI - yG) * dzdy;

            for (unsigned int xI = (unsigned int) xL; xI <= xR; xI++) {
                double z_inv =
                        z_inv_constant + (((double) xI - xG) * dzdx) +
                        z_inv_y_term;
                double z_inv_accurate =
                        zG_inverted + (((double) xI - xG) * dzdx) +
                        z_inv_y_term;
                if (z_inv < z_buffer(xI, yI)) {
                    z_buffer(xI, yI) = z_inv;
                    col::Color final_color = color_values;
                    if (pnt_spec_lights.size() > 0) {
                        apply_point_light(color_values, final_color,
                                          pnt_spec_lights, fig, eyeMatrix,
                                          z_inv, z_inv_accurate,
                                          xI, yI, dx, dy, d, n);
                    }
                    img::Color color(final_color.red_int_value(),
                                     final_color.green_int_value(),
                                     final_color.blue_int_value());
                    image(xI, yI) = color;
                }
            }
        }
    }

    void draw_transparent(img::EasyImage& image, col::ZBuffer& z_buffer,
                          col::OpacityMatrix& opacity_matrix,
                          Matrix& eyeMatrix, Vector3D& A, Vector3D& B, Vector3D& C,
                          const double& d, const double& dx, const double& dy,
                          fig::Figure* fig, col::Lights& lights) {
        double xA = -((d * A.x) / (A.z)) + dx;
        double yA = -((d * A.y) / (A.z)) + dy;
        double xB = -((d * B.x) / (B.z)) + dx;
        double yB = -((d * B.y) / (B.z)) + dy;
        double xC = -((d * C.x) / (C.z)) + dx;
        double yC = -((d * C.y) / (C.z)) + dy;

        double xG = (xA + xB + xC) / 3.0;
        double yG = (yA + yB + yC) / 3.0;
        double zG_inverted =
                (1.0 / (3.0 * A.z)) + (1.0 / (3.0 * B.z)) + (1.0 / (3.0 * C.z));

        int yMin = (int) std::round(std::min(std::min(yA, yB), yC) + 0.5);
        int yMax = (int) std::round(std::max(std::max(yA, yB), yC) - 0.5);

        double dzdx, dzdy;
        Vector3D n;
        {
            Vector3D u = B - A;
            Vector3D v = C - A;
            Vector3D w = Vector3D::cross(u, v);
            double k = Vector3D::dot(w, A);
            dzdx = -(w.x / (d * k));
            dzdy = -(w.y / (d * k));
            n = Vector3D::normalise(w);
        }

        // Color
        col::Color color_values = {0, 0, 0};
        std::set<col::Light*> pnt_spec_lights;
        apply_face_light(color_values, pnt_spec_lights, lights, fig, n);

        for (unsigned int yI = (unsigned int) yMin; yI <= yMax; yI++) {
            int xL, xR;
            {
                double xL_AB = inf;
                double xL_AC = inf;
                double xL_BC = inf;
                double xR_AB = -inf;
                double xR_AC = -inf;
                double xR_BC = -inf;
                // PQ == AB
                if ((((yI - yA) * (yI - yB)) <= 0) and yA != yB) {
                    double xI = xA + ((xB - xA) * ((yI - yA) / (yB - yA)));
                    xL_AB = xI;
                    xR_AB = xI;
                }
                // PQ == AC
                if ((((yI - yA) * (yI - yC)) <= 0) and yA != yC) {
                    double xI = xA + ((xC - xA) * ((yI - yA) / (yC - yA)));
                    xL_AC = xI;
                    xR_AC = xI;
                }
                // PQ == BC
                if ((((yI - yB) * (yI - yC)) <= 0) and yB != yC) {
                    double xI = xB + ((xC - xB) * ((yI - yB) / (yC - yB)));
                    xL_BC = xI;
                    xR_BC = xI;
                }
                xL = (int) std::round(
                        std::min(std::min(xL_AB, xL_AC), xL_BC) + 0.5);
                xR = (int) std::round(
                        std::max(std::max(xR_AB, xR_AC), xR_BC) - 0.5);
            }

            double z_inv_constant = 1.0001 * zG_inverted;
            double z_inv_y_term = ((double) yI - yG) * dzdy;

            for (unsigned int xI = (unsigned int) xL; xI <= xR; xI++) {
                double z_inv =
                        z_inv_constant + (((double) xI - xG) * dzdx) +
                        z_inv_y_term;
                if (z_inv < z_buffer(xI, yI)) {
                    z_buffer(xI, yI) = z_inv;
                    col::Color final_color = color_values;
                    if (pnt_spec_lights.size() > 0) {
                        apply_point_light(color_values, final_color,
                                          pnt_spec_lights, fig, eyeMatrix,
                                          z_inv, 0, xI, yI, dx, dy, d, n);
                    }
                    if (fig->opacity == 1) {
                        img::Color color(final_color.red_int_value(),
                                         final_color.green_int_value(),
                                         final_color.blue_int_value());
                        image(xI, yI) = color;
                        opacity_matrix(xI, yI) = fig->opacity;
                    } else {
                        col::Color old_color = {image(xI, yI).red / 255.0,
                                                image(xI, yI).green / 255.0,
                                                image(xI, yI).blue / 255.0};
                        col::Color new_color = final_color * fig->opacity +
                                                 old_color * (1 - fig->opacity);
                        img::Color color(new_color.red_int_value(),
                                         new_color.green_int_value(),
                                         new_color.blue_int_value());
                        image(xI, yI) = color;
                        opacity_matrix(xI, yI) = fig->opacity;
                    }
                } else if (opacity_matrix(xI, yI) < 1) {
                    col::Color final_color = color_values;
                    if (pnt_spec_lights.size() > 0) {
                        apply_point_light(color_values, final_color,
                                          pnt_spec_lights, fig, eyeMatrix,
                                          z_inv, 0, xI, yI, dx, dy, d, n);
                    }
                    col::Color new_color = {image(xI, yI).red / 255.0,
                                            image(xI, yI).green / 255.0,
                                            image(xI, yI).blue / 255.0};
                    col::Color combined_color =
                            new_color * opacity_matrix(xI, yI) +
                            final_color * (1 - opacity_matrix(xI, yI));
                    img::Color color(combined_color.red_int_value(),
                                     combined_color.green_int_value(),
                                     combined_color.blue_int_value());
                    image(xI, yI) = color;
                }
            }
        }
    }

    void apply_face_light(col::Color& color_values,
                          std::set<col::Light*>& pnt_spec_lights,
                          col::Lights& lights,
                          fig::Figure* fig,
                          Vector3D& n) {
        for (auto& light : lights[col::Normal]) {
            if (light->is_ambient()) {
                color_values += (fig->ambient_reflection * light->ambient());
            }
            if (light->is_diffuse_inf()) {
                Vector3D l = Vector3D::normalise(-light->get_direction());
                double cos_alpha = Vector3D::dot(n, l);
                if (cos_alpha > 0)
                    color_values += (fig->diffuse_reflection *
                                     light->diffuse() * cos_alpha);
            }
            if (light->is_diffuse_pnt())
                pnt_spec_lights.insert(light);
            if (light->is_specular())
                pnt_spec_lights.insert(light);
        }
        for (auto& light : lights[col::Shadow   ]) {
            if (light->is_ambient()) {
                color_values += (fig->ambient_reflection * light->ambient());
            }
            if (light->is_diffuse_inf()) {
                Vector3D l = Vector3D::normalise(-light->get_direction());
                double cos_alpha = Vector3D::dot(n, l);
                if (cos_alpha > 0)
                    color_values += (fig->diffuse_reflection *
                                     light->diffuse() * cos_alpha);
            }
            if (light->is_diffuse_pnt())
                pnt_spec_lights.insert(light);
            if (light->is_specular())
                pnt_spec_lights.insert(light);
        }
    }

    void apply_point_light(col::Color& color_values,
                           col::Color& final_color,
                           std::set<col::Light*>& pnt_spec_lights,
                           fig::Figure* fig, Matrix& eyeMatrix,
                           double z_inv, double z_inv_accurate,
                           unsigned int xI, unsigned int yI,
                           double dx, double dy, double d,
                           Vector3D& n) {
        Vector3D l = Vector3D::vector(0, 0, 0);
        double zEye = 1.0 / z_inv;
        double xEye = (xI - dx) * (-zEye) / d;
        double yEye = (yI - dy) * (-zEye) / d;
        Vector3D pEye = Vector3D::point(xEye, yEye, zEye);
        Vector3D cam = Vector3D::normalise(
                Vector3D::point(0, 0, 0) - pEye);
        for (auto light : pnt_spec_lights) {
            if (light->is_diffuse_inf()) {
                l = Vector3D::normalise(
                        -light->get_direction());
            } else {
                l = Vector3D::normalise(
                        light->get_location() - pEye);
            }
            double cos_alpha = Vector3D::dot(n, l);
            if (cos_alpha > 0) {
                if (light->is_shadow()) {
                    double dL, dXL, dYL;
                    light->get_values(dL, dXL, dYL);
                    zEye = 1.0 / z_inv_accurate;
                    xEye = (xI - dx) * (-zEye) / d;
                    yEye = (yI - dy) * (-zEye) / d;
                    pEye = Vector3D::point(xEye, yEye, zEye);
                    Vector3D pShadow =
                            pEye * Matrix::inv(eyeMatrix) * light->get_eye();
                    double xShadow = (pShadow.x * dL / -pShadow.z) + dXL;
                    double yShadow = (pShadow.y * dL / -pShadow.z) + dYL;
                    double z_inv_Shadow;
                    {
                        double alphaX = xShadow - std::floor(xShadow);
                        double alphaY = yShadow - std::floor(yShadow);
                        double z_inv_A = light->get_shadow_mask_at(
                                std::floor(xShadow), std::floor(yShadow));
                        double z_inv_B = light->get_shadow_mask_at(
                                std::ceil(xShadow), std::floor(yShadow));
                        double z_inv_C = light->get_shadow_mask_at(
                                std::floor(xShadow), std::ceil(yShadow));
                        double z_inv_D = light->get_shadow_mask_at(
                                std::ceil(xShadow), std::ceil(yShadow));
                        double z_inv_E =
                                ((1.0 - alphaX) * z_inv_A) + (alphaX * z_inv_B);
                        double z_inv_F =
                                ((1.0 - alphaX) * z_inv_C) + (alphaX * z_inv_D);
                        z_inv_Shadow = ((1.0 - alphaY) * z_inv_E) + (alphaY * z_inv_F);
                    }
                    if (std::fabs(z_inv_Shadow - 1.0/pShadow.z) <= 1e-5) {
                        if (light->is_diffuse_pnt())
                            final_color +=
                                    (fig->diffuse_reflection *
                                     light->diffuse()) *
                                    cos_alpha;
                        if (light->is_specular()) {
                            Vector3D r = Vector3D::normalise(
                                    (2 * cos_alpha * n) - l);
                            double cos_beta = Vector3D::dot(r, cam);
                            if (cos_beta > 0) {
                                double pow = std::pow(cos_beta,
                                                      fig->reflection_coefficient);
                                final_color +=
                                        (fig->specular_reflection *
                                         light->specular()) * pow;
                            }
                        }
                    }
                }
                else {
                    if (light->is_diffuse_pnt())
                        final_color +=
                                (fig->diffuse_reflection *
                                 light->diffuse()) *
                                cos_alpha;
                    if (light->is_specular()) {
                        Vector3D r = Vector3D::normalise(
                                (2 * cos_alpha * n) - l);
                        double cos_beta = Vector3D::dot(r, cam);
                        if (cos_beta > 0) {
                            double pow = std::pow(cos_beta,
                                                  fig->reflection_coefficient);
                            final_color +=
                                    (fig->specular_reflection *
                                     light->specular()) * pow;
                        }
                    }
                }
            }
        }
    }

    void setup_shadow_mask(Vector3D& A, Vector3D& B, Vector3D& C, col::Light* light) {
        double d, dx, dy;
        light->get_values(d, dx, dy);
        double xA = -((d * A.x) / (A.z)) + dx;
        double yA = -((d * A.y) / (A.z)) + dy;
        double xB = -((d * B.x) / (B.z)) + dx;
        double yB = -((d * B.y) / (B.z)) + dy;
        double xC = -((d * C.x) / (C.z)) + dx;
        double yC = -((d * C.y) / (C.z)) + dy;

        double xG = (xA + xB + xC) / 3.0;
        double yG = (yA + yB + yC) / 3.0;
        double zG_inverted =
                (1.0 / (3.0 * A.z)) + (1.0 / (3.0 * B.z)) + (1.0 / (3.0 * C.z));

        int yMin = (int) std::round(std::min(std::min(yA, yB), yC) + 0.5);
        int yMax = (int) std::round(std::max(std::max(yA, yB), yC) - 0.5);

        double dzdx, dzdy;
        {
            Vector3D u = B - A;
            Vector3D v = C - A;
            Vector3D w = Vector3D::cross(u, v);
            double k = Vector3D::dot(w, A);
            dzdx = -(w.x / (d * k));
            dzdy = -(w.y / (d * k));
        }

        for (unsigned int yI = (unsigned int) yMin; yI <= yMax; yI++) {
            int xL, xR;
            {
                double xL_AB = inf;
                double xL_AC = inf;
                double xL_BC = inf;
                double xR_AB = -inf;
                double xR_AC = -inf;
                double xR_BC = -inf;
                // PQ == AB
                if ((((yI - yA) * (yI - yB)) <= 0) and yA != yB) {
                    double xI = xA + ((xB - xA) * ((yI - yA) / (yB - yA)));
                    xL_AB = xI;
                    xR_AB = xI;
                }
                // PQ == AC
                if ((((yI - yA) * (yI - yC)) <= 0) and yA != yC) {
                    double xI = xA + ((xC - xA) * ((yI - yA) / (yC - yA)));
                    xL_AC = xI;
                    xR_AC = xI;
                }
                // PQ == BC
                if ((((yI - yB) * (yI - yC)) <= 0) and yB != yC) {
                    double xI = xB + ((xC - xB) * ((yI - yB) / (yC - yB)));
                    xL_BC = xI;
                    xR_BC = xI;
                }
                xL = (int) std::round(
                        std::min(std::min(xL_AB, xL_AC), xL_BC) + 0.5);
                xR = (int) std::round(
                        std::max(std::max(xR_AB, xR_AC), xR_BC) - 0.5);
            }
            double z_inv_y_term = ((double) yI - yG) * dzdy;
            for (unsigned int xI = (unsigned int) xL; xI <= xR; xI++) {
                double z_inv = zG_inverted +
                               (((double) xI - xG) * dzdx) + z_inv_y_term;
                if (z_inv < light->get_shadow_mask_at(xI, yI))
                    light->set_shadow_mask_at(xI, yI, z_inv);
            }
        }
    }
}
