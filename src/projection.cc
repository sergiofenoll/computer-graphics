//
// Created by sergio on 30/04/17.
//

#include "projection.hh"

namespace prj {
    Point::Point() {}

    Point::Point(const double& x, const double& y) :
    x(x), y(y){}

    Point::Point(const Vector3D& point, const double& d) {
        (*this).x = (d * point.x) / (-point.z);
        (*this).y = (d * point.y) / (-point.z);
        (*this).z = point.z;
    }

    void Point::project(const Vector3D& point, const double& d) {
        (*this).x = (d * point.x) / -(point.z);
        (*this).y = (d * point.y) / -(point.z);
        (*this).z = point.z;
    }

    Line::Line() {}

    Line::Line(Point p1, Point p2, col::Color color) {
        (*this).p1=p1;
        (*this).p2=p2;
        (*this).color=color;
    }

    Lines project_figures(fig::Figures& figures) {
        Lines lines;
        for (fig::Figure* figure : figures) {
            for (fig::Face& face : figure->faces) {
                Points points;
                for (unsigned int index : face.point_indexes){
                    Point p(figure->points[index]);
                    points.push_back(p);
                }
                for (int i=0; i<points.size(); i++) {
                    Line line;
                    line.p1 = points[i];
                    line.p2 = points[(i + 1) % points.size()];
                    line.color = figure->get_ambient();
                    lines.push_back(line);
                }
            }
        }
        return lines;
    }

    void recursivePrintString(
            LParser::LSystem2D& l_system, std::string& print_string,
            unsigned int& cur_it, unsigned int& max_it, double& alph_angle, double& delt_angle,
            Lines& lines, std::stack<double>& brStack, col::Color& lc, double& x, double& y){
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
                                         lines, brStack, lc, x, y);
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
                    Point p1 (x, y);
                    x += std::cos(alph_angle);
                    y += std::sin(alph_angle);
                    Point p2 (x, y);
                    Line line(p1, p2, lc);
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
            LParser::LSystem3D& l_system, std::string& print_string,
            unsigned int& cur_it, unsigned int& max_it, double& delt_angle,
            fig::Figure& figure, std::stack<double>& doubleStack, std::stack<Vector3D>& vectorStack,
            double& x, double& y, double& z, Vector3D& H, Vector3D& L, Vector3D& U){
        for (char &c : print_string){
            // Iterate over characters of string
            if (cur_it < max_it) {
                // Recursion ahead!
                if (c=='+') {
                    Vector3D H_New; H_New = H * std::cos(delt_angle) + L * std::sin(delt_angle);
                    Vector3D L_New; L_New = -H * std::sin(delt_angle) + L * std::cos(delt_angle);
                    H = H_New;
                    L = L_New;
                } else if (c=='-') {
                    Vector3D H_New; H_New = H * std::cos(-delt_angle) + L * std::sin(-delt_angle);
                    Vector3D L_New; L_New = -H * std::sin(-delt_angle) + L * std::cos(-delt_angle);
                    H = H_New;
                    L = L_New;
                } else if (c == '^') {
                    Vector3D H_New; H_New = H * std::cos(delt_angle) + U * std::sin(delt_angle);
                    Vector3D U_New; U_New = -H * std::sin(delt_angle) + U * std::cos(delt_angle);
                    H = H_New;
                    U = U_New;
                } else if (c == '&') {
                    Vector3D H_New; H_New = H * std::cos(-delt_angle) + U * std::sin(-delt_angle);
                    Vector3D U_New; U_New = -H * std::sin(-delt_angle) + U * std::cos(-delt_angle);
                    H = H_New;
                    U = U_New;
                } else if (c == '\\') {
                    Vector3D L_New; L_New = L * std::cos(delt_angle) - U * std::sin(delt_angle);
                    Vector3D U_New; U_New = L * std::sin(delt_angle) + U * std::cos(delt_angle);
                    L = L_New;
                    U = U_New;
                } else if (c == '/') {
                    Vector3D L_New; L_New = L * std::cos(-delt_angle) - U * std::sin(-delt_angle);
                    Vector3D U_New; U_New = L * std::sin(-delt_angle) + U * std::cos(-delt_angle);
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
                    Vector3D H_New; H_New = H * std::cos(delt_angle) + L * std::sin(delt_angle);
                    Vector3D L_New; L_New = -H * std::sin(delt_angle) + L * std::cos(delt_angle);
                    H = H_New;
                    L = L_New;
                } else if (c=='-') {
                    Vector3D H_New; H_New = H * std::cos(-delt_angle) + L * std::sin(-delt_angle);
                    Vector3D L_New; L_New = -H * std::sin(-delt_angle) + L * std::cos(-delt_angle);
                    H = H_New;
                    L = L_New;
                } else if (c == '^') {
                    Vector3D H_New; H_New = H * std::cos(delt_angle) + U * std::sin(delt_angle);
                    Vector3D U_New; U_New = -H * std::sin(delt_angle) + U * std::cos(delt_angle);
                    H = H_New;
                    U = U_New;
                } else if (c == '&') {
                    Vector3D H_New; H_New = H * std::cos(-delt_angle) + U * std::sin(-delt_angle);
                    Vector3D U_New; U_New = -H * std::sin(-delt_angle) + U * std::cos(-delt_angle);
                    H = H_New;
                    U = U_New;
                } else if (c == '\\') {
                    Vector3D L_New; L_New = L * std::cos(delt_angle) - U * std::sin(delt_angle);
                    Vector3D U_New; U_New = L * std::sin(delt_angle) + U * std::cos(delt_angle);
                    L = L_New;
                    U = U_New;
                } else if (c == '/') {
                    Vector3D L_New; L_New = L * std::cos(-delt_angle) - U * std::sin(-delt_angle);
                    Vector3D U_New; U_New = L * std::sin(-delt_angle) + U * std::cos(-delt_angle);
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
                    unsigned int p1Ind = figure.points.size() - 1;
                    x += H.x;
                    y += H.y;
                    z += H.z;
                    Vector3D p2 = Vector3D::point(x, y, z);
                    figure.points.push_back(p2);
                    unsigned int p2Ind = p1Ind + 1;
                    fig::Face f({p1Ind, p2Ind});
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
}