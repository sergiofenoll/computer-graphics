#include "easy_image.hh"
#include "ini_configuration.hh"
#include "l_parser.hh"
#include "Line2D.hh"
#include "3DFigures.hh"

#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <stdlib.h>
#include <time.h>
#include <stack>

const double PI = std::atan(1.0) * 4;

// Helper functions

inline int roundToInt(double d){
    return d<0 ? (int)std::ceil(d-0.5) : (int)std::floor(d+0.5);
}

img::EasyImage draw2DLines(Lines2D &lines, const int &size, const img::Color &bc=img::Color(255, 255, 255)){
    double maxX = lines[0].p1.x;
    double maxY = lines[0].p1.y;
    double minX = lines[0].p1.x;
    double minY = lines[0].p1.y;
    // Iterate over all lines looking for xmax and ymax
    for(Line2D &l : lines){
        if(l.p1.x>maxX) maxX = l.p1.x;
        if(l.p1.x<minX) minX = l.p1.x;
        if(l.p2.x>maxX) maxX = l.p2.x;
        if(l.p2.x<minX) minX = l.p2.x;
        if(l.p1.y>maxY) maxY = l.p1.y;
        if(l.p1.y<minY) minY = l.p1.y;
        if(l.p2.y>maxY) maxY = l.p2.y;
        if(l.p2.y<minY) minY = l.p2.y;
    }
    // Ranges
    double xRange = maxX - minX;
    double yRange = maxY - minY;
    // if ((int)xRange==0) xRange = 1;
    // if ((int)yRange==0) yRange = 1;
    // Size of image (width, height)
    double maxXY = std::max(xRange, yRange);
    double imageX = size * (xRange / maxXY);
    double imageY = size * (yRange / maxXY);
    // Scaling factor d
    double d = 0.95 * (imageX / xRange);
    // Centers
    double DCx = d * ((minX + maxX) / 2);
    double DCy = d * ((minY + maxY) / 2);
    double dx = (imageX / 2) - DCx;
    double dy = (imageY / 2) - DCy;
    // Multply every coordinate with d,
    // add dx and dy
    // and draw lines
    img::EasyImage image(roundToInt(imageX), roundToInt(imageY), bc);
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
        unsigned int x0 = roundToInt(l.p1.x);
        unsigned int y0 = roundToInt(l.p1.y);
        unsigned int x1 = roundToInt(l.p2.x);
        unsigned int y1 = roundToInt(l.p2.y);
        image.draw_line(x0, y0, x1, y1, l.color);
    }
    return image;
}

Point2D doProjection(const Vector3D &point, const double &d=1.0){
    Point2D p;
    p.x = (d * point.x) / -(point.z);
    p.y = (d * point.y) / -(point.z);
    return p;
}

Lines2D doProjection(Figures3D &figures){
    Lines2D lines;
    for (Figure &fig : figures){
        for (Face &face : fig.faces){
            Point2D p1;
            Point2D p2;
            p1 = doProjection(fig.points[face.point_indexes[0]]);
            p2 = doProjection(fig.points[face.point_indexes[1]]);
            Line2D line;
            line.p1 = p1;
            line.p2 = p2;
            line.color = fig.color;
            lines.push_back(line);
        }
    }
    return lines;
}

// Image generating functions

img::EasyImage IntroColorRectangle(const ini::Configuration &configuration){
    const ini::Section imgProp = configuration["ImageProperties"];
    const unsigned int w = imgProp["width"].as_int_or_die();
    const unsigned int h = imgProp["height"].as_int_or_die();
    img::EasyImage image(w, h);
    // Iterate over pixels and colour them in
    for (unsigned int x=0; x<w; x++) {
        for (unsigned int y=0; y<h; y++) {
            image(x, y).red = x;
            image(x, y).blue = y;
            image(x, y).green = (x + y) % 256;
        }
    }
    return image;
}

img::EasyImage IntroBlocks(const ini::Configuration &configuration){
    const ini::Section imgProp = configuration["ImageProperties"];
    const unsigned int w = imgProp["width"].as_int_or_die();
    const unsigned int h = imgProp["height"].as_int_or_die();
    img::EasyImage image(w, h);
    const ini::Section blockProp = configuration["BlockProperties"];
    const unsigned int nrXBlocks = blockProp["nrXBlocks"].as_int_or_die();
    const unsigned int nrYBlocks = blockProp["nrYBlocks"].as_int_or_die();
    std::vector<double> colorWhite = blockProp["colorWhite"].as_double_tuple_or_die();
    std::vector<double> colorBlack = blockProp["colorBlack"].as_double_tuple_or_die();
    const bool invertColors = blockProp["invertColors"].as_bool_or_die();
    const unsigned int wBlock = w / nrXBlocks;
    const unsigned int hBlock = h / nrYBlocks;
    for(unsigned int x=0; x<w; x++){
        for(unsigned int y=0; y<h; y++){
            const unsigned int xCoord = x / wBlock;
            const unsigned int yCoord = y / hBlock;
            if(invertColors) {
                if ((xCoord+yCoord) % 2 == 0) {
                    // ZWART
                    image(x, y).red = roundToInt(colorBlack[0] * 255);
                    image(x, y).blue = roundToInt(colorBlack[1] * 255);
                    image(x, y).green = roundToInt(colorBlack[2] * 255);
                } else {
                    // WIT
                    image(x, y).red = roundToInt(colorWhite[0] * 255);
                    image(x, y).blue = roundToInt(colorWhite[1] * 255);
                    image(x, y).green = roundToInt(colorWhite[2] * 255);
                }
            } else {
                if ((xCoord+yCoord) % 2 == 0) {
                    // WIT
                    image(x, y).red = roundToInt(colorWhite[0] * 255);
                    image(x, y).blue = roundToInt(colorWhite[1] * 255);
                    image(x, y).green = roundToInt(colorWhite[2] * 255);
                } else {
                    // ZWART
                    image(x, y).red = roundToInt(colorBlack[0] * 255);
                    image(x, y).blue = roundToInt(colorBlack[1] * 255);
                    image(x, y).green = roundToInt(colorBlack[2] * 255);
                }
            }
        }
    }
    return image;
}

img::EasyImage IntroLines(const ini::Configuration &configuration){
    const ini::Section imgProp = configuration["ImageProperties"];
    const unsigned int w = imgProp["width"].as_int_or_die();
    const unsigned int h = imgProp["height"].as_int_or_die();
    img::EasyImage image(w, h);
    const ini::Section lineProp = configuration["LineProperties"];
    const std::string figure = lineProp["figure"].as_string_or_die();
    const std::vector<double> bg = lineProp["backgroundcolor"].as_double_tuple_or_die();
    const std::vector<double> col = lineProp["lineColor"].as_double_tuple_or_die();
    const unsigned int nrLines = lineProp["nrLines"].as_int_or_die();
    unsigned int startXCoor = 0;
    unsigned int startYCoor = 0;
    unsigned int endXCoor = 0;
    unsigned int endYCoor = 0;
    double hInterval = (double) h /(double) (nrLines-1);
    double wInterval = (double) w / (double) (nrLines-1);

    image.clear(img::Color(bg[0]*255, bg[1]*255, bg[2]*255));

    if (figure == "QuarterCircle"){
        for(unsigned int i=0; i<nrLines; i++){
            image.draw_line(0, startYCoor, endXCoor, h-1, img::Color(col[0]*255, col[1]*255, col[2]*255));
            startYCoor = startYCoor + roundToInt(wInterval);
            endXCoor = endXCoor + roundToInt(hInterval);
        }
        image.draw_line(0, h-1, w-1, h-1, img::Color(col[0]*255, col[1]*255, col[2]*255));
    } else if(figure == "Eye"){
        for(unsigned int i=0; i<nrLines; i++){
            image.draw_line(startXCoor, 0, w-1, endYCoor, img::Color(col[0]*255, col[1]*255, col[2]*255));
            image.draw_line(0, startYCoor, endXCoor, h-1, img::Color(col[0]*255, col[1]*255, col[2]*255));
            startXCoor = startXCoor + roundToInt(wInterval);
            startYCoor = startYCoor + roundToInt(wInterval);
            endXCoor = endXCoor + roundToInt(hInterval);
            endYCoor = endYCoor + roundToInt(hInterval);
        }
        image.draw_line(0, h-1, w-1, h-1, img::Color(col[0]*255, col[1]*255, col[2]*255));
        image.draw_line(w-1, 0, w-1, h-1, img::Color(col[0]*255, col[1]*255, col[2]*255));
    } else if(figure == "Diamond"){
        unsigned int hMid = roundToInt((h-1)/2);
        unsigned int wMid = roundToInt((w-1)/2);
        unsigned int startX = startXCoor;
        unsigned int startYTop = h-1;
        unsigned int startYBot = 0;
        unsigned int endX = wMid;
        unsigned int endYTop = hMid;
        unsigned int endYBot = hMid;
        for(unsigned int i=0; i<nrLines; i++){
            image.draw_line(startX, hMid, wMid, endYTop, img::Color(col[0]*255, col[1]*255, col[2]*255));
            image.draw_line(startX, hMid, wMid, endYBot, img::Color(col[0]*255, col[1]*255, col[2]*255));
            image.draw_line(wMid, startYTop, endX, hMid, img::Color(col[0]*255, col[1]*255, col[2]*255));
            image.draw_line(wMid, startYBot, endX, hMid, img::Color(col[0]*255, col[1]*255, col[2]*255));
            image.draw_line(wMid, hMid, w-1, hMid, img::Color(col[0]*255, col[1]*255, col[2]*255));
            startX += roundToInt(wInterval/2);
            startYTop -= roundToInt(hInterval/2);
            startYBot += roundToInt(hInterval/2);
            endX += roundToInt(wInterval/2);
            endYTop += roundToInt(hInterval/2);
            endYBot -= roundToInt(hInterval/2);
        }
    }
    return image;
}

img::EasyImage L2DSystem(const ini::Configuration &configuration){
    // TODO: Implement stochastic L systems
    std::string type = configuration["General"]["type"].as_string_or_die();
    const double size = configuration["General"]["size"].as_double_or_die();
    const std::vector<double> bcVector = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
    img::Color bc(bcVector[0]*255, bcVector[1]*255, bcVector[2]*255);
    const std::vector<double> lcVector = configuration["2DLSystem"]["color"].as_double_tuple_or_die();
    img::Color lc(lcVector[0]*255, lcVector[1]*255, lcVector[2]*255);
    std::string inf = configuration["2DLSystem"]["inputfile"];
    // Create 2D L-system
    LParser::LSystem2D l_system;
    std::ifstream input_stream(inf);
    if (type=="2DLSystemStoch") l_system.setStoch(true);
    input_stream >> l_system;
    input_stream.close();
    // Get values from 2D L-system
    std::set<char> alph = l_system.get_alphabet();
    std::string init = l_system.get_initiator();
    const unsigned int iter = l_system.get_nr_iterations();
    double alph_angle = (l_system.get_starting_angle() / 360.0) * 2*PI;
    const double delt_angle = (l_system.get_angle() / 360.0) * 2*PI;
    // Create string
    std::string print_string = init;
    std::string temp_string = "";
    for(uint8_t i=0; i<iter; i++){
        temp_string = "";
        for(const char c : print_string){
            if(alph.find(c)!=alph.end()){
                temp_string += l_system.get_replacement(c);
            } else if (c=='+'||c=='-'||c=='('||c==')') {
                temp_string += c;
            }
        }
        print_string = temp_string;
    }
    temp_string = "";
    // Draw lines based on print_string
    Lines2D lines;
    double x = 0; // Starting x coordinate
    double y = 0; // Starting y coordinate
    // For brackets
    std::stack<double> brStack; // Brackets Stack
    for (char c : print_string) {
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
            Line2D line(p1, p2, lc);
            lines.push_back(line);
        } else {
            x += std::cos(alph_angle);
            y += std::sin(alph_angle);
        }
    }
    img::EasyImage image = draw2DLines(lines, size, bc);
    return image;
};

img::EasyImage Wireframe(const ini::Configuration &configuration){
    ini::Section general = configuration["General"];
    double size = general["size"].as_double_or_die();
    std::vector<double> bcDub = general["backgroundcolor"].as_double_tuple_or_die();
    img::Color bc(bcDub[0]*255, bcDub[1]*255, bcDub[2]*255);
    int nrFig = general["nrFigures"].as_int_or_die();
    std::vector<double> e = general["eye"].as_double_tuple_or_die();
    Vector3D eye = Vector3D::point(e[0], e[1], e[2]);
    // Vector3D trans = Vector3D::vector(0, 0, -e[2]);
    Figures3D figures;
    for (int i=0; i<nrFig; i++){
        ini::Section figure = configuration["Figure" + std::to_string(i)];
        std::string type = figure["type"].as_string_or_die();
        double rotXang = (figure["rotateX"].as_double_or_die() / 360.0) * 2*PI;
        double rotYang = (figure["rotateY"].as_double_or_die() / 360.0) * 2*PI;
        double rotZang = (figure["rotateZ"].as_double_or_die() / 360.0) * 2*PI;
        double scale = figure["scale"].as_double_or_die();
        std::vector<double> cen = figure["center"].as_double_tuple_or_die();
        Vector3D center = Vector3D::vector(cen[0], cen[1], cen[2]);
        std::vector<double> col = figure["color"].as_double_tuple_or_die();
        img::Color color(col[0]*255, col[1]*255, col[2]*255);
        int nrPoints = figure["nrPoints"].as_int_or_die();
        int nrLines = figure["nrLines"].as_int_or_die();

        std::vector<Vector3D> points;
        std::vector<Face> faces;

        for (int j=0; j<nrPoints; j++){
            std::vector<double> p = figure["point" + std::to_string(j)].as_double_tuple_or_die();
            Vector3D point = Vector3D::point(p[0], p[1], p[2]);
            points.push_back(point);
        }

        for (int j=0; j<nrLines; j++){
            std::vector<int> p = figure["line" + std::to_string(j)].as_int_tuple_or_die();
            Face f;
            for (int k : p){
                f.point_indexes.push_back(k);
            }
            faces.push_back(f);
        }
        Matrix sclM = scaleFigure(scale);
        Matrix rotX = rotateX(rotXang);
        Matrix rotY = rotateY(rotYang);
        Matrix rotZ = rotateZ(rotZang);
        // Matrix trsM = translate(center);
        Matrix eyeT = eyePointTrans(eye);
        Matrix V = sclM * rotX * rotY * rotZ * /* trsM * */ eyeT;
        Figure fig;
        fig.points = points;
        fig.faces = faces;
        fig.color = color;
        applyTransformation(fig, V);
        figures.push_back(fig);
    }
    Lines2D lines;
    lines = doProjection(figures);
    return draw2DLines(lines, size, bc);
}

img::EasyImage generate_image(const ini::Configuration &configuration){
    std::string type = configuration["General"]["type"].as_string_or_die();
    if(type=="IntroColorRectangle") {
        return IntroColorRectangle(configuration);
    } else if(type=="IntroBlocks"){
        return IntroBlocks(configuration);
    } else if(type=="IntroLines"){
        return IntroLines(configuration);
    } else if(type=="2DLSystem"||type=="2DLSystemStoch"){
        return L2DSystem(configuration);
    } else if (type=="Wireframe")
        return Wireframe(configuration);
    img::EasyImage image;
    return image;
}

int main(int argc, char const* argv[]){
    std::srand(time(NULL)); // Random seed
    int retVal = 0;
    try{
        for(int i = 1; i < argc; ++i){
            ini::Configuration conf;
            try{
                std::ifstream fin(argv[i]);
                fin >> conf;
                fin.close();
            } catch(ini::ParseException& ex){
                std::cerr << "Error parsing file: " << argv[i] << ": " << ex.what() << std::endl;
                retVal = 1;
                continue;
            }
            img::EasyImage image = generate_image(conf);
                if(image.get_height() > 0 && image.get_width() > 0){
                    std::string fileName(argv[i]);
                    std::string::size_type pos = fileName.rfind('.');
                    if(pos == std::string::npos){
                        //filename does not contain a '.' --> append a '.bmp' suffix
                        fileName += ".bmp";
                    } else{
                        fileName = fileName.substr(0,pos) + ".bmp";
                    } try{
                        std::ofstream f_out(fileName.c_str(),std::ios::trunc | std::ios::out | std::ios::binary);
                        f_out << image;
                    } catch(std::exception& ex){
                        std::cerr << "Failed to write image to file: " << ex.what() << std::endl;
                        retVal = 1;
                    }
                } else{
                    std::cout << "Could not generate image for " << argv[i] << std::endl;
                }
        }
    }
    catch(const std::bad_alloc &exception){
	//When you run out of memory this exception is thrown. When this happens the return value of the program MUST be '100'.
	//Basically this return value tells our automated test scripts to run your engine on a pc with more memory.
	//(Unless of course you are already consuming the maximum allowed amount of memory)
	//If your engine does NOT adhere to this requirement you risk losing points because then our scripts will
	//mark the test as failed while in reality it just needed a bit more memory
        std::cerr << "Error: insufficient memory" << std::endl;
        retVal = 100;
    }
    return retVal;
}
