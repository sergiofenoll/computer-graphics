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

// Helper functions

inline int roundToInt(double d){
    return d<0 ? std::ceil(d-0.5) : std::floor(d+0.5);
}

img::Color colorGradient(img::Color &a, img::Color &b, double &t){
    return img::Color(
        roundToInt(a.red + (b.red - a.red) * t),
        roundToInt(a.green + (b.green - a.green) * t),
        roundToInt(a.blue + (b.blue - a.blue) * t));
}

void gradientBg(img::EasyImage &image, img::Color &a, img::Color &b){
    for (unsigned int x=0; x<image.get_width(); x++){
        for (unsigned int y=0; y<image.get_height(); y++){
            double p = x / (image.get_width() - 1.0);
            image(x,y).red = roundToInt((1 - p) * a.red + p * b.red + 0.5);
            image(x,y).green = roundToInt((1 - p) * a.green + p * b.green + 0.5);
            image(x,y).blue = roundToInt((1 - p) * a.blue + p * b.blue + 0.5);
        }
    }
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
    ini::Section general = configuration["General"];
    std::string type = general["type"].as_string_or_die();
    const double size = general["size"].as_double_or_die();
    std::vector<double> aspectratio = {0, 1};
    general["aspectratio"].as_double_tuple_if_exists(aspectratio);
    double ar = aspectratio[0] / aspectratio[1];

    const std::vector<double> bcVector = general["backgroundcolor"].as_double_tuple_or_die();
    img::Color bc(bcVector[0]*255, bcVector[1]*255, bcVector[2]*255);
    std::vector<double> bcVectorGr = {0, 0, 0};
    bool isGrB = general["gradient"].as_double_tuple_if_exists(bcVectorGr);
    img::Color bcGr(bcVectorGr[0]*255, bcVectorGr[1]*255, bcVectorGr[2]*255);

    const std::vector<double> lcVector = configuration["2DLSystem"]["color"].as_double_tuple_or_die();
    img::Color lc(lcVector[0]*255, lcVector[1]*255, lcVector[2]*255);
    std::vector<double> lcVectorGr = {0, 0, 0};
    bool isGrL = configuration["2DLSystem"]["gradient"].as_double_tuple_if_exists(lcVectorGr);
    img::Color lcGr(lcVectorGr[0]*255, lcVectorGr[1]*255, lcVectorGr[2]*255);

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
    unsigned int max_it = l_system.get_nr_iterations();
    double alph_angle = (l_system.get_starting_angle() / 360.0) * 2*PI;
    double delt_angle = (l_system.get_angle() / 360.0) * 2*PI;

    Lines2D lines;
    std::stack<double> brStack;
    unsigned int cur_it = 0;
    double x = 0;
    double y = 0;
    recursivePrintString(l_system, init, cur_it, max_it,
                         alph_angle, delt_angle, lines, brStack,
                         lc, lcGr, isGrL, x, y);
    img::EasyImage image = draw2DLines(lines, size, bc, bcGr, isGrB, ar);
    return image;
};

img::EasyImage Wireframe(const ini::Configuration &configuration){
    ini::Section general = configuration["General"];

    // ##### ZBuffered Wireframe #####
    bool isZBuf = false;
    if (general["type"].as_string_or_die() == "ZBufferedWireframe") isZBuf = true;

    // ##### ZBuffering with triangles #####
    bool isTriang = false;
    if (general["type"].as_string_or_die() == "ZBuffering") isTriang = true;

    // ##### Image size #####
    double size = general["size"].as_double_or_die();

    // ##### Aspect ratio #####
    std::vector<double> aspectratio = {0, 1};
    general["aspectratio"].as_double_tuple_if_exists(aspectratio);
    double ar = aspectratio[0] / aspectratio[1];

    // ##### Background colour #####
    const std::vector<double> bcVector = general["backgroundcolor"].as_double_tuple_or_die();
    img::Color bc(bcVector[0]*255, bcVector[1]*255, bcVector[2]*255);

    // ##### Background gradient #####
    std::vector<double> bcVectorGr = {0, 0, 0};
    bool isGrB = general["gradient"].as_double_tuple_if_exists(bcVectorGr);
    img::Color bcGr(bcVectorGr[0]*255, bcVectorGr[1]*255, bcVectorGr[2]*255);

    int nrFig = general["nrFigures"].as_int_or_die();

    std::vector<double> e = general["eye"].as_double_tuple_or_die();
    Vector3D eye = Vector3D::point(e[0], e[1], e[2]);

    Figures3D figures;
    for (int i=0; i<nrFig; i++){
        ini::Section figure = configuration["Figure" + std::to_string(i)];
        std::string type = figure["type"].as_string_or_die();

        // ##### Rotation matrices #####
        double rotXang = (figure["rotateX"].as_double_or_die() / 360.0) * 2*PI;
        double rotYang = (figure["rotateY"].as_double_or_die() / 360.0) * 2*PI;
        double rotZang = (figure["rotateZ"].as_double_or_die() / 360.0) * 2*PI;
        double scale = figure["scale"].as_double_or_die();

        // ##### Figure center #####
        std::vector<double> cen = figure["center"].as_double_tuple_or_die();
        Vector3D center = Vector3D::vector(cen[0], cen[1], cen[2]);

        // ##### Figure colour #####
        std::vector<double> col = figure["color"].as_double_tuple_or_die();
        img::Color color(col[0]*255, col[1]*255, col[2]*255);

        // ##### Figure gradient #####
        std::vector<double> colGr = {0, 0, 0};
        bool isGr = figure["gradient"].as_double_tuple_if_exists(colGr);
        img::Color colorGr(colGr[0]*255, colGr[1]*255, colGr[2]*255);

        Figure fig;
        if (type=="LineDrawing"){
            fig = createLineDrawing(configuration, i);
        } else if (type == "Cube"){
            fig = createCube();
        } else if (type=="Tetrahedron"){
            fig = createTetrahedron();
        } else if (type=="Octahedron"){
            fig = createOctahedron();
        } else if (type=="Icosahedron"){
            fig = createIcosahedron();
        } else if (type=="Dodecahedron"){
            fig = createDodecahedron();
        } else if (type=="Icosahedron"){
            fig = createIcosahedron();
        } else if (type=="Dodecahedron"){
            fig = createDodecahedron();
        } else if (type=="Cylinder"){
            double h = figure["height"].as_double_or_die();
            int n = figure["n"].as_int_or_die();
            fig = createCylinder(h, n);
        } else if (type=="Cone"){
            double h = figure["height"].as_double_or_die();
            int n = figure["n"].as_int_or_die();
            fig = createCone(h, n);
        } else if(type=="Sphere"){
            int n = figure["n"].as_int_or_die();
            double radius;
            fig = createSphere(radius, n);
        } else if(type=="Torus"){
            double r = figure["r"].as_double_or_die();
            double R = figure["R"].as_double_or_die();
            int m = figure["m"].as_int_or_die();
            int n = figure["n"].as_int_or_die();
            fig = createTorus(r, R, m, n);
        } else if(type=="SpecialSphere"){
            double r = figure["r"].as_double_or_die();
            double R = figure["R"].as_double_or_die();
            int m = figure["m"].as_int_or_die();
            int n = figure["n"].as_int_or_die();
            fig = createSpecialSphere(r, R, m, n);
        } else if(type=="3DLSystem") {
            std::string inf = figure["inputfile"].as_string_or_die();
            // Create 2D L-system
            LParser::LSystem3D l_system;
            std::ifstream input_stream(inf);
            if (type=="3DLSystemStoch") l_system.setStoch(true);
            input_stream >> l_system;
            input_stream.close();
            // Get values from 3D L-system
            std::set<char> alph = l_system.get_alphabet();
            std::string init = l_system.get_initiator();
            unsigned int max_it = l_system.get_nr_iterations();
            double delt_angle = (l_system.get_angle() / 360.0) * 2*PI;

            std::stack<double> doubleStack;
            std::stack<Vector3D> vectorStack;
            unsigned int cur_it = 0;
            double x = 0;
            double y = 0;
            double z = 0;
            Vector3D H = Vector3D::vector(1, 0, 0);
            Vector3D L = Vector3D::vector(0, 1, 0);
            Vector3D U = Vector3D::vector(0, 0, 1);

            recursivePrintString(l_system, init, cur_it, max_it,
                                 delt_angle, fig,
                                 doubleStack, vectorStack,
                                 x, y, z, H, L, U);
        }
        fig.color = color;
        fig.colorGr = colorGr;
        fig.isGradient = isGr;
        Matrix sclM = scaleFigure(scale);
        Matrix rotX = rotateX(rotXang);
        Matrix rotY = rotateY(rotYang);
        Matrix rotZ = rotateZ(rotZang);
        Matrix trsM = translate(center);
        Matrix eyeT = eyePointTrans(eye);
        Matrix V = sclM * rotX * rotY * rotZ * trsM  * eyeT;
        applyTransformation(fig, V);
        figures.push_back(fig);
    }
    Lines2D lines;
    lines = doProjection(figures);
    if (isZBuf) {
        return drawZBufLines(lines, size, bc, bcGr, isGrB, ar);
    } else if (isTriang) {
        return drawTriangLines(figures, size, bc, bcGr, isGrB, ar);
    } else {
        return draw2DLines(lines, size, bc, bcGr, isGrB, ar);
    }
}

img::EasyImage generate_image(const ini::Configuration &configuration){
    std::string type = configuration["General"]["type"].as_string_or_die();
    if(type=="IntroColorRectangle") {
        return IntroColorRectangle(configuration);
    } else if(type=="IntroBlocks"){
        return IntroBlocks(configuration);
    } else if(type=="IntroLines"){
        return IntroLines(configuration);
    } else if(type=="2DLSystem" or type=="2DLSystemStoch"){
        return L2DSystem(configuration);
    } else if (type=="Wireframe" or type=="ZBufferedWireframe" or type == "ZBuffering"){
        return Wireframe(configuration);
    }
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
