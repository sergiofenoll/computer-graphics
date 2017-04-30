#include <fstream>
#include <vector>
#include <cmath>
#include <set>
#include <stack>
#include <iostream>
#include "easy_image.hh"
#include "ini_configuration.hh"
#include "vector.hh"
#include "l_parser.hh"
#include "figures.hh"

img::EasyImage IntroColorRectangle(const ini::Configuration& configuration){
    const ini::Section imgProp = configuration["ImageProperties"];
    const unsigned int w = imgProp["width"].as_int_or_die();
    const unsigned int h = imgProp["height"].as_int_or_die();
    img::EasyImage image(w, h);
    // Iterate over pixels and colour them in
    for (unsigned int x=0; x<w; x++) {
        for (unsigned int y=0; y<h; y++) {
            image(x, y).red = (uint8_t) x;
            image(x, y).blue = (uint8_t) y;
            image(x, y).green = (uint8_t) (x + y) % 256;
        }
    }
    return image;
}

img::EasyImage IntroBlocks(const ini::Configuration& configuration){
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
                    image(x, y).red = (uint8_t) std::round(colorBlack[0] * 255);
                    image(x, y).blue = (uint8_t) std::round(colorBlack[1] * 255);
                    image(x, y).green = (uint8_t) std::round(colorBlack[2] * 255);
                } else {
                    // WIT
                    image(x, y).red = (uint8_t) std::round(colorWhite[0] * 255);
                    image(x, y).blue = (uint8_t) std::round(colorWhite[1] * 255);
                    image(x, y).green = (uint8_t) std::round(colorWhite[2] * 255);
                }
            } else {
                if ((xCoord+yCoord) % 2 == 0) {
                    // WIT
                    image(x, y).red = (uint8_t) std::round(colorWhite[0] * 255);
                    image(x, y).blue = (uint8_t) std::round(colorWhite[1] * 255);
                    image(x, y).green = (uint8_t) std::round(colorWhite[2] * 255);
                } else {
                    // ZWART
                    image(x, y).red = (uint8_t) std::round(colorBlack[0] * 255);
                    image(x, y).blue = (uint8_t) std::round(colorBlack[1] * 255);
                    image(x, y).green = (uint8_t) std::round(colorBlack[2] * 255);
                }
            }
        }
    }
    return image;
}

img::EasyImage IntroLines(const ini::Configuration& configuration){
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

    image.clear(img::Color((uint8_t) (bg[0]*255), (uint8_t) (bg[1]*255), (uint8_t) (bg[2]*255)));

    unsigned int zero = 0;
    unsigned int height = h - 1;
    unsigned int width = w - 1;

    if (figure == "QuarterCircle"){
        for(unsigned int i=0; i<nrLines; i++){
            image.draw_line(zero, startYCoor, endXCoor, height, img::Color(col[0]*255, col[1]*255, col[2]*255));
            startYCoor = startYCoor + (unsigned int) std::round(wInterval);
            endXCoor = endXCoor + (unsigned int) std::round(hInterval);
        }
        image.draw_line(zero, height, width, height, img::Color(col[0]*255, col[1]*255, col[2]*255));
    } else if(figure == "Eye"){
        for(unsigned int i=0; i<nrLines; i++){
            image.draw_line(startXCoor, zero, width, endYCoor, img::Color(col[0]*255, col[1]*255, col[2]*255));
            image.draw_line(zero, startYCoor, endXCoor, height, img::Color(col[0]*255, col[1]*255, col[2]*255));
            startXCoor = startXCoor + (unsigned int) std::round(wInterval);
            startYCoor = startYCoor + (unsigned int) std::round(wInterval);
            endXCoor = endXCoor + (unsigned int) std::round(hInterval);
            endYCoor = endYCoor + (unsigned int) std::round(hInterval);
        }
        image.draw_line(zero, height, width, height, img::Color(col[0]*255, col[1]*255, col[2]*255));
        image.draw_line(width, zero, width, height, img::Color(col[0]*255, col[1]*255, col[2]*255));
    } else if(figure == "Diamond"){
        unsigned int hMid = (unsigned int) std::round((h-1)/2);
        unsigned int wMid = (unsigned int) std::round((w-1)/2);
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
            image.draw_line(wMid, hMid, width, hMid, img::Color(col[0]*255, col[1]*255, col[2]*255));
            startX += std::round(wInterval/2);
            startYTop -= std::round(hInterval/2);
            startYBot += std::round(hInterval/2);
            endX += std::round(wInterval/2);
            endYTop += std::round(hInterval/2);
            endYBot -= std::round(hInterval/2);
        }
    }
    return image;
}

img::EasyImage LSystem(const ini::Configuration& configuration){
    ini::Section general = configuration["General"];
    std::string type = general["type"].as_string_or_die();
    const int size = general["size"].as_int_or_die();

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
    img::EasyImage image = draw2DLines(lines, size, bc, bcGr, isGrB);
    return image;
};

img::EasyImage Wireframe(const ini::Configuration& configuration){
    ini::Section general = configuration["General"];
    FigureType fig_type = Wires;

    // ##### ZBuffered Wireframe #####
    if (general["type"].as_string_or_die() == "ZBufferedWireframe") fig_type = ZBuff;

    // ##### ZBuffering with triangles #####
    if (general["type"].as_string_or_die() == "ZBuffering") fig_type = Trian;

    // ##### Lighted ZBuffering #####
    if (general["type"].as_string_or_die() == "LightedZBuffering") fig_type = Trian;

    // ##### Image size #####
    const unsigned int size = (const unsigned int) general["size"].as_double_or_die();

    // ##### Background colour #####
    const std::vector<double> bcVector = general["backgroundcolor"].as_double_tuple_or_die();
    img::Color bc((uint8_t) (bcVector[0] * 255), (uint8_t) (bcVector[1] * 255), (uint8_t) (bcVector[2] * 255));

    // ##### Background gradient #####
    std::vector<double> bcVectorGr = {0, 0, 0};
    bool isGrB = general["gradient"].as_double_tuple_if_exists(bcVectorGr);
    img::Color bcGr((uint8_t) (bcVectorGr[0]*255), (uint8_t) (bcVectorGr[1]*255), (uint8_t) (bcVectorGr[2]*255));

    const unsigned int nrFig = (const unsigned int) general["nrFigures"].as_int_or_die();

    std::vector<double> e = general["eye"].as_double_tuple_or_die();
    Vector3D eye = Vector3D::point(e[0], e[1], e[2]);

    // ##### Lights #####
    int nrLights = 0; general["nrLights"].as_int_if_exists(nrLights);
    Lights lights;
    for (unsigned int i = 0; i < nrLights; i++) {
        ini::Section light = configuration["Light" + std::to_string(i)];
        bool infinity;
        bool isDiffuse = light["infinity"].as_bool_if_exists(infinity);
        Light l;
        if (isDiffuse) {
            if (infinity) {
                std::vector<double> dir; light["direction"].as_double_tuple_if_exists(dir);
                Vector3D direction = Vector3D::vector(dir[0], dir[1], dir[2]);
                l.isAmbient = light["ambientLight"].as_double_tuple_if_exists(l.ambientLight);
                l.isDiffuseInf = light["diffuseLight"].as_double_tuple_if_exists(l.diffuseLight);
                l.isSpecular = light["specularLight"].as_double_tuple_if_exists(l.specularLight);
                l.isDiffusePnt = false;
                l.lightDirection = direction;
                lights.push_back(l);
            }
            else {
                std::vector<double> loc; light["location"].as_double_tuple_if_exists(loc);
                Vector3D location = Vector3D::point(loc[0], loc[1], loc[2]);
                l.isAmbient = light["ambientLight"].as_double_tuple_if_exists(l.ambientLight);
                l.isDiffusePnt = light["diffuseLight"].as_double_tuple_if_exists(l.diffuseLight);
                l.isSpecular = light["specularLight"].as_double_tuple_if_exists(l.specularLight);
                l.isDiffuseInf = false;
                l.location = location;
                lights.push_back(l);
            }
        }
        else {
            l.isAmbient = light["ambientLight"].as_double_tuple_if_exists(l.ambientLight);
            l.isDiffuseInf = false;
            l.isDiffusePnt = false;
            // light["diffuseLight"].as_double_tuple_if_exists(l.diffuseLight);
            l.isSpecular = light["specularLight"].as_double_tuple_if_exists(l.specularLight);
            lights.push_back(l);
        }
    }

    Figures3D figures;
    for (unsigned int i = 0; i < nrFig; i++) {
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
        std::vector<double> color;
        bool isCol = figure["color"].as_double_tuple_if_exists(color);

        // ##### Figure gradient #####
        std::vector<double> colGr = {0, 0, 0};
        bool isGr = figure["gradient"].as_double_tuple_if_exists(colGr);
        img::Color colorGr((uint8_t) (colGr[0]*255), (uint8_t) (colGr[1]*255), (uint8_t) (colGr[2]*255));

        Figure* fig;
        bool isFract = false;
        double fractalScale;
        unsigned int max_iter;
        if (type=="LineDrawing"){
            // fig = createLineDrawing(configuration, i);
        }
        else if (type == "Cube"){
            fig = new Cube();
        }
        else if (type=="Tetrahedron"){
            fig = new Tetrahedron();
        }
        else if (type=="Octahedron"){
            fig = new Octahedron();
        }
        else if (type=="Icosahedron"){
            fig = new Icosahedron();
        }
        else if (type=="Dodecahedron"){
            fig = new Dodecahedron();
        }
        else if (type=="Icosahedron"){
            fig = new Icosahedron();
        }
        else if (type=="Dodecahedron"){
            fig = new Dodecahedron();
        }
        else if (type=="Cylinder"){
            double h = figure["height"].as_double_or_die();
            const unsigned int n = (const unsigned int) figure["n"].as_int_or_die();
            fig = new Cylinder(h, n);
        }
        else if (type=="Cone"){
            double h = figure["height"].as_double_or_die();
            const unsigned int n = (const unsigned int) figure["n"].as_int_or_die();
            fig = new Cone(h, n);
        }
        else if(type=="Sphere"){
            const unsigned int n = (const unsigned int) figure["n"].as_int_or_die();
            double radius;
            fig = new Sphere(radius, n);
        }
        else if(type=="Torus"){
            double r = figure["r"].as_double_or_die();
            double R = figure["R"].as_double_or_die();
            const unsigned int m = (const unsigned int) figure["m"].as_int_or_die();
            const unsigned int n = (const unsigned int) figure["n"].as_int_or_die();
            fig = new Torus(r, R, m, n);
        }
        else if(type=="SpecialSphere"){
            double r = figure["r"].as_double_or_die();
            double R = figure["R"].as_double_or_die();
            const unsigned int m = (const unsigned int) figure["m"].as_int_or_die();
            const unsigned int n = (const unsigned int) figure["n"].as_int_or_die();
            // fig = createSpecialSphere(r, R, m, n);
        }
        else if(type=="3DLSystem") {
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
                                 delt_angle, (*fig),
                                 doubleStack, vectorStack,
                                 x, y, z, H, L, U);
        }
        else if (type == "FractalCube") {
            fractalScale = figure["fractalScale"].as_double_or_die();
            max_iter = (unsigned int) figure["nrIterations"].as_int_or_die();
            fig = new Cube();
            isFract = true;
        }
        else if (type == "FractalDodecahedron") { 
            fractalScale = figure["fractalScale"].as_double_or_die();
            max_iter = (unsigned int) figure["nrIterations"].as_int_or_die();
            fig = new Dodecahedron();
            isFract = true;
        }
        else if (type == "FractalIcosahedron") {
            fractalScale = figure["fractalScale"].as_double_or_die();
            max_iter = (unsigned int) figure["nrIterations"].as_int_or_die();
            fig = new Icosahedron();
            isFract = true;
        }
        else if (type == "FractalOctahedron") {
            fractalScale = figure["fractalScale"].as_double_or_die();
            max_iter = (unsigned int) figure["nrIterations"].as_int_or_die();
            fig = new Octahedron();
            isFract = true;
        }
        else if (type == "FractalTetrahedron") {
            fractalScale = figure["fractalScale"].as_double_or_die();
            max_iter = (unsigned int) figure["nrIterations"].as_int_or_die();
            fig = new Tetrahedron();
            isFract = true;
        }
        else if (type == "BuckyBall") {
            fig = new Buckyball();
        }
        else if (type == "FractalBuckyBall") {
            fractalScale = figure["fractalScale"].as_double_or_die();
            max_iter = (unsigned int) figure["nrIterations"].as_int_or_die();
            fig = new Buckyball();
            isFract = true;
        }
        else if (type == "MengerSponge") {
            max_iter = (unsigned int) figure["nrIterations"].as_int_or_die();
            fig = new Cube();
            // fig.generateMengerSponge(fig, max_iter);
        }
        else {
            std::cerr << "ERROR: Unknown figure.\n";
        }

        figure["color"].as_fig_color_if_exists(fig->ambient_reflection);
        figure["ambient_reflection"].as_fig_color_if_exists(fig->ambient_reflection);
        figure["diffuseReflection"].as_fig_color_if_exists(fig->diffuse_reflection);
        figure["specular_reflection"].as_fig_color_if_exists(fig->specular_reflection);
        figure["reflectionCoefficient"].as_double_if_exists(fig->reflection_coefficient);
        // fig.colorGr = colorGr;
        fig->isGradient = isGr;
        Matrix sclM = scaleFigure(scale);
        Matrix rotX = rotateX(rotXang);
        Matrix rotY = rotateY(rotYang);
        Matrix rotZ = rotateZ(rotZang);
        Matrix trsM = translate(center);
        Matrix eyeT = eyePointTrans(eye);
        Matrix V = sclM * rotX * rotY * rotZ * trsM  * eyeT;
        fig->apply_transformation(V);
        unsigned int cur_iter = 1;
        if (isFract) recursiveGenerateFractal((*fig), figures, fractalScale, cur_iter, max_iter);
        else figures.push_back(fig);
    }
    Matrix eyeT = eyePointTrans(eye);
    for (auto& light : lights) {
        if (light.isDiffuseInf) {
            light.lightDirection *= eyeT;
            light.location *= eyeT;
        }
        else if (light.isDiffusePnt) {
            light.location *= eyeT;
        }
    }
    return draw3DLines(figures, fig_type, size, lights, bc, bcGr, isGrB);
}

img::EasyImage generate_image(const ini::Configuration& configuration){
    std::string type = configuration["General"]["type"].as_string_or_die();
    if(type=="IntroColorRectangle") {
        return IntroColorRectangle(configuration);
    } else if(type=="IntroBlocks"){
        return IntroBlocks(configuration);
    } else if(type=="IntroLines"){
        return IntroLines(configuration);
    } else if(type=="2DLSystem" or type=="2DLSystemStoch"){
        return LSystem(configuration);
    } else if (type == "Wireframe" or 
               type == "ZBufferedWireframe" or
               type == "ZBuffering" or 
               type == "LightedZBuffering"){
        return Wireframe(configuration);
    }
    img::EasyImage image;
    return image;
}

int main(int argc, char const* argv[]){
    std::srand((unsigned int) time(NULL)); // Random seed
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
