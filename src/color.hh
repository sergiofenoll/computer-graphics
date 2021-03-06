//
// Created by sergio on 30/04/17.
//

#ifndef COMPGRAPHX_COLOR_HH
#define COMPGRAPHX_COLOR_HH

#include <vector>
#include <map>
#include <set>
#include <assert.h>
#include "easy_image.hh"
#include "vector.hh"

namespace col {

    /**
     * @brief A class that stores RGB color values.
     *
     * The values are stored as doubles between 0 and 1.
     */
    class Color {
    private:
        /**
         * @brief The intensity of the red color component.
         */
        double red;

        /**
         * @brief The intensity of the green color component.
         */
        double green;

        /**
         * @brief The intensity of the blue color component.
         */
        double blue;
    public:
        /**
         * @brief Default constructor
         */
        Color();

        /**
         * @brief Construct a color with the given intensities.
         * @param red The red color component.
         * @param green The green color component.
         * @param blue The blue color component.
         */
        Color(const double& red, const double& green, const double& blue);

        /**
         * @brief Equality operator overload by comparing the red, green and blue intensities of this color and @param rhs.
         * @param rhs The color that is compared to this.
         * @return true if both colors are the same, false otherwise.
         */
        bool operator==(const Color& rhs) const;

        /**
         * @brief Inequality operator overload by comparing the red, green and blue intensities of this color and @param rhs.
         * @param rhs The color that is compared to this.
         * @return true if both colors are different, false otherwise.
         */
        bool operator!=(const Color& rhs) const;

        /**
         * @brief Assignment operator overload for a Color object.
         * @param rhs The Color that is copied.
         * @return A reference to this color.
         */
        Color& operator=(const Color& rhs);

        /**
         * @brief Assignment operator overload for a vector of doubles.
         * @param rhs The vector of size 3 containing doubles whose values will be used as red, green and blue intensities.
         * @return A reference to this color.
         */
        Color& operator=(const std::vector<double>& rhs);

        Color& operator+=(const Color& rhs);

        /**
         * @brief Addition operator overload by adding the red, green and blue intensities of this color and @param rhs.
         * @param rhs The color that will be added to this color.
         * @return A reference to this color.
         */
        friend Color operator+(Color lhs, const Color& rhs);

        Color& operator*=(const Color& rhs);

        friend Color operator*(Color lhs, const Color& rhs);

        Color& operator*=(const double& rhs);

        /**
         * @brief Mulitplication operator overload by mulitplying the red, green and blue intensities by @param rhs.
         * @param rhs The double that will be mulitplied to this color.
         * @return A reference to this color.
         */
        friend Color operator*(Color lhs, const double& rhs);

        void set_red_value(const double& red);

        double get_red_value() const;

        uint8_t red_int_value() const;

        void set_green_value(const double& green);

        double get_green_value() const;

        uint8_t green_int_value() const;

        void set_blue_value(const double& blue);

        double get_blue_value() const;

        uint8_t blue_int_value() const;
    };

    class ZBuffer {
    private:
        std::vector<std::vector<double>> z_buf;
    public:
        ZBuffer();
        ZBuffer(const unsigned int& width, const unsigned int& height);
        double& operator()(unsigned int x, unsigned int y);
    };

    class OpacityMatrix {
    // tfw std::vector<bool> is stupid and doesn't work as you'd expect
    private:
        std::vector<std::vector<double>> matrix;
    public:
        OpacityMatrix();
        OpacityMatrix(const unsigned int& width, const unsigned int& height);
        double& operator()(unsigned int x, unsigned int y);
    };

    enum LightType {
        Normal,
        Shadow
    };

    class Light {
    private:
        Color ambient_light;
        Color diffuse_light;
        Color specular_light;
    public:
        Light();

        Color& ambient();

        Color& diffuse();

        Color& specular();

        virtual void set_direction(Vector3D& direction);

        virtual Vector3D get_direction();

        virtual void set_location(Vector3D& location);

        virtual Vector3D get_location();

        virtual void set_mask_size(int mask_size);

        virtual int get_mask_size();

        virtual void set_shadow_mask(ZBuffer& shadow_mask);

        virtual ZBuffer get_shadow_mask();

        virtual void set_shadow_mask_at(unsigned int x, unsigned int y, double val);

        virtual double get_shadow_mask_at(unsigned int x, unsigned int y);

        virtual void set_eye(Matrix& eye);

        virtual Matrix get_eye();

        virtual void set_values(double& d, double& dx, double& dy);

        virtual void get_values(double& d, double& dx, double& dy);

        bool is_ambient();

        virtual bool is_diffuse_inf();

        virtual bool is_diffuse_pnt();

        bool is_specular();

        virtual bool is_shadow();
    };

    class InfLight : public Light {
    private:
        Vector3D direction_vector;
    public:
        virtual void set_direction(Vector3D& direction);

        virtual Vector3D get_direction();

        virtual bool is_diffuse_inf();

        virtual bool is_diffuse_pnt();
    };

    class PntLight : public Light {
    private:
        Vector3D location_vector;

        int mask_size;

        ZBuffer shadow_mask;

        Matrix eye;

        double d, dx, dy;
    public:
        virtual void set_location(Vector3D& location);

        virtual Vector3D get_location();

        virtual void set_mask_size(int mask_size);

        virtual int get_mask_size();

        virtual void set_shadow_mask(ZBuffer& shadow_mask);

        virtual ZBuffer get_shadow_mask();

        virtual void set_shadow_mask_at(unsigned int x, unsigned int y, double val);

        virtual double get_shadow_mask_at(unsigned int x, unsigned int y);

        virtual void set_eye(Matrix& eye);

        virtual Matrix get_eye();

        virtual void set_values(double& d, double& dx, double& dy);

        virtual void get_values(double& d, double& dx, double& dy);

        virtual bool is_diffuse_inf();

        virtual bool is_diffuse_pnt();

        virtual bool is_shadow();
    };

    typedef std::map<LightType, std::vector<Light*>> Lights;
};


#endif //COMPGRAPHX_COLOR_HH
