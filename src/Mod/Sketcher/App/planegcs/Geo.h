/***************************************************************************
 *   Copyright (c) 2011 Konstantinos Poulios <logari81@gmail.com>          *
 *                                                                         *
 *   This file is part of the FreeCAD CAx development system.              *
 *                                                                         *
 *   This library is free software; you can redistribute it and/or         *
 *   modify it under the terms of the GNU Library General Public           *
 *   License as published by the Free Software Foundation; either          *
 *   version 2 of the License, or (at your option) any later version.      *
 *                                                                         *
 *   This library  is distributed in the hope that it will be useful,      *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU Library General Public License for more details.                  *
 *                                                                         *
 *   You should have received a copy of the GNU Library General Public     *
 *   License along with this library; see the file COPYING.LIB. If not,    *
 *   write to the Free Software Foundation, Inc., 59 Temple Place,         *
 *   Suite 330, Boston, MA  02111-1307, USA                                *
 *                                                                         *
 ***************************************************************************/

#ifndef PLANEGCS_GEO_H
#define PLANEGCS_GEO_H

#include <cmath>
#include "Util.h"

#if __cplusplus < 202002L
// see https://en.cppreference.com/w/cpp/numeric/lerp already present in c++20
namespace std
{
inline double lerp(double s, double e, double t)
{
    if (t == 0.0)
        return s;
    if (t == 1.0)
        return e;
    return std::fma(e - s, t, s);// see https://en.cppreference.com/w/cpp/numeric/math/fma
}
}// namespace std
#endif

namespace GCS
{
inline int signum(double x) { return x < 0 ? -1 : (x > 0 ? 1 : 0); }

    class Point
    {
    public:
        Point() : x(nullptr), y(nullptr) {}
        Point(double *px, double *py) : x(px), y(py) {}
        double *x;
        double *y;
    };

    using VEC_P = std::vector<Point>;

    ///Class DeriVector2 holds a vector value and its derivative on the
    ///parameter that the derivatives are being calculated for now. x,y is the
    ///actual vector (v). dx,dy is a derivative of the vector by a parameter
    ///(dv/dp). The easiest way to fill the vector in is by passing a point and
    ///a derivative parameter pointer to its constructor. x,y are read from the
    ///pointers in Point, and dx,dy are set to either 0 or 1 depending on what
    ///pointers of Point match the supplied pointer. The derivatives can be set
    ///manually as well. The class also provides a bunch of methods to do math
    ///on it (and derivatives are calculated implicitly).
    ///
    class DeriVector2
    {
    public:
        DeriVector2(){x=0; y=0; dx=0; dy=0;}
        DeriVector2(double x, double y) {this->x = x; this->y = y; this->dx = 0; this->dy = 0;}
        DeriVector2(double x, double y, double dx, double dy) {this->x = x; this->y = y; this->dx = dx; this->dy = dy;}
        DeriVector2(const Point &p, const double* derivparam);
        double x, dx;
        double y, dy;

        double length() const { return std::hypot(x, y); }
        double sqLength() const { return x * x + y * y; }
        double length(double &dlength) const //returns length and writes length deriv into the dlength argument.
        {
            double l = length();
            if (l == 0) {
                dlength = 1.0;
                return l;
            }
            dlength = (x * dx + y * dy) / l;
            return l;
        }

        // unlike other vectors in FreeCAD, this normalization creates a new vector instead of modifying existing one.
        DeriVector2 getNormalized() const
        {
            double dLen;
            double l = length(dLen);
            return (l == 0) ? DeriVector2(0, 0, dx, dy) : divD(l, dLen);
        }
        double scalarProd(const DeriVector2 &v2, double* dprd=nullptr) const;//calculates scalar product of two vectors and returns the result. The derivative of the result is written into argument dprd.
        // Returns the signed modulus of the cross product (this X v2). The derivative is stored in the dprd (if not-null provided)
        double crossProduct(const DeriVector2& v2, double* dprd = nullptr) const;

        DeriVector2 sum(const DeriVector2 &v2) const {//adds two vectors and returns result
            return DeriVector2(x + v2.x, y + v2.y,
                               dx + v2.dx, dy + v2.dy);}
        DeriVector2 subtr(const DeriVector2 &v2) const {//subtracts two vectors and returns result
            return DeriVector2(x - v2.x, y - v2.y,
                               dx - v2.dx, dy - v2.dy);}
        DeriVector2 mult(double val) const {
            return DeriVector2(x*val, y*val, dx*val, dy*val);}//multiplies the vector by a number. Derivatives are scaled.
        DeriVector2 multD(double val, double dval) const {//multiply vector by a variable with a derivative.
            return DeriVector2(x*val, y*val, dx*val+x*dval, dy*val+y*dval);}
        DeriVector2 divD(double val, double dval) const; //divide vector by a variable with a derivative
        DeriVector2 rotate90ccw() const {return DeriVector2(-y,x,-dy,dx);}
        DeriVector2 rotate90cw() const {return DeriVector2(y,-x,dy,-dx);}

        DeriVector2 rotateCCW_D(double angle, double dAngle)
        {
            double co = std::cos(angle), si = std::sin(angle);
            double dco = -si * dAngle, dsi = co * dAngle;
            double _x = x * co - y * si, _y = x * si + y * co;
            double _dx = (dx * co + dco * x) - (dy * si + y * dsi);
            double _dy = (dx * si + x * dsi) + (dy * si + y * dsi);
            return DeriVector2(_x, _y, _dx, _dy);
        }

        DeriVector2 rotateCW_D(double angle, double dAngle)
        {
            double co = std::cos(angle), si = std::sin(angle);
            double dco = -si * dAngle, dsi = co * dAngle;
            double _x = x * co + y * si, _y = - x * si + y * co;
            double _dx = (dx * co + dco * x) + (dy * si + y * dsi);
            double _dy = -(dx * si + x * dsi) + (dy * si + y * dsi);
            return DeriVector2(_x, _y, _dx, _dy);
        }

        // use lerp, one less mulpiplication per dimension
        //DeriVector2 linCombi(double m1, const DeriVector2 &v2, double m2) const {//linear combination of two vectors
        //    return DeriVector2(x*m1 + v2.x*m2, y*m1 + v2.y*m2,
        //                       dx*m1 + v2.dx*m2, dy*m1 + v2.dy*m2);}

        static DeriVector2 lerp(const DeriVector2& s, const DeriVector2& e, double t) {
            return DeriVector2(
                std::lerp(s.x, e.x, t), std::lerp(s.y, e.y, t), 
                std::lerp(s.dx, e.dx, t), std::lerp(s.dy, e.dy, t));
        }
    };

    ///////////////////////////////////////
    // Geometries
    ///////////////////////////////////////

    class Curve //a base class for all curve-based objects (line, circle/arc, ellipse/arc)
    {
    public:
        virtual ~Curve(){}
        //returns normal vector. The vector should point to the left when one
        // walks along the curve from start to end. Ellipses and circles are
        // assumed to be walked counterclockwise, so the vector should point
        // into the shape.
        //derivparam is a pointer to a curve parameter (or point coordinate) to
        // compute the derivative for. The derivative is returned through dx,dy
        // fields of DeriVector2.
        virtual DeriVector2 CalculateNormal(const Point &p, const double* derivparam = nullptr) const = 0;

        /**
        * @brief getRepresentativeSize : offers hints on how large the curve is; suggestion: implement
        *    it to use elements of the bounding box or other specific value derived from the 
        *    sizes defining the geometry of the curve.
        */
        virtual double getRepresentativeSize() const = 0;

        /**
         * @brief Value: returns point (vector) given the value of parameter
         * @param u: value of parameter
         * @param du: derivative of parameter by derivparam
         * @param derivparam: pointer to sketch parameter to calculate the derivative for
         * @return
         */
        virtual DeriVector2 Value(double u, double du, const double* derivparam = nullptr) const;

        //adds curve's parameters to pvec (used by constraints)
        virtual int PushOwnParams(VEC_pD &pvec) = 0;
        //recunstruct curve's parameters reading them from pvec starting from index cnt.
        //cnt will be incremented by the same value as returned by PushOwnParams()
        virtual void ReconstructOnNewPvec (VEC_pD &pvec, int &cnt) = 0;
        virtual Curve* Copy() = 0; //DeepSOIC: I haven't found a way to simply copy a curve object provided pointer to a curve object.
    };

    class Line: public Curve
    {
    public:
        Line(){}
        ~Line() override{}
        Point p1;
        Point p2;
        DeriVector2 CalculateNormal(const Point &p, const double* derivparam = nullptr) const override;
        DeriVector2 Value(double u, double du, const double* derivparam = nullptr) const override;
        int PushOwnParams(VEC_pD &pvec) override;
        void ReconstructOnNewPvec (VEC_pD &pvec, int &cnt) override;
        Line* Copy() override;

        double getRepresentativeSize() const override { return getRepresentativeSize(p1, p2); }

        static double getRepresentativeSize(double dx, double dy)
        {
            // taxicab/manhattan distance on extent
            return std::abs(dx) + std::abs(dy);
        }
        static double getRepresentativeSize(double sx, double sy, double ex, double ey)
        {
            // taxicab/manhattan distance between the ends
            return getRepresentativeSize(ex - sx, ey - sy);
        }
        static double getRepresentativeSize(const Point& p1, const Point& p2)
        {
            return getRepresentativeSize(p1.x ? *p1.x : 2, p1.y ? *p1.y : 2,
                                          p2.x ? *p2.x : 1, p2.y ? *p2.y : 1);
        }

    };

    class Circle: public Curve
    {
    public:
        Circle(){rad = nullptr;}
        ~Circle() override{}
        Point center;
        double *rad;
        DeriVector2 CalculateNormal(const Point &p, const double* derivparam = nullptr) const override;
        DeriVector2 Value(double u, double du, const double* derivparam = nullptr) const override;
        int PushOwnParams(VEC_pD &pvec) override;
        void ReconstructOnNewPvec (VEC_pD &pvec, int &cnt) override;
        Circle* Copy() override;

        virtual double getRepresentativeSize() const override { return rad ? std::abs(*rad) : 1; }
    };

    class Arc: public Circle
    {
    public:
        Arc()
        {
            startAngle = nullptr;
            endAngle = nullptr;
            rad = nullptr;
        }
        ~Arc() override {}
        double* startAngle;
        double* endAngle;
        //double *rad; //inherited
        Point start;
        Point end;
        //Point center; //inherited
        int PushOwnParams(VEC_pD& pvec) override;
        void ReconstructOnNewPvec(VEC_pD& pvec, int& cnt) override;
        Arc* Copy() override;

        virtual double getRepresentativeSize() const override
        {
            double r = rad ? std::abs(*rad) : 1;
            double sa = startAngle ? *startAngle : 0;
            double ea = endAngle ? *endAngle : 0;
            return std::abs((ea - sa) * r);
        }
    };

    class MajorRadiusConic: public Curve
    {
    public:
        ~MajorRadiusConic() override{}
        virtual double getRadMaj(const DeriVector2 &center, const DeriVector2 &f1, double b, double db, double &ret_dRadMaj) const = 0;
        virtual double getRadMaj(double* derivparam, double &ret_dRadMaj) const = 0;
        virtual double getRadMaj() const = 0;
        //DeriVector2 CalculateNormal(Point &p, double* derivparam = 0) = 0;
    };

    class Ellipse: public MajorRadiusConic
    {
    public:
        Ellipse(){ radmin = nullptr;}
        ~Ellipse() override{}
        Point center;
        Point focus1;
        double *radmin;
        double getRadMaj(const DeriVector2 &center, const DeriVector2 &f1, double b, double db, double &ret_dRadMaj) const override;
        double getRadMaj(double* derivparam, double &ret_dRadMaj) const override;
        double getRadMaj() const override;
        DeriVector2 CalculateNormal(const Point &p, const double* derivparam = nullptr) const override;
        DeriVector2 Value(double u, double du, const double* derivparam = nullptr) const override;
        int PushOwnParams(VEC_pD &pvec) override;
        void ReconstructOnNewPvec (VEC_pD &pvec, int &cnt) override;
        Ellipse* Copy() override;

        virtual double getRepresentativeSize() const override {
            // circle returns radius, ellipse returns avg radius
            return (radmin ? std::abs(*radmin) : 1) + std::abs(getRadMaj())/2;
        }
    };

    class ArcOfEllipse: public Ellipse
    {
    public:
        ArcOfEllipse(){startAngle=nullptr;endAngle=nullptr;radmin = nullptr;}
        ~ArcOfEllipse() override{}
        double *startAngle;
        double *endAngle;
        //double *radmin; //inherited
        Point start;
        Point end;
        //Point center;  //inherited
        //double *focus1.x; //inherited
        //double *focus1.y; //inherited
        int PushOwnParams(VEC_pD &pvec) override;
        void ReconstructOnNewPvec (VEC_pD &pvec, int &cnt) override;
        ArcOfEllipse* Copy() override;
    };

    class Hyperbola: public MajorRadiusConic
    {
    public:
        Hyperbola(){ radmin = nullptr;}
        ~Hyperbola() override{}
        Point center;
        Point focus1;
        double *radmin;
        double getRadMaj(const DeriVector2 &center, const DeriVector2 &f1, double b, double db, double &ret_dRadMaj) const override;
        double getRadMaj(double* derivparam, double &ret_dRadMaj) const override;
        double getRadMaj() const override;
        DeriVector2 CalculateNormal(const Point &p, const double* derivparam = nullptr) const override;
        DeriVector2 Value(double u, double du, const double* derivparam = nullptr) const override;
        int PushOwnParams(VEC_pD &pvec) override;
        void ReconstructOnNewPvec (VEC_pD &pvec, int &cnt) override;
        Hyperbola* Copy() override;
        virtual double getRepresentativeSize() const override
        {
            return (radmin ? std::abs(*radmin) : 1) + std::abs(getRadMaj())/2.0;
        }
    };

    class ArcOfHyperbola: public Hyperbola
    {
    public:
        ArcOfHyperbola(){startAngle=nullptr;endAngle=nullptr;radmin = nullptr;}
        ~ArcOfHyperbola() override{}
        // parameters
        double *startAngle;
        double *endAngle;
        Point start;
        Point end;
        // interface helpers
        int PushOwnParams(VEC_pD &pvec) override;
        void ReconstructOnNewPvec (VEC_pD &pvec, int &cnt) override;
        ArcOfHyperbola* Copy() override;
     };

    class Parabola: public Curve
    {
    public:
        Parabola(){ }
        ~Parabola() override{}
        Point vertex;
        Point focus1;
        DeriVector2 CalculateNormal(const Point &p, const double* derivparam = nullptr) const override;
        DeriVector2 Value(double u, double du, const double* derivparam = nullptr) const override;
        int PushOwnParams(VEC_pD &pvec) override;
        void ReconstructOnNewPvec (VEC_pD &pvec, int &cnt) override;
        Parabola* Copy() override;
        virtual double getRepresentativeSize() const override
        { // the taxicab semiperimeter of the latus rectum & directrix square (this is why the 4*)
            return 4*Line::getRepresentativeSize(vertex, focus1);
        }
    };

    class ArcOfParabola: public Parabola
    {
    public:
        ArcOfParabola(){startAngle=nullptr;endAngle=nullptr;}
        ~ArcOfParabola() override{}
        // parameters
        double *startAngle;
        double *endAngle;
        Point start;
        Point end;
        // interface helpers
        int PushOwnParams(VEC_pD &pvec) override;
        void ReconstructOnNewPvec (VEC_pD &pvec, int &cnt) override;
        ArcOfParabola* Copy() override;
    };

    class BSpline: public Curve
    {
    public:
        BSpline(){periodic=false;degree=2;}
        ~BSpline() override{}
        // parameters
        VEC_P poles; // TODO: use better data structures so poles.x and poles.y
        VEC_pD weights;
        VEC_pD knots;
        // dependent parameters (depends on previous parameters,
        // but an "arcrules" constraint alike would be required to gain the commodity of simple coincident
        // with endpoint constraints)
        Point start;
        Point end;
        // not solver parameters
        VEC_I mult;
        int degree;
        bool periodic;
        VEC_I knotpointGeoids; // geoids of knotpoints as to index Geom array
        VEC_D flattenedknots; // knot vector with repetitions for multiplicity and "padding" for periodic spline
        // interface helpers
        DeriVector2 CalculateNormal(const Point &p, const double* derivparam = nullptr) const override;
        DeriVector2 Value(double u, double du, const double* derivparam = nullptr) const override;
        int PushOwnParams(VEC_pD &pvec) override;
        void ReconstructOnNewPvec (VEC_pD &pvec, int &cnt) override;
        BSpline* Copy() override;
        /// finds the value B_i(x) such that spline(x) = sum(poles[i] * B_i(x))
        /// i is index of control point
        /// x is the point at which combination is needed
        /// k is the range in `flattenedknots` that contains x
        double getLinCombFactor(double x, size_t k, size_t i);
        void setupFlattenedKnots();

        // Curent implementation returns the taxicab distance of the spatial extent of the poles
        double getRepresentativeSize() const override;
    };

} //namespace GCS

#endif // PLANEGCS_GEO_H
