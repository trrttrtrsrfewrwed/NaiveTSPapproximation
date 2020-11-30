#pragma once

#include <cmath>
#include <istream>


struct Vector {
    Vector(double x, double y): x_(x), y_(y) {}
    double x_;
    double y_;

    double Norm() {
        return std::sqrt(x_ * x_ + y_ * y_);
    }
};

Vector operator -(Vector first, Vector second) {
    return {first.x_ - second.x_, first.y_ - second.y_};
}

std::istream& operator >>(std::istream &in, Vector &vector){
    in >> vector.x_ >> vector.y_;
    return in;
}