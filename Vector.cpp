#include <cmath>
#include "Vector.h"

double Vector::Norm() const {
    return std::sqrt(x_ * x_ + y_ * y_);
}

Vector::Vector(double x, double y) : x_(x), y_(y) {
}