#ifndef UNTITLED1_VECTOR_H
#define UNTITLED1_VECTOR_H

#include <istream>

struct Vector {
    Vector() = default;

    Vector(double x, double y);

    double Norm() const;

    friend Vector operator-(Vector first, Vector second);
    friend std::istream &operator>>(std::istream &in, Vector &vector);

private:
    double x_;
    double y_;
};

inline Vector operator-(Vector first, Vector second) {
    return {first.x_ - second.x_, first.y_ - second.y_};
}

inline std::istream &operator>>(std::istream &in, Vector &vector) {
    in >> vector.x_ >> vector.y_;
    return in;
}

#endif //UNTITLED1_VECTOR_H
