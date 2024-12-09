#ifndef __LINE_H__
#define __LINE_H__
#include "Vec4.h"
#include "Color.h"

class Line
{
public:
    Vec4 v0, v1;
    Color c0, c1 ;

    Line(Vec4 v0, Vec4 v1, Color c0, Color c1);
    Line(const Line &other);
    void swapLine(Line&line);
    friend std::ostream &operator<<(std::ostream &os, const Line &c);
};

#endif