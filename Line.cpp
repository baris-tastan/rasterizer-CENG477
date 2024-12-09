#include "Line.h"


Line::Line(Vec4 v0, Vec4 v1, Color c0, Color c1)
{
    this->v0 = v0;
    this->v1 = v1;
    this->c0 = c0;
    this->c1 = c1;
}

Line::Line(const Line &other)
{
    this->v0 = other.v0;
    this->v1 = other.v1;
    this->c0 = other.c0;
    this->c1 = other.c1;
}

void Line::swapLine(Line&line){
    Vec4 temp,temp2;
    Color tempc,tempc2;
    temp=this->v0;
    temp2=this->v1;
    tempc=this->c0;
    tempc2=this->c1;
    this->v0 =temp2;
    this->v1 = temp;
    this->c0 = tempc2;
    this->c1 = tempc;
}