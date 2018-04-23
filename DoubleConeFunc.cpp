//
// Created by Jonathan Pereyra on 4/17/18.
//

#include "DoubleConeFunc.h"

DoubleConeFunc::DoubleConeFunc(GLfloat a, GLfloat b){
    this->a = a;
    this->b = b;
}

GLfloat DoubleConeFunc::function(GLfloat x, GLfloat y, GLfloat z){
    return (a*(x*x)+b*(y*y)-(z*z));
}

bool DoubleConeFunc::isInside(GLfloat x, GLfloat y, GLfloat z){
    return function(x, y, z) <= 0;
}


