//
// Created by Jonathan Pereyra on 4/17/18.
//

#include "HyperbolicParaboloidFunc.h"

HyperbolicParaboloidFunc::HyperbolicParaboloidFunc(GLfloat a, GLfloat b){
    this->a = a;
    this->b = b;
}

GLfloat HyperbolicParaboloidFunc::function(GLfloat x, GLfloat y, GLfloat z){
    return (a*(x*x)+b*(y*y)-(z));
}

bool HyperbolicParaboloidFunc::isInside(GLfloat x, GLfloat y, GLfloat z){
    return function(x, y, z) >= 0;
}