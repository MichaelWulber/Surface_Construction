////
// Created by Jonathan Pereyra on 4/17/18.
//

#include "ImplicitFunc.h"

#ifndef HYPERBOLICPARABOLOID_H_
#define HYPERBOLICPARABOLOID_H_


class HyperbolicParaboloidFunc: public ImplicitFunc {
private:
    GLfloat a;
    GLfloat b;
public:
    HyperbolicParaboloidFunc(GLfloat a, GLfloat b);
    bool isInside(GLfloat x, GLfloat y, GLfloat z);
    GLfloat function(GLfloat x, GLfloat y, GLfloat z);
};

#endif /* HYPERBOLICPARABOLOID_H_ */
