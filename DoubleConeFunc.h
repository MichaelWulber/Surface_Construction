//
// Created by Jonathan Pereyra on 4/17/18.
//

#include "ImplicitFunc.h"

#ifndef DOUBLECONE_H_
#define DOUBLECONE_H_


class DoubleConeFunc: public ImplicitFunc {
private:
    GLfloat a;
    GLfloat b;
public:
    DoubleConeFunc(GLfloat a, GLfloat b);
    bool isInside(GLfloat x, GLfloat y, GLfloat z);
    GLfloat function(GLfloat x, GLfloat y, GLfloat z);
};

#endif /* DOUBLECONE_H_ */
