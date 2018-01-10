#ifndef VECTOR_DATATYPE_H
#define VECTOR_DATATYPE_H

#include <malloc.h>


typedef struct Vector {

    /* how many elements? */
    unsigned int size;

    /* the actual vector */
    float *v;

} Vector;

/* vector functions */
/* build a new vector with a given size */
Vector buildVector(unsigned int s);

/* build a new vector with a given size and value */
Vector buildVectorWithValue(unsigned int s, float value);

/* delete a given vector */
void deleteVector(Vector vector);

/* copy a vector */
Vector copyVector(Vector vector);

/* Copy the vector A to the vector B */
void copyVectoAToB(Vector *a, Vector *b);

/* inner product */
float innerProduct(Vector a, Vector b);

/* cross product */
Vector crossProduct(Vector a, Vector b);

/* show a vector */
void showVector(Vector vector);

/* fill the entire vector */
void fillVector(Vector vector);

#endif
