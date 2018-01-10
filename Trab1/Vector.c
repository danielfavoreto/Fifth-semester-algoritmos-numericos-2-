#include "Vector.h"

/* vector functions */
/* build a new vector with a given size */
Vector buildVector(unsigned int s) {

    /* build a vector */
    Vector vector;

    /* assign the size value */
    vector.size = s;

    /* alloc the vector */
    vector.v = (float*) malloc(s*sizeof(float));
    if (NULL == vector.v) {

        printf("\nMemmory Error! Could not allocate the vector\n");

    }

    /* return the vector */
    return(vector);

}

/* build a new vector with a given size */
Vector buildVectorWithValue(unsigned int s, float value) {

    /* build a vector */
    Vector vector;

    /* helper */
    unsigned int i;

    /* assign the size value */
    vector.size = s;

    /* alloc the vector */
    vector.v = (float*) malloc(s*sizeof(float));
    if (NULL == vector.v) {

        printf("\nMemmory Error! Could not allocate the vector\n");

    }

    for (i = 0; i < s; i++) {

        vector.v[i] = value;

    }

    /* return the vector */
    return(vector);

}

/* delete a given vector */
void deleteVector(Vector vector) {

    /* remove the vector */
    free(vector.v);

    return;

}

/* copy a vector */
Vector copyVector(Vector vector) {

    /* helpers */
    unsigned int i;
    Vector newVector;

    /* build a new vector with the same size */
    newVector = buildVector(vector.size);

    /* copy each element */
    for (i = 0; i < vector.size; i++) {

        /* copy the element */
        newVector.v[i] = vector.v[i];

    }

    /* return the new vector */
    return(newVector);

}

/* Copy the vector A to the vector B */
void copyVectoAToB(Vector *a, Vector *b) {

    unsigned int i;

    for (i = 0; i < a->size; i++) {

        /* copy the current element */
        b->v[i] = a->v[i];

    }

    return;

}

/* inner product */
float innerProduct(Vector a, Vector b) {

    /* helpers */
    unsigned int i;
    float result = 0.0f;

    if (a.size != b.size) {

        printf("Error! The vectors does not have the same sizes!\n");

    } else {

        /* multiply each element */
        for (i = 0; i < a.size; i++) {

            /* update the result */
            result += a.v[i]*b.v[i];

        }

    }

    /* return the desired result */
    return(result);

}

/* cross product */
Vector crossProduct(Vector a, Vector b) {

    /* build the new vector */
    Vector newVector = buildVector(3);

    if (3 != a.size || 3 != b.size) {

        /* error! */
        printf("Error! The input vector are not 3D vectors!\n");

    } else {

        /* get the cross product */
        newVector.v[0] = a.v[1]*b.v[2] - a.v[2]*b.v[1];
        newVector.v[1] = a.v[2]*b.v[0] - a.v[0]*b.v[2];
        newVector.v[2] = a.v[0]*b.v[1] - a.v[1]*b.v[0];

    }

    /* return the new Vector */
    return newVector;

}

/* show a vector */
void showVector(Vector vector) {

    unsigned int i;

    printf("[ ");
    for (i = 0; i < vector.size; i++) {

        printf("%.3f, ", vector.v[i]);

    }
    printf("]\n");

}

/* fill the entire vector */
void fillVector(Vector vector) {

    /* helpers */
    unsigned int i;

    printf("\nFilling the vector with size: %d\n", vector.size);

    for (i = 0; i < vector.size; i++) {

        printf("\nvector[%d] = ", i + 1);
//        vector.v[i] = GetFloat();

    }

    printf("\n");

    return;

}
