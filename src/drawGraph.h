#ifndef DGRAPH
#define DGRAPH

#include "Ad.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <GL/glut.h>
#include <GL/gl.h>

using namespace std;

void keyPressed(unsigned char key, int x, int y);
void drawFunc();
void initFunc();
void idleFunc();
void drawGraph(Ad *ad);

#endif