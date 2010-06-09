#include "drawGraph.h"
#include <GL/glut.h>
#include <GL/gl.h>

/*
para desenhar o grafo eu uso a classe Ad
eu abro duas janelas uma rodando o AD* e outra o A*
para desenhar os v�rtices eu uso md->drawx[i] e md->drawy[i] ( que s�o as coordenadas do arquivo .co)
para desenhar as arestas eu uso uma array de pairs md->edge_array[i]
para pintar eu uso um array de cores md->colorastar[i] e md->color[i]
*/

Ad* md;

float sx, sy;
double minx, miny, maxx, maxy;
int color[12][3];
int arr[300000];

void drawFuncAstar ()
{
	int topo = 0;

glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	
	glLoadIdentity();
	glTranslatef(sx-1.0,sy-1.0,0);
	//glScalef( sx, sy, 1.0 );


	glPointSize(1.0);

	glColor3ub(0,0,0);
	//DRAW EDGES
	glBegin(GL_LINES);
	for(int i = 0; i < md->num_edges;i++){
		int u = md->edge_array[i].first, v = md->edge_array[i].second;
		double xu = (-minx + md->drawx[u])* 100.0 / (-minx+maxx),
			   yu = (-miny + md->drawy[u])* 100.0 / (-miny+maxy),
			   xv = (-minx + md->drawx[v])* 100.0 / (-minx+maxx),
			   yv = (-miny + md->drawy[v])* 100.0 / (-miny+maxy);
		glVertex3f(xu,yu,0);
		glVertex3f(xv,yv,0);
	}
	glEnd();

	// DRAW NODES
	glBegin(GL_POINTS);
	for(int i = 1; i <= md->num_nodes;i++){
		double x = (-minx + md->drawx[i])* 100.0 / (-minx+maxx),
			   y = (-miny + md->drawy[i])* 100.0 / (-miny+maxy);
		if(md->colorastar[i] == Ad::BLACK)
			arr[topo++] = i;
		glColor3ub(color[md->colorastar[i]][0],color[md->colorastar[i]][1],color[md->colorastar[i]][2]);
		glVertex3f(x,y,0);
	}
	glEnd();
	
	glPointSize(6.0);
	if(topo > 0){
		glBegin(GL_POINTS);
		for(int j = 0; j < topo; j++){
			//cout << "ok ";
			int i = arr[j];
			double x = (-minx + md->drawx[i])* 100.0 / (-minx+maxx),
				   y = (-miny + md->drawy[i])* 100.0 / (-miny+maxy);

			glColor3ub(color[md->colorastar[i]][0],color[md->colorastar[i]][1],color[md->colorastar[i]][2]);
			glVertex3f( x,y,0);
		}
		glEnd();
	}

    glutSwapBuffers();  
}
void drawFunc ()
{
	int topo = 0;

    glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	
	glLoadIdentity();
	glTranslatef(sx-1.0,sy-1.0,0);
	//glScalef( sx, sy, 1.0 );


	glPointSize(1.0);

	glColor3ub(120,120,120);
	//DRAW EDGES
	glBegin(GL_LINES);
	for(int i = 0; i < md->num_edges;i++){
		int u = md->edge_array[i].first, v = md->edge_array[i].second;
		double xu = (-minx + md->drawx[u])* 100.0 / (-minx+maxx),
			   yu = (-miny + md->drawy[u])* 100.0 / (-miny+maxy),
			   xv = (-minx + md->drawx[v])* 100.0 / (-minx+maxx),
			   yv = (-miny + md->drawy[v])* 100.0 / (-miny+maxy);
		glVertex3f(xu,yu,0);
		glVertex3f(xv,yv,0);
	}
	glEnd();

	// DRAW NODES
	glBegin(GL_POINTS);
	for(int i = 1; i <= md->num_nodes;i++){
		double x = (-minx + md->drawx[i])* 100.0 / (-minx+maxx),
			   y = (-miny + md->drawy[i])* 100.0 / (-miny+maxy);

		if(md->color[i] == Ad::BLACK)
			arr[topo++] = i;
		glColor3ub(color[md->color[i]][0],color[md->color[i]][1],color[md->color[i]][2]);
		glVertex3f( x,y,0);
	}
	glEnd();
	glPointSize(6.0);
	if(topo > 0){
		glBegin(GL_POINTS);
		for(int j = 0; j < topo; j++){
			//cout << "ok ";
			int i = arr[j];
			double x = (-minx + md->drawx[i])* 100.0 / (-minx+maxx),
				   y = (-miny + md->drawy[i])* 100.0 / (-miny+maxy);

			glColor3ub(color[md->color[i]][0],color[md->color[i]][1],color[md->color[i]][2]);
			glVertex3f( x,y,0);
		}
		glEnd();
	}

    glutSwapBuffers();  
}

void newBoundaries()
{
    
/*	sx = sy = 1.0;
	minx = miny =  99999999;
	maxx = maxy = -99999999;
	
	for(int i = 1; i <= md->num_nodes;i++) {
		//if( md->color[i] != Ad::GRAY){
		  minx = min(minx,md->drawx[i]);
		  miny = min(miny,md->drawy[i]);
		  maxx = max(maxx,md->drawx[i]);
		  maxy = max(maxy,md->drawy[i]);
		}
	}
	// Pega a maior diferenca e usa
	if((maxx - minx) > (maxy - miny))
	  miny = maxy - (maxx - minx);
	else
	  minx = maxx - (maxy - miny);
	  
	
	  
	  
	glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	//glOrtho(minx, maxx+6, miny, maxy+6, -20.0, 20.0);
	    //glOrtho( -85530767, 75530767, 40000000, 42085396, -20020.0, 20.0);
	glOrtho(-10.0, 110.0, -10.0, 110.0, -20.0, 20.0);
	//glViewport(minx,miny,maxx-minx,maxy-miny);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
*/}

void initFunc()
{
	color[Ad::GRAY][0] = 230;
	color[Ad::GRAY][1] = 230;
	color[Ad::GRAY][2] = 230;
	color[Ad::BLACK][0] = 0;
	color[Ad::BLACK][1] = 255;
	color[Ad::BLACK][2] = 0;
	color[Ad::RED][0] = 150;
	color[Ad::RED][1] = 150;
	color[Ad::RED][2] = 150;
	color[Ad::YELLOW][0] = 30;
	color[Ad::YELLOW][1] = 30;
	color[Ad::YELLOW][2] = 30;

	sx = sy = 1.0;
	minx = miny =  99999999;
	maxx = maxy = -99999999;
	
	for(int i = 1; i <= md->num_nodes;i++) {
	  minx = min(minx,md->drawx[i]);
	  miny = min(miny,md->drawy[i]);
	  maxx = max(maxx,md->drawx[i]);
	  maxy = max(maxy,md->drawy[i]);
	}
	
	glClearColor( 1.0f, 1.0f, 1.0f, 0.0f );
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity();
	//glOrtho(minx, maxx+6, miny, maxy+6, -20.0, 20.0);
	    //glOrtho( -85530767, 75530767, 40000000, 42085396, -20020.0, 20.0);
	glOrtho( -10.0, 110.0, -10.0, 110.0, -20.0, 20.0);
	glMatrixMode (GL_MODELVIEW);
	glLoadIdentity();
}
void initAstar(){
glClearColor( 1.0f, 1.0f, 1.0f, 0.0f );
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity();
	glOrtho( -10.0, 110.0, -10.0, 110.0, -20.0, 20.0);
    glMatrixMode (GL_MODELVIEW);
    glLoadIdentity();
}
void idleFunc() 
{
 /*   if ( rotacionaX ) {
         glRotatef( 0.7f, 1.0f, 0.0f, 0.0f );
         glutPostRedisplay();
    } 
    if ( rotacionaY) {
         glRotatef( 0.7f, 0.0f, 1.0f, 0.0f );
         glutPostRedisplay();
    }*/
	glutPostRedisplay();
}

void keyPressed(unsigned char key, int x, int y) {  
      //static int a = 0; static int b = 0; static int c = 10; static int d = 10;
	switch(key){
		case 'q': md->graphicTest();
			  newBoundaries();
			break;
		//case 'w': a+=10; break;
		//case 'a': b+=10; break;
		//case 's': c+=10; break;
		//case 'd': d+=10; break;
			
		default: break;
	}
	//glViewport(a,b,c,d);
}  

void drawGraph(Ad *ad){
	int argc = 1;
	char** argv = new char*[1]();
	char* ch =(char*) "ara";
	argv[0] = ch;

	md = ad;

	glutInit(&argc, argv);      
	glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);  

	//glutInitWindowSize(700,700);
	//glutInitWindowPosition(600,0);
	//glutCreateWindow("A*");
	//glutDisplayFunc(drawFuncAstar);     
	//initAstar();

	//janela
	glutInitWindowSize(700,700);
	glutInitWindowPosition(0,0);
	glutCreateWindow("Ad*");

	// local initialization
	initFunc();
	glutKeyboardFunc(keyPressed);
	// configure callbacks
	glutDisplayFunc(drawFunc);
	glutIdleFunc(idleFunc);
	//glutMouseFunc(mouseFunc);     

	glutMainLoop();
}
