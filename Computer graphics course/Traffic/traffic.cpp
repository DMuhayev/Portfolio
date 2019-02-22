#include <iostream>
#include <GL/glut.h>
#include <vector>
#include <tuple>
#include <cmath>
#include <time.h>
#include <thread>
#include <chrono>


using namespace std;

//int size = 0;
int SoW = 800;
int iter = 0;
int count_of_idiots = 0;

class point
{
public:
	float x;
	float y;

	point (float i = 0, float j = 0)
	{
		x = i;
		y = j;
	}

};

point make_point (const float i, const float j)
{
	point a(i,j);
	return a;
}

point start(0.5, -1);
point centre(0, 0.4);

class car
{
public:
	float vel;
	point pos;
	point vec;
	bool is_stupid;


	car (float vel1 = 0, point pos1 = make_point(0, 0), point vec1 = make_point(0, 0), bool is = false)
	{
		vel = vel1;
		pos = pos1;
		vec = vec1;
		is_stupid = is;
	}
};

vector <car> traf;

//car first(0.02, start, make_point(0, 1));


car make_car(float v = 0.015, point p = start, point g = make_point(0, 1))
{
	car tmp(v, p, g);

	return tmp;
}
/*
	прямые: x = 0.5 и x = -0.5
	окружность: x^2 + (y - 0.4)^2 = 0.25
*/

void Reshape(int width, int height)
{
  glViewport(0, 0, width, height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(-1, 1, -1, 1);
  glMatrixMode(GL_MODELVIEW);
}

point normalize(point x)
{
	float alpha = sqrt(x.x*x.x + x.y*x.y);
	point tmp = make_point(x.x/alpha, x.y/alpha);
}

float distance(point x1, point x2)
{

	return sqrt( (x1.x-x2.x)*(x1.x-x2.x) + (x1.y-x2.y)*(x1.y-x2.y) );
}

bool is_place_free(void)
{
	bool o = true;

	for (int i = 0; i<traf.size(); i++)
	{
		if ( (traf[i].pos.x == 0.5) && (traf[i].pos.y < -0.9)) o = false;
	} 

	return o;
}

void GOGO(void)
{

	//while(true)
	//{
	srand(time(NULL));

	if (count_of_idiots == 0)
		{
			if ( rand()%5 < 1)
			{
				int random;
				int s = traf.size();

				for (int j1 = 0; j1 < 20; j1++)
				{
					random += rand()%s ;
				}

				random = (int) random/20;
		
				count_of_idiots = 1;
			
				traf[random].is_stupid = true;
				
				//}
			}
		}

	for (int i = 0; i < traf.size(); i++)
	{


		if (traf[i].pos.y < -1) 
		{
			if (traf.size() == 1) traf.clear(); 
			else 
			{
				vector <car> tmpv;

				for (int j = 0; j < traf.size(); j++)
				{
					if (j!=i) tmpv.push_back(traf[j]);
				}
				traf.clear();
				traf.swap(tmpv);
			}
		}

		if (!traf[i].is_stupid)
		{

			if ( (traf[i].vel < 0.02) && (distance(traf[i-1].pos, traf[i].pos) > 0.05)) traf[i].vel += 0.0005;

			if ( (i != 0) && (distance(traf[i-1].pos, traf[i].pos) < 0.1) && (traf[i-1].vel < traf[i].vel) )
			{
				if (distance(traf[i-1].pos, traf[i].pos) < 0.05) traf[i].vel = traf[i-1].vel; else
				if (traf[i].vel >= 0.004) traf[i].vel -= 0.004;	else traf[i].vel = 0;		
			} 
		
		}
		else
		{
			if (traf[i].vel >= 0.005) traf[i].vel -= 0.005;	
			else
			{
				traf[i].is_stupid = false;
				count_of_idiots = 0;
			}

		}


		//1st part of field
		if ( (traf[i].pos.y <= 0.4) && (traf[i].pos.x >= 0))
		{
			//cout<<"traf[i]\n";
			traf[i].vec = make_point(0, 1);
			traf[i].pos = make_point(traf[i].pos.x + traf[i].vel*traf[i].vec.x, traf[i].pos.y + traf[i].vel*traf[i].vec.y);
		}

		//2nd part
		if ( (traf[i].pos.y <= 0.4) && (traf[i].pos.x < 0))
		{
			//cout<<"sec\n";
			traf[i].vec = make_point(0, -1);
			traf[i].pos = make_point(traf[i].pos.x + traf[i].vel*traf[i].vec.x, traf[i].pos.y + traf[i].vel*traf[i].vec.y);
		}

		//3rd part
		if ((traf[i].pos.y > 0.4) && (traf[i].pos.x >= 0))
		{
			
			//cout<<"third\n";
			traf[i].vec = make_point(centre.y - traf[i].pos.y, traf[i].pos.x - centre.x);
			traf[i].vec = normalize(traf[i].vec);
			point tmp = make_point(traf[i].pos.x + traf[i].vel*traf[i].vec.x, traf[i].pos.y + traf[i].vel*traf[i].vec.y);
			float tmp1 = sqrt( (tmp.x - centre.x)*(tmp.x - centre.x) + (tmp.y - centre.y)*(tmp.y - centre.y) );
			float x1 = (0.5*(tmp.x - centre.x)/tmp1 ) + centre.x;
			float y1 = (0.5*(tmp.y - centre.y)/tmp1 ) + centre.y;
			traf[i].pos = make_point(x1, y1);
		}

		//4th part
		if ((traf[i].pos.y > 0.4) && (traf[i].pos.x < 0))
		{
			//cout<<"third\n";
			traf[i].vec = make_point(centre.y - traf[i].pos.y, traf[i].pos.x - centre.x);
			traf[i].vec = normalize(traf[i].vec);
			point tmp = make_point(traf[i].pos.x + traf[i].vel*traf[i].vec.x, traf[i].pos.y + traf[i].vel*traf[i].vec.y);
			float tmp1 = sqrt( (tmp.x - centre.x)*(tmp.x - centre.x) + (tmp.y - centre.y)*(tmp.y - centre.y) );
			float x1 = (0.5*(tmp.x - centre.x)/tmp1 ) + centre.x;
			float y1 = (0.5*(tmp.y - centre.y)/tmp1 ) + centre.y;
			traf[i].pos = make_point(x1, y1);

		}

	}
	//}


}

void Draw(void)
{
	glMatrixMode(GL_MODELVIEW);
  	glLoadIdentity();
  	glClear(GL_COLOR_BUFFER_BIT);
  	glColor3f(0, 0, 0);

  	float PI = 3.1415926535898;
	int circle_points = 100;

	glEnable(GL_LINE_SMOOTH);
	glLineWidth(50.0);


	glColor3f(0, 0, 0);
	glBegin(GL_LINES);
			glVertex2f((float) 	-1, (float) 0.5);
			glVertex2f((float)  0.4, (float) 0.5);
			glVertex2f((float) 	-1, (float) -0.5);
			glVertex2f((float)  0.4, (float) -0.5);
	glEnd();

	glBegin(GL_LINE_STRIP);
		for( int i = 0; i < circle_points; i++)
      {
            float angle = 2*PI*i/circle_points;
            if (angle >= 3*PI/2)
            {
            	glVertex2f( cos(angle)/2 + 0.4, sin(angle)/2 );
            }
      }
      for( int i = 0; i < circle_points; i++)
      {
            float angle = 2*PI*i/circle_points;

            if (angle <= PI/2)
            {
            	glVertex2f( cos(angle)/2 + 0.4, sin(angle)/2 );
            }
      }
	glEnd(); 

	glColor3f(1, 1, 1);
	glPointSize(10.0);

	glBegin(GL_POINTS);
		for (int i = 0; i < traf.size(); i++)
		{
			glColor3f(1, 1, 1);
			if (traf[i].is_stupid) glColor3f(1, 0, 0); 
			glVertex2f(traf[i].pos.y, traf[i].pos.x);
		}
	glEnd();


//helping shit
/*
	glDisable(GL_LINE_SMOOTH);
	glColor3f(0, 0, 1);
	glLineWidth(1.0);
	glBegin(GL_LINES);
			glVertex2f((float) 	0.4, (float) 1);
			glVertex2f((float)  0.4, (float) -1);
			glVertex2f((float) 	-1, (float) 0);
			glVertex2f((float)  0.4, (float) 0);
	glEnd();
*/
 	 glFlush();  
}

void Keyboard( unsigned char key, int x, int y)
{


	if (key=='g')
	{
		GOGO();

		if (iter == 0)	
		{
			if (is_place_free())
			{
				traf.push_back( make_car() );
				iter = 4;
			}
		}

		else iter--;

	}

	glutPostRedisplay();
}

void IdleFunc(){
	glutPostRedisplay();
}

int main(int argc, char *argv[])
{

	glutInit(&argc, argv);
	glutInitWindowSize(SoW, SoW);
	glutInitWindowPosition(100, 100);

	glutInitDisplayMode(GLUT_RGB);
	glutCreateWindow("Traffic");

	glutReshapeFunc(Reshape);

	glutDisplayFunc(Draw);
	glutKeyboardFunc(Keyboard);
	glutIdleFunc(IdleFunc);

	glClearColor(0.5, 0.5, 0.2, 0);

	glutMainLoop();
	return 0;
}