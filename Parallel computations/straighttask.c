//This is simple straight realization of task
//my variant of task is 2 - where border conditions is: first type for x and y, and periodic for z
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <stdlib.h> 

typedef struct {	//step point stucture as we need three steps for approximation
	double before;
	double now;
	double after;
} step_point;

const double PI  =3.141592653589793238463;

int n = 32; //size of area (matrix) n^3
int init_type = 3; //type of initial condition: 1 is for a0*exp(-b*(r-r0)^2), 2 is for point source
double a0 = 1; //some constants for first type of initial condition
double b = 0.1;
double r0 = 1;
double tau = 0.1; //time step
int num_time_steps = 20; //guess
double L = 10;
double h; //grid step
long int all_time, write_time;

double init(int i, int j, int k);
double init_C1(int i, int j, int k);
double init_C2(int i, int j, int k);
double init_C3(int i, int j, int k);
double realU(int i, int j, int k, int t);
double realU_C1(int i, int j, int k, int t);
double realU_C2(int i, int j, int k, int t);
double realU_C3(int i, int j, int k, int t);

// ------------------------------------------------------------
int main( int argc, char** argv )
{
	h = L/n;
	if (argc != 2) 
	{
		printf("wrong input!\n");
		return -1;
	}
	else sscanf(argv[1], "%d", &n);

	all_time = time(NULL);
	//////////////////////////////////prepare first steps of approximation/////////////////////////////////////
	step_point Grid[n][n][n];

	FILE *res_file;
	res_file = fopen("1StraightRes.txt", "w");
	fprintf(res_file, "|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n");
	fprintf(res_file, "n = %d\n", n);
	for (int i=0; i<n; i++)
	{
		for (int j=0; j<n; j++)
		{
			for (int k=0; k<n; k++)
			{
				bool is_border = (i == 0) || (j == 0) || (k == 0) || (i == (n-1)) || (j == (n-1)) || (k == (n-1));
			
				if (is_border)
				{
					if ( (i == 0) || (i == (n-1)) ) 
					{
						Grid[i][j][k].before = 0;
						Grid[i][j][k].now = 0;
					}
					else if ( (j == 0) || (j == (n-1)) ) 
					{
						Grid[i][j][k].before = 0;
						Grid[i][j][k].now = 0;
					}
					else if (k == (n-1) )
					{
						Grid[i][j][k].before = (Grid[i][j][1].before + Grid[i][j][k-1].before)/2;
						Grid[i][j][k].now = (Grid[i][j][1].now + Grid[i][j][k-1].now)/2;
						Grid[i][j][0].before = Grid[i][j][k].before;
						Grid[i][j][0].now = Grid[i][j][k].now;
					}
				}
				else
				{

					Grid[i][j][k].before = init(i, j, k);
				
					double Laplasian = 		(init(i-1, j, k) - 2*init(i, j, k) + init(i+1, j, k))/(h*h);
					Laplasian = Laplasian + (init(i, j-1, k) - 2*init(i, j, k) + init(i, j+1, k))/(h*h);
					Laplasian = Laplasian + (init(i, j, k-1) - 2*init(i, j, k) + init(i, j, k+1))/(h*h);

					Grid[i][j][k].now = Grid[i][j][k].before + tau*tau*Laplasian/2;

				}
			}
		}
	}

	//////////////////////////////////////MAIN APPROXIMATION////////////////////////////////////
	for (int step = 0; step < num_time_steps; step++)
	{

		for (int i=0; i<n; i++)
		{
			for (int j=0; j<n; j++)
			{
				for (int k=0; k<n; k++)
				{
					bool is_border = (i == 0) || (j == 0) || (k == 0) || (i == (n-1)) || (j == (n-1)) || (k == (n-1));
				
					if (is_border)
					{
						if ( (i == 0) || (i == (n-1)) ) 
						{
							Grid[i][j][k].after = 0;
						}
						else if ( (j == 0) || (j == (n-1)) )
						{
							Grid[i][j][k].after = 0;
						}
						else if (k == (n-1) )
						{
							Grid[i][j][k].after = (Grid[i][j][1].after + Grid[i][j][k-1].after)/2;
							Grid[i][j][0].after = Grid[i][j][k].after;
						}
					}
					else
					{
						double Laplasian = 		(Grid[i-1][j][k].now - 2*Grid[i][j][k].now + Grid[i+1][j][k].now);
						Laplasian = Laplasian + (Grid[i][j-1][k].now - 2*Grid[i][j][k].now + Grid[i][j+1][k].now);
						Laplasian = Laplasian + (Grid[i][j][k-1].now - 2*Grid[i][j][k].now + Grid[i][j][k+1].now);
						Laplasian = Laplasian/(h*h);
						Grid[i][j][k].after = 2*Grid[i][j][k].now - Grid[i][j][k].before + tau*tau*Laplasian;
					}
				}
			}
		}
		//write to file
		write_time = time(NULL);
		
		if (step == num_time_steps - 1)	
		{
			//fprintf(res_file, "step = %d\n", step);
			double max = 0;

			for (int i=0; i<n; i++)
			{
				for (int j=0; j<n; j++)
				{
					for (int k=0; k<n; k++)
					{
						double tmp = fabs( realU(i, j, k, step) - Grid[i][j][k].before );
						if (tmp>max) max = tmp;
					}
				}
			}	
			fprintf(res_file, "diff = %f", max);
		}			
		
		write_time = time(NULL) - write_time;
		all_time = all_time - write_time;

		//preparation for next step
		for (int i=0; i<n; i++)
		{
			for (int j=0; j<n; j++)
			{
				for (int k=0; k<n; k++)
				{
					Grid[i][j][k].before = Grid[i][j][k].now;
					Grid[i][j][k].now = Grid[i][j][k].after;
					//Grid[i][j][k].after = 0;
				}
			}
		}
	}
	
	all_time = time(NULL) - all_time;
	fseek(res_file, 0, SEEK_SET);
	fprintf(res_file, "alg work time = %ld\n", all_time);
	fclose(res_file);
	return 0;
}
///////////////////////////////////END OF MAIN///////////////////////////////////////////////

double init(int i, int j, int k)
{
	switch (init_type) {
	case 1:
		return init_C1(i, j, k);
	case 2:
		return init_C2(i, j, k);
	case 3:
		return init_C3(i, j, k);
	}
}

double init_C1(int i, int j, int k)
{
	double r = h*sqrt(i*i + j*j + k*k);

	return a0*exp(-b*(r-r0)*(r-r0));
}

double init_C2(int i, int j, int k)
{
	int  p = (int) n/2;//some values for second type of initial condition
	double volume = n*n*n*h*h*h;

	zreturn volume*exp( (-1)*( (i*h-p)*(i*h-p) + (j*h-p)*(j*h-p) + (k*h-p)*(k*h-p) )/(4*0.01) );
}

double init_C3(int i, int j, int k)
{
	return sin(2*PI*i/n)*sin(2*PI*j/n)*sin(2*PI*k/n);
}

double realU(int i, int j, int k, int t)
{
	switch (init_type) {
	case 1:
		return realU_C1(i, j, k, t);
	case 2:
		return realU_C2(i, j, k, t);
	case 3:
		return realU_C3(i, j, k, t);
	}
}

double realU_C1(int i, int j, int k, int t)
{
	double r = h*sqrt(i*i + j*j + k*k);

	return a0*exp(-b*(r-r0)*(r-r0) + t*tau*t*tau);
}

double realU_C2(int i, int j, int k, int t)
{
	int p = (int) n/2;//some values for second type of initial condition
	double volume = n*n*n*h*h*h;

	return volume*exp( (-1)*((i*h-p)*(i*h-p) + (j*h-p)*(j*h-p) + (k*h-p)*(k*h-p))/(4*t*tau) );
}

double realU_C3(int i, int j, int k, int t)
{
	return sin(2*PI*i/n)*sin(2*PI*j/n)*sin(2*PI*k/n)*cos(2*PI*sqrt(3/(n*n*h*h))*t*tau);
}