//This is MPI realization of task
//my variant of task is 2 - where border conditions is: first type for x and y, and periodic for z
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <mpi.h>
#include <time.h>

typedef struct {	//step point stucture as we need three steps for approximation
	double before;
	double now;
	double after;
} step_point;

const double PI  = 3.141592653589793238463;

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
void FindBlockBreak(int *x, int *y, int *z);
void FindBlockCoord(int rank,int *x, int *y, int *z, int x_div, int y_div, int z_div);
int FindBlockRank(int x, int y, int z, int x_div, int y_div);

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

	MPI_Init(&argc, &argv);
	int Comm_size;
	int Comm_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &Comm_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &Comm_rank);

	FILE *res_file;

	if (Comm_rank == 0) 
	{
		all_time = time(NULL);
		res_file = fopen("1MpiRes.txt", "w");
		fprintf(res_file, "|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n");
		fprintf(res_file, "n = %d\n", n);
	}

	MPI_Status status;

	int x_div, y_div, z_div;

	x_div = Comm_size;
	y_div = 1;
	z_div = 1;

	FindBlockBreak(&x_div, &y_div, &z_div);

	int x_step = (int) n/x_div;//these ones represents num of axis separations
	int y_step = (int) n/y_div;
	int z_step = (int) n/z_div;

	int x, y, z; //this block coords
	FindBlockCoord(Comm_rank, &x, &y, &z, x_div, y_div, z_div);

	int x_size, y_size, z_size;

	if (x != x_div - 1) x_size = x_step;
	else x_size = x_step + n%x_div;
	if (y != y_div - 1) y_size = y_step;
	else y_size = y_step + n%y_div;
	if (z != z_div - 1) z_size = z_step;
	else z_size = z_step + n%z_div;

	//////////////////////////////////prepare first steps of approximation/////////////////////////////////////
	//printf("rank = %d: init started\n", Comm_rank);
	step_point Grid[x_size][y_size][z_size];

	for (int i=0; i<x_size; i++)
	{
		int glob_i = i + x*x_step; //global coord of point

		for (int j=0; j<y_size; j++)
		{
			int glob_j = j + y*y_step;

			for (int k=0; k<z_size; k++)
			{
				
				int glob_k = k + z*z_step;

				Grid[i][j][k].after = 0;

				if (((x == 0) && (i == 0)) || ((x == x_div-1) && (i == x_size-1)))
				{
					Grid[i][j][k].before = 0;
					Grid[i][j][k].now = 0;
				}
				else if (((y == 0) && (j == 0)) || ((y == y_div-1) && (j == y_size-1)))
				{
					Grid[i][j][k].before = 0;
					Grid[i][j][k].now = 0;
				}
				else if ((z_div-1 == 0) && (k == z_size - 1))
				{
					Grid[i][j][k].before = (Grid[i][j][1].before + Grid[i][j][k-1].before)/2;
					Grid[i][j][k].now = (Grid[i][j][1].now + Grid[i][j][k-1].now)/2;
					Grid[i][j][0].before = Grid[i][j][k].before;
					Grid[i][j][0].now = Grid[i][j][k].now;
				}
				else if (((k == z_size - 1) && (z == 0)) || ((k == z_size - 1) && (z == z_div - 1)))
				{

					if ((k == z_size - 1) && (z == 0))
					{
						//send [i][j][1]
						int partner_rank = FindBlockRank(x, y, z_div - 1, x_div, y_div);
						MPI_Send(&Grid[i][j][1].before, 1, MPI_DOUBLE, partner_rank, 3, MPI_COMM_WORLD);
						MPI_Send(&Grid[i][j][1].now, 1, MPI_DOUBLE, partner_rank, 4, MPI_COMM_WORLD);
						//evaluate
						if (init_type == 1)
						{
							Grid[i][j][k].before = init_C1(glob_i, glob_j, glob_k);
						
							double Laplasian = 		(init_C1(glob_i-1, glob_j, glob_k) - 2*init_C1(glob_i, glob_j, glob_k) + init_C1(glob_i+1, glob_j, glob_k))/(h*h);
							Laplasian = Laplasian + (init_C1(glob_i, glob_j-1, glob_k) - 2*init_C1(glob_i, glob_j, glob_k) + init_C1(glob_i, glob_j+1, glob_k))/(h*h);
							Laplasian = Laplasian + (init_C1(glob_i, glob_j, glob_k-1) - 2*init_C1(glob_i, glob_j, glob_k) + init_C1(glob_i, glob_j, glob_k+1))/(h*h);

							Grid[i][j][k].now = Grid[i][j][k].before + tau*tau*Laplasian/2;
						}
						else
						{
							Grid[i][j][k].before = init_C2(glob_i, glob_j, glob_k);

							double Laplasian = 		(init_C2(glob_i-1, glob_j, glob_k) - 2*init_C2(glob_i, glob_j, glob_k) + init_C2(glob_i+1, glob_j, glob_k))/(h*h);
							Laplasian = Laplasian + (init_C2(glob_i, glob_j-1, glob_k) - 2*init_C2(glob_i, glob_j, glob_k) + init_C2(glob_i, glob_j+1, glob_k))/(h*h);
							Laplasian = Laplasian + (init_C2(glob_i, glob_j, glob_k-1) - 2*init_C2(glob_i, glob_j, glob_k) + init_C2(glob_i, glob_j, glob_k+1))/(h*h);

							Grid[i][j][k].now = Grid[i][j][k].before + tau*tau*Laplasian/2;
						}
						//recive [i][j][0]
						MPI_Recv(&Grid[i][j][0].before, 1, MPI_DOUBLE, partner_rank, 5, MPI_COMM_WORLD, &status);
						MPI_Recv(&Grid[i][j][0].now, 1, MPI_DOUBLE, partner_rank, 6, MPI_COMM_WORLD, &status);
					}
					else if ((k == z_size - 1) && (z == z_div - 1))
					{	
						//printf("rank = %d: here\n", Comm_rank);
						int partner_rank = FindBlockRank(x, y, 0, x_div, y_div);
						//recive [i][j][glob1]
						double tmpprevbefore;
						double tmpprevnow;
						MPI_Recv(&tmpprevbefore, 1, MPI_DOUBLE, partner_rank, 3, MPI_COMM_WORLD, &status);
						MPI_Recv(&tmpprevnow, 1, MPI_DOUBLE, partner_rank, 4, MPI_COMM_WORLD, &status);
						//evaluate
						Grid[i][j][k].before = (tmpprevbefore + Grid[i][j][k-1].before)/2;
						Grid[i][j][k].now = (tmpprevnow + Grid[i][j][k-1].now)/2;
						//send [i][j][glob0]
						MPI_Send(&Grid[i][j][k].before, 1, MPI_DOUBLE, partner_rank, 5, MPI_COMM_WORLD);
						MPI_Send(&Grid[i][j][k].now, 1, MPI_DOUBLE, partner_rank, 6, MPI_COMM_WORLD);
					}
				}
				else
				{
					if (init_type == 1)
					{
						Grid[i][j][k].before = init_C1(glob_i, glob_j, glob_k);
					
						double Laplasian = 		(init_C1(glob_i-1, glob_j, glob_k) - 2*init_C1(glob_i, glob_j, glob_k) + init_C1(glob_i+1, glob_j, glob_k))/(h*h);
						Laplasian = Laplasian + (init_C1(glob_i, glob_j-1, glob_k) - 2*init_C1(glob_i, glob_j, glob_k) + init_C1(glob_i, glob_j+1, glob_k))/(h*h);
						Laplasian = Laplasian + (init_C1(glob_i, glob_j, glob_k-1) - 2*init_C1(glob_i, glob_j, glob_k) + init_C1(glob_i, glob_j, glob_k+1))/(h*h);

						Grid[i][j][k].now = Grid[i][j][k].before + tau*tau*Laplasian/2;
					}
					else
					{
						Grid[i][j][k].before = init_C2(glob_i, glob_j, glob_k);

						double Laplasian = 		(init_C2(glob_i-1, glob_j, glob_k) - 2*init_C2(glob_i, glob_j, glob_k) + init_C2(glob_i+1, glob_j, glob_k))/(h*h);
						Laplasian = Laplasian + (init_C2(glob_i, glob_j-1, glob_k) - 2*init_C2(glob_i, glob_j, glob_k) + init_C2(glob_i, glob_j+1, glob_k))/(h*h);
						Laplasian = Laplasian + (init_C2(glob_i, glob_j, glob_k-1) - 2*init_C2(glob_i, glob_j, glob_k) + init_C2(glob_i, glob_j, glob_k+1))/(h*h);

						Grid[i][j][k].now = Grid[i][j][k].before + tau*tau*Laplasian/2;
					}
				}

			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	//printf("rank = %d: init ended\n", Comm_rank);
	//////////////////////////////////////MAIN APPROXIMATION////////////////////////////////////

	for (int step = 0; step < num_time_steps; step++)
	{

		for (int i=0; i<x_size; i++)
		{
			for (int j=0; j<y_size; j++)
			{
				for (int k=0; k<z_size; k++)
				{
					//if (Grid[i][j][k].now != 0) printf("rank's = %d Grid[%d][%d][%d].now is %f\n", Comm_rank, i, j, k, Grid[i][j][k].now);
					bool is_border = (i == 0) || (j == 0) || (k == 0);
					is_border = is_border || (i == (x_size-1)) || (j == (y_size-1)) || (k == (z_size-1));

					if (is_border)
					{
						//sendrecv all!
						//записывать лапласиан, пересылать между крайними с конца элементами
						if (i == x_size - 1)
						{
							double tmp;
							if (x != x_div-1)
							{
								int up_partner_rank = FindBlockRank(x+1, y, z, x_div, y_div);
								MPI_Sendrecv(&Grid[i][j][k].now, 1, MPI_DOUBLE, up_partner_rank, 7, &tmp, 1, MPI_DOUBLE, up_partner_rank, 7, MPI_COMM_WORLD, &status);
								Grid[i][j][k].after = Grid[i][j][k].after + Grid[i-1][j][k].now - 2*Grid[i][j][k].now + tmp;
							}
							else
							{
								Grid[i][j][k].after = 0;
							}

							if (x != 0)
							{
								int down_partner_rank = FindBlockRank(x-1, y, z, x_div, y_div);
								MPI_Sendrecv(&Grid[0][j][k].now, 1, MPI_DOUBLE, down_partner_rank, 7, &tmp, 1, MPI_DOUBLE, down_partner_rank, 7, MPI_COMM_WORLD, &status);
								Grid[0][j][k].after = Grid[0][j][k].after + Grid[1][j][k].now - 2*Grid[0][j][k].now + tmp;
							}
							else
							{
								Grid[0][j][k].after = 0;
							}

							if ((j != 0) && (j != y_size-1))
							{
								Grid[i][j][k].after = Grid[i][j][k].after +  Grid[i][j-1][k].now - 2*Grid[i][j][k].now + Grid[i][j+1][k].now;
								Grid[0][j][k].after = Grid[0][j][k].after +  Grid[0][j-1][k].now - 2*Grid[0][j][k].now + Grid[0][j+1][k].now;
							}
							if ((k != 0) && (k != z_size-1))
							{
								Grid[i][j][k].after = Grid[i][j][k].after +  Grid[i][j][k-1].now - 2*Grid[i][j][k].now + Grid[i][j][k+1].now;
								Grid[0][j][k].after = Grid[0][j][k].after +  Grid[0][j][k-1].now - 2*Grid[0][j][k].now + Grid[0][j][k+1].now;
							}
						}
//printf("rank = %d: x ended\n", Comm_rank);
						if (j == y_size - 1)
						{
							double tmp;
							if (y != y_div-1)
							{
								int up_partner_rank = FindBlockRank(x, y+1, z, x_div, y_div);
								MPI_Sendrecv(&Grid[i][j][k].now, 1, MPI_DOUBLE, up_partner_rank, 7, &tmp, 1, MPI_DOUBLE, up_partner_rank, 7, MPI_COMM_WORLD, &status);
								Grid[i][j][k].after = Grid[i][j][k].after + Grid[i][j-1][k].now - 2*Grid[i][j][k].now + tmp;
							}
							else
							{
								Grid[i][j][k].after = 0;
							}

							if (y != 0)
							{
								int down_partner_rank = FindBlockRank(x, y-1, z, x_div, y_div);
								MPI_Sendrecv(&Grid[i][0][k].now, 1, MPI_DOUBLE, down_partner_rank, 7, &tmp, 1, MPI_DOUBLE, down_partner_rank, 7, MPI_COMM_WORLD, &status);
								Grid[i][0][k].after = Grid[i][0][k].after + Grid[i][1][k].now - 2*Grid[i][0][k].now + tmp;
							}
							else
							{
								Grid[i][0][k].after = 0;
							}

							if ((i != 0) && (i != x_size-1))
							{
								Grid[i][j][k].after = Grid[i][j][k].after + Grid[i-1][j][k].now - 2*Grid[i][j][k].now + Grid[i+1][j][k].now;
								Grid[i][0][k].after = Grid[i][0][k].after + Grid[i-1][0][k].now - 2*Grid[i][0][k].now + Grid[i+1][0][k].now;

								if ((k != 0) && (k != z_size-1))
								{
									Grid[i][j][k].after = Grid[i][j][k].after + Grid[i][j][k-1].now - 2*Grid[i][j][k].now + Grid[i][j][k+1].now;
									Grid[i][0][k].after = Grid[i][0][k].after + Grid[i][0][k-1].now - 2*Grid[i][0][k].now + Grid[i][0][k+1].now;
								}
							}
						}
//printf("rank = %d: y ended\n", Comm_rank);
						bool only_z = ((x == x_div-1) && (i == x_size-1)) || ((x == 0) && (i == 0)) ||
									  ((y == y_div-1) && (j == y_size-1)) || ((y == 0) && (j == 0));
						only_z = !only_z && (k == z_size-1);
						//if ((i == 0) && (j == 2) && (k == 31)) printf("rank = %d: only_z = %d\n", Comm_rank, only_z);
						if (only_z)
						{
							double tmp;

							if (z != 0)
							{
								int down_partner_rank = FindBlockRank(x, y, z-1, x_div, y_div);
								MPI_Sendrecv(&Grid[i][j][0].now, 1, MPI_DOUBLE, down_partner_rank, 7, &tmp, 1, MPI_DOUBLE, down_partner_rank, 7, MPI_COMM_WORLD, &status);
								Grid[i][j][0].after = Grid[i][j][0].after + Grid[i][j][1].now - 2*Grid[i][j][0].now + tmp;
							}

							if (z != z_div-1)
							{
								int up_partner_rank = FindBlockRank(x, y, z+1, x_div, y_div);
								MPI_Sendrecv(&Grid[i][j][k].now, 1, MPI_DOUBLE, up_partner_rank, 7, &tmp, 1, MPI_DOUBLE, up_partner_rank, 7, MPI_COMM_WORLD, &status);
								Grid[i][j][k].after = Grid[i][j][k].after + Grid[i][j][k-1].now - 2*Grid[i][j][k].now + tmp;
							}

							if ((i != 0) && (i != x_size-1) && (j != 0) && (j != y_size-1) && (z != 0) && (z != z_div-1))
							{
								Grid[i][j][k].after = Grid[i][j][k].after + Grid[i-1][j][k].now - 2*Grid[i][j][k].now + Grid[i+1][j][k].now;
								Grid[i][j][0].after = Grid[i][j][0].after + Grid[i-1][j][0].now - 2*Grid[i][j][0].now + Grid[i+1][j][0].now;
								Grid[i][j][k].after = Grid[i][j][k].after + Grid[i][j-1][k].now - 2*Grid[i][j][k].now + Grid[i][j+1][k].now;
								Grid[i][j][0].after = Grid[i][j][0].after + Grid[i][j-1][0].now - 2*Grid[i][j][0].now + Grid[i][j+1][0].now;
							}
						
						}

//printf("rank = %d: z ended\n", Comm_rank);
					}
					else
					{
						double Laplasian = 		(Grid[i-1][j][k].now - 2*Grid[i][j][k].now + Grid[i+1][j][k].now);
						Laplasian = Laplasian + (Grid[i][j-1][k].now - 2*Grid[i][j][k].now + Grid[i][j+1][k].now);
						Laplasian = Laplasian + (Grid[i][j][k-1].now - 2*Grid[i][j][k].now + Grid[i][j][k+1].now);

						Grid[i][j][k].after = Laplasian;//2*Grid[i][j][k].now - Grid[i][j][k].before + tau*tau*Laplasian;
					}

				}
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);

		for (int i=0; i<x_size; i++)
		{
			for (int j=0; j<y_size; j++)
			{
				for (int k=0; k<z_size; k++)
				{
					if ((k != 0) && (k != z_size)) Grid[i][j][k].after = 2*Grid[i][j][k].now - Grid[i][j][k].before + tau*tau*Grid[i][j][k].after/(h*h);
					
					bool only_z = ((x == x_div-1) && (i == x_size-1)) || ((x == 0) && (i == 0)) ||
								  ((y == y_div-1) && (j == y_size-1)) || ((y == 0) && (j == 0));
					only_z = !only_z && (k == z_size-1);
					//if ((i == 0) && (j == 2) && (k == 31)) printf("rank = %d: only_z = %d\n", Comm_rank, only_z);
					if (only_z)
					{

						if (z == 0)
						{
							if (z_div-1 != 0)
							{
								int partner_rank = FindBlockRank(x, y, z_div - 1, x_div, y_div);
								MPI_Send(&Grid[i][j][1].after, 1, MPI_DOUBLE, partner_rank, 3, MPI_COMM_WORLD);
								MPI_Recv(&Grid[i][j][0].after, 1, MPI_DOUBLE, partner_rank, 4, MPI_COMM_WORLD, &status);
							}
							else
							{
			//if ((i == 0) && (j == 2) && (k == 4)) printf("rank = %d: step = %d, Gij1 = %f, Gijk_1 = %f\n", Comm_rank, step, Gij1, Gijk_1);
								Grid[i][j][k].after = (Grid[i][j][1].after + Grid[i][j][k-1].after)/2;
								Grid[i][j][0].after = Grid[i][j][k].after;
							}
						}

						if ((z == z_div-1) && (z_div-1 != 0))
						{					
							double tmp;					

							int partner_rank = FindBlockRank(x, y, 0, x_div, y_div);
							MPI_Recv(&tmp, 1, MPI_DOUBLE, partner_rank, 3, MPI_COMM_WORLD, &status);
							Grid[i][j][k].after = (tmp + Grid[i][j][k-1].after)/2;
							MPI_Send(&Grid[i][j][k].after, 1, MPI_DOUBLE, partner_rank, 4, MPI_COMM_WORLD);

						}
					}

				}
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		//printf(" endstep\n");
		//write to file
		if (Comm_rank == 0) write_time = time(NULL);
		
		//if (step == num_time_steps - 1)
		{
			if (Comm_rank == 0)
			{
				double max = 0;
				//fprintf(res_file, "step = %d\n", step);
				for (int i=0; i<x_size; i++)
				{
					for (int j=0; j<y_size; j++)
					{
						for (int k=0; k<z_size; k++)
						{
							double tmp = fabs(realU(i, j, k, step) - Grid[i][j][k].before );
							if (tmp > max) max = tmp;
						}
					}
				}

				for (int rank = 1; rank < Comm_size; rank++)
				{
					//printf("rank = %d is writing to file\n", rank);
					int p_x, p_y, p_z; //partner's coords

					FindBlockCoord(rank, &p_x, &p_y, &p_z, x_div, y_div, z_div);

					int p_x_size, p_y_size, p_z_size; //partner's sizes of block 

					if (p_x != x_div - 1) p_x_size = x_step;
					else p_x_size = x_step + n%x_div;
					if (p_y != y_div - 1) p_y_size = y_step;
					else p_y_size = y_step + n%y_div;
					if (p_z != z_div - 1) p_z_size = z_step;
					else p_z_size = z_step + n%z_div;

					step_point tmp_Grid[p_x_size][p_y_size][p_z_size]; //partners part of grid

					int p_count = p_x_size*p_y_size*p_z_size*3;

					MPI_Recv(&tmp_Grid, p_count, MPI_DOUBLE, rank, 51, MPI_COMM_WORLD, &status);

					for (int i=0; i<p_x_size; i++)
					{
						for (int j=0; j<p_y_size; j++)
						{
							for (int k=0; k<p_z_size; k++)
							{	
								double tmp = fabs(realU(i, j, k, step) - tmp_Grid[i][j][k].before );
								if (tmp > max) max = tmp;
							}
						}
					}
				}
				fprintf(res_file, "step=%d: diff = %f\n", step, max);
				printf("step=%d: diff = %f\n", step, max);
			}
			else
			{
				int count = x_size*y_size*z_size*3;
				MPI_Send(&Grid, count, MPI_DOUBLE, 0, 51, MPI_COMM_WORLD);
			}

		}
		
		if (Comm_rank == 0)
		{
			write_time = time(NULL) - write_time;
			all_time = all_time - write_time;
			write_time = time(NULL) - all_time;
			printf("step=%d: alg work time = %ld\n", step, write_time);
		}

		//printf(" endwrite\n");		
		//preparation for next step
		for (int i=0; i<x_size; i++)
		{
			for (int j=0; j<y_size; j++)
			{
				for (int k=0; k<z_size; k++)
				{
					Grid[i][j][k].before = Grid[i][j][k].now;
					Grid[i][j][k].now = Grid[i][j][k].after;
					Grid[i][j][k].after = 0;
				}
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	//printf(" mainend\n");		
	if (Comm_rank == 0) 
	{
		all_time = time(NULL) - all_time;
		fseek(res_file, 0, SEEK_SET);
		fprintf(res_file, "alg work time = %ld\n", all_time);
		
		fclose(res_file);
	}
	MPI_Finalize();
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

	return volume*exp( (-1)*( (i*h-p)*(i*h-p) + (j*h-p)*(j*h-p) + (k*h-p)*(k*h-p) )/(4*0.01) );
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

void FindBlockBreak(int *x, int *y, int *z)
{

	int x_div = *x;
	int y_div = *y;
	int z_div = *z;


	int sum = 6666;
	int err = 2;

	while (err > 1)
	{
		int divider = 0;
		int tmp = sum;

		for (int i = 2; i<x_div; i++)
		{
			if (x_div % i == 0)
			{
				divider = i;
				break;
			}
		}

		if (divider == 0) break;

		x_div = (int) x_div / divider;
		z_div *= divider;

		if (z_div > y_div)
		{
			z_div += y_div;
			y_div = z_div - y_div;
			z_div = z_div - y_div;
		}

		if (y_div > x_div)
		{
			y_div += x_div;
			x_div = y_div - x_div;
			y_div = y_div - x_div;
		}

		if (z_div > y_div)
		{
			z_div += y_div;
			y_div = z_div - y_div;
			z_div = z_div - y_div;
		}

		sum = x_div + y_div + z_div;
		err = tmp - sum;

		*x = x_div;
		*y = y_div;
		*z = z_div;
	}

}

void FindBlockCoord(int rank,int *x, int *y, int *z, int x_div, int y_div, int z_div)
{
	*z = (int) rank/(x_div*y_div);
	*y = rank%(x_div*y_div);
	*x = *y % x_div;
	*y = (int) *y/x_div;
}

int FindBlockRank(int x, int y, int z, int x_div, int y_div)
{
	return z*x_div*y_div + y*x_div + x;
}