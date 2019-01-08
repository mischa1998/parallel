#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <fstream>
#include <mpi.h>

using namespace std;

int main(int argc, char **argv)
{
	int size, rank, mx, matr_size = 0;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int n = 0;
	int power = atoi(argv[1]);
	if (rank > 0)
	{
		MPI_Recv(&n, 1, MPI_INT, rank-1, 100, MPI_COMM_WORLD, &status);
	}
	ifstream file_1;
	ofstream file_2;
	double a;
	if (rank == 0) {
		file_1.open("matr.txt");
		file_2.open("matr1.txt");
		ofstream file_3;
		file_3.open("b.txt");
		file_1 >> matr_size;
		vector<double> in_b(matr_size);
		vector<vector<double> > in_matrix(matr_size);
		for(size_t i = 0; i < matr_size; i++) {
			for(size_t j = 0; j < matr_size; j++) {
				file_1 >> a;
				in_matrix[i].push_back(a); 
			}
		}
		for(size_t i = 0; i < matr_size; i++) {
			file_1 >> in_b[i];
		}
		vector<vector<double> > out_matr(matr_size); 
		for(size_t i = 0; i < matr_size; i++) {
			in_b[i] = in_b[i] / in_matrix[i][i];
			for(size_t j = 0; j < matr_size; j++) {
				a = -in_matrix[i][j] / in_matrix[i][i];
				out_matr[i].push_back(a);
			}
		}
		for(size_t i = 0; i < matr_size; i++) {
			for(size_t j = 0; j < matr_size; j++) {
				file_2 << out_matr[i][j] << " ";
			}
			file_2 << endl;
		}
		for(size_t i = 0; i < matr_size; i++) {
			file_3 << in_b[i] << " ";
		}
		file_1.close();
		file_2.close();
		file_3.close();
	}
	MPI_Bcast(&matr_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (rank < (size-1)) {
		MPI_Send(&n, 1, MPI_INT, rank + 1, 100, MPI_COMM_WORLD);
	}
	mx = matr_size / size;
	int old_mx = mx;
	if ( rank == size-1)
	{
		mx = matr_size - (size-1)*mx;
	}
	cout << "Hello! " << matr_size << " " << mx << endl;
	fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);

	vector<vector<double> >matr(mx);
	vector<double> b(mx);
	for(size_t i = 0; i < mx; i++) {
		matr[i].resize(matr_size);
	}
	if (rank > 0)
	{
		MPI_Recv(&n, 1, MPI_INT, rank-1, 100, MPI_COMM_WORLD, &status);
	}
	ifstream input, vec;
	input.open("matr1.txt", ios_base::in);
	vec.open("b.txt", ios_base::in);
	double *vec_old = new double [matr_size];

	for(size_t i = 0; i < matr_size; i++) {
		vec >> vec_old[i]; 
	}
	vec.close();
	vec.open("b.txt", ios_base::in);
	if (rank < size - 1) {
		for(size_t i = 0; i < rank*mx; i++) {
			for(size_t j = 0; j < matr_size; j++) {
				input >> a;
			}
		}
		for(size_t i = 0; i < mx; i++) {
			for(size_t j = 0; j < matr_size; j++) {
				input >> a;
				matr[i][j] = a;
			}
		}
		for(size_t i = 0; i < rank*mx; i++) {
			vec >> a;
		}
		for(size_t i = 0; i < mx; i++) {
			vec >> b[i];
		}
	} else {
		for(size_t i = 0; i < matr_size - mx; i++) {
			for(size_t j = 0; j < matr_size; j++) {
				input >> a;
			}
		}
		for(size_t i = 0; i < mx; i++) {
			for(size_t j = 0; j < matr_size; j++) {
				input >> a;
				matr[i][j] = a;
			}
		}
		for(size_t i = 0; i < matr_size - mx; i++) {
			vec >> a;
		}
		for(size_t i = 0; i < mx; i++) {
			vec >> b[i];
		}
	}
	input.close();
	vec.close();
	if (rank < (size-1)) {
		MPI_Send(&n, 1, MPI_INT, rank + 1, 100, MPI_COMM_WORLD);
	}
	
	double *sendbuf = new double [size];
	double *our_vec = new double [mx];
	double *our = new double [mx];
	double *vec_new = new double [matr_size];
	int *mas_count_s = new int [size]; // размер отправляемых сообщений
	int *mas_disp_s = new int [size]; // смещение отправляемых сообщений
	int *mas_count_r = new int [size]; // размер принимаемых сообщений
	int *mas_disp_r = new int [size]; // смещение принимаемых сообщений
	for(size_t i = 0; i < size; i++) {
		mas_count_s[i] = mx;
		mas_count_r[i] = mx;
		mas_disp_s[i] = 0;
	}
	mas_count_r[size - 1] = matr_size - (size-1)*old_mx;
	for(size_t i = 0; i < size - 1; i++) {
		mas_disp_r[i] = i * old_mx;
	}
	mas_disp_r[size - 1] = (size-1)*old_mx;

	for(size_t i = 0; i < mx; i++) {
		our_vec[i] = 0;
		our[i] = 0;
	}

	double norm = 1;
	double epsilon = pow(10., -power);
	int iter = 0;
	do {
		norm = 0;
		for(size_t i = 0; i < mx; i++) {
			for(size_t j = 0; j < matr_size; j++) {
				our_vec[i] += matr[i][j] * vec_old[j];
			}
			our_vec[i] += b[i];
		}
		MPI_Alltoallv(&our_vec[0], mas_count_s, mas_disp_s, MPI_DOUBLE, &vec_new[0], mas_count_r, mas_disp_r, MPI_DOUBLE, MPI_COMM_WORLD);
		norm = 0;
		for(size_t i = 0; i < matr_size; i++) {
			a = fabs(vec_new[i] - vec_old[i]);
			if (a > norm) {
				norm = a;
			} 
		}
		//MPI_Allgather(&norm, 1, MPI_DOUBLE, sendbuf, 1, MPI_DOUBLE, MPI_COMM_WORLD);
		for(size_t i = 0; i < matr_size; i++) {
			vec_old[i] = vec_new[i];
		}
		iter++;
	} while (norm > epsilon);

	if (rank > 0)
	{
		MPI_Recv(&n, 1, MPI_INT, rank-1, 100, MPI_COMM_WORLD, &status);
	}
	cout << "rank " << rank << endl;
	cout << "iterations " << iter << endl;
	for(size_t i = 0; i < matr_size; i++) {
		cout << vec_new[i] << " ";
	}
	cout << endl;
	if (rank < (size-1)) {
		MPI_Send(&n, 1, MPI_INT, rank + 1, 100, MPI_COMM_WORLD);
	}
	MPI_Finalize();
	delete vec_old;
	delete sendbuf;
	delete our_vec;
	delete our;
	delete vec_new;
	delete mas_count_s; // размер отправляемых сообщений
	delete mas_disp_s; // смещение отправляемых сообщений
	delete mas_count_r; // размер принимаемых сообщений
	delete mas_disp_r;
	return 0;
}
