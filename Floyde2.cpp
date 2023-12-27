#include <mpi.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <chrono>
using namespace std;

void printMatrix(vector<int> matr, int n) {
	ofstream outf;
	outf.open("out.txt");
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if(matr[i * n + j] == INT_MAX)
				outf << "INF" << " ";
			else
				outf << matr[i * n + j] << " ";
		}
		outf << endl;
	}
}


int main(int argc, char** argv) {
	int rank, size;
	MPI_Status st;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	if (rank == 0) {
		//int m = 0, n = 0;
		//vector<int> matrix;
		//ifstream in;
		//in.open("in.txt");
		//in >> n; // вершин
		//MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
		//matrix.resize(n * n);
		//in >> m; // рёбер
		//int p = n / size;
		//for (int i = 0; i < m; i++) {
		//	int from = 0, to = 0, w = 0;
		//	in >> from; from--;
		//	in >> to; to--;
		//	in >> w;
		//	matrix[from * n + to] = w;
		//}
		//in.close();
		//for (int i = 0; i < n; i++)
		//	for (int j = 0; j < n; j++)
		//		if (matrix[i * n + j] == 0)
		//			if (i != j)
		//				matrix[i * n + j] = INT_MAX;
		int n;
		cout << "Enter the number of vertices:" << endl;
		cin >> n;
		int p = n / size;
		//int* Matr = new int[n * n];
		vector<int> Matr(n * n, 0);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				Matr[i * n + j] = 0;
			}
		}
		for (int i = 0; i < n; i++) {
			if (i != n - 1) {
				Matr[i * n + i + 1] = 1;
				Matr[(i + 1) * n + i] = INT_MAX; // генерируем матрицу смежностей
			}
			if (i % 3 == 0 && i != 0) {
				Matr[i * n + i - 2] = 1;
				Matr[(i - 2) * n + i] = INT_MAX;
			}
		}
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i == j) continue;
				if (Matr[i * n + j] == 0) Matr[i * n + j] = INT_MAX;
			}
		}

		auto start = chrono::high_resolution_clock::now();
		for (int i = 1; i < size; i++)
			MPI_Send(&Matr[i * p * n], p * n, MPI_INT, i, i, MPI_COMM_WORLD);
		for (int k = 0; k < n; k++) {
			MPI_Bcast(&Matr[k * n], n, MPI_INT, (int)(k / p), MPI_COMM_WORLD);
			for (int i = 0; i < p; i++)
				for (int j = 0; j < n; j++)
					if (Matr[i * n + k] < INT_MAX && Matr[k * n + j] < INT_MAX)
					{
						if (Matr[i * n + k] + Matr[k * n + j] < Matr[i * n + j])
							Matr[i * n + j] = Matr[i * n + k] + Matr[k * n + j];
						//matrix[i * n + j] = min(matrix[i * n + j], matrix[i * n + k] + matrix[k * n + j]);
					}	
		}
		for (int i = 1; i < size; i++)
			MPI_Recv(&Matr[i * p * n], p * n, MPI_INT, i, (size + i) * p, MPI_COMM_WORLD, &st);
		double time = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start).count();
		cout << "Time: " << setprecision(15) << time / 1000 << " ms" << endl;
		printMatrix(Matr, n);
		//ofstream out;
		//out.open("out.txt");
		//for (int i = 0; i < n; i++) {
		//	out << setfill(' ') << setw(4) << Matr[i * n];
		//	for (int j = 1; j < n; j++)
		//		out << setfill(' ') << setw(5) << Matr[i * n + j];
		//	out << endl;
		//}
		//out.close();
		
	}
	else {
		int n = 0;
		vector<int> list;
		vector<int> veck;
		MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
		int p = n / size;
		list.resize(p * n);
		veck.resize(n);
		MPI_Recv(&list[0], p * n, MPI_INT, 0, rank, MPI_COMM_WORLD, &st);
		for (int k = 0; k < n; k++) {
			if ((int)(k / p) == rank) {
				MPI_Bcast(&list[k % p * n], n, MPI_INT, (int)(k / p), MPI_COMM_WORLD);
				for (int i = 0; i < n; i++)
					veck[i] = list[k % p * n + i];
			}
			else {
				MPI_Bcast(&veck[0], n, MPI_INT, (int)(k / p), MPI_COMM_WORLD);
			}
			for (int i = 0; i < p; i++)
				for (int j = 0; j < n; j++)
					list[i * n + j] = min(list[i * n + j], list[i * n + k] + veck[j]);
		}
		MPI_Send(&list[0], p * n, MPI_INT, 0, (size + rank) * p, MPI_COMM_WORLD);
	}
	MPI_Finalize();
	return 0;
}
