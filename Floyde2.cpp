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
			if (matr[i * n + j] == INT_MAX)
				outf << "INF" << "\t";
			else
				outf << matr[i * n + j] << "\t";
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
		int n;
		cout << "Enter the number of vertices:" << endl;
		cin >> n;
		MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
		int p = n / size;
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

		auto timerStart = chrono::steady_clock::now();
		for (int i = 1; i < size; i++)
			MPI_Send(&Matr[i * p * n], p * n, MPI_INT, i, i, MPI_COMM_WORLD);
		for (int k = 0; k < n; k++) {
			MPI_Bcast(&Matr[k * n], n, MPI_INT, (int)(k / p), MPI_COMM_WORLD);
			for (int i = 0; i < p; i++)
				for (int j = 0; j < n; j++)
					if (Matr[i * n + k] < INT_MAX && Matr[k * n + j] < INT_MAX)
							Matr[i * n + j] = min(Matr[i * n + k] + Matr[k * n + j], Matr[i*n+j]);
		}
		for (int i = 1; i < size; i++)
			MPI_Recv(&Matr[i * p * n], p * n, MPI_INT, i, (size + i) * p, MPI_COMM_WORLD, &st);
		auto timerEnd = chrono::steady_clock::now();
		cout << chrono::duration_cast<chrono::milliseconds>(timerEnd - timerStart).count() << " ms";
		//printMatrix(Matr, n);
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
					if (list[i * n + k] < INT_MAX && veck[j] < INT_MAX)
							list[i * n + j] = min(list[i * n + k] + veck[j], list[i * n + j]);
		}
		MPI_Send(&list[0], p * n, MPI_INT, 0, (size + rank) * p, MPI_COMM_WORLD);
	}
	MPI_Finalize();
	return 0;
}
