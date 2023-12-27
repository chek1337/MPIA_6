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
			MPI_Send(&Matr[i * p * n], p * n, MPI_INT, i, i, MPI_COMM_WORLD); // Отправка блоков матрицы в другие процессы
		for (int k = 0; k < n; k++) { // Обход к-ых строчек
			int root = (int)(k / p);
			MPI_Bcast(&Matr[k * n], n, MPI_INT, root, MPI_COMM_WORLD); // если root == текущему ранку, то он отсылает всем остальным процессам k-ую строчку, иначе ждет от другого процесса к-ую строку
			for (int i = 0; i < p; i++) // Проход строчек в блоке
				for (int j = 0; j < n; j++) // Проход столбцов в строчке
					if (Matr[i * n + k] < INT_MAX && Matr[k * n + j] < INT_MAX) // Сам Флойд
							Matr[i * n + j] = min(Matr[i * n + k] + Matr[k * n + j], Matr[i*n+j]);
		}
		for (int i = 1; i < size; i++)
			MPI_Recv(&Matr[i * p * n], p * n, MPI_INT, i, (size + i) * p, MPI_COMM_WORLD, &st); // Принятие блоков от других процессов
		auto timerEnd = chrono::steady_clock::now();
		cout << chrono::duration_cast<chrono::milliseconds>(timerEnd - timerStart).count() << " ms";
		//printMatrix(Matr, n);
	}
	else 
	{
		int n=0;
		MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
		int strInRow = n / size;
		vector<int> block(strInRow*n, 0);
		vector<int> kRow(n, 0);

		MPI_Recv(&block[0], strInRow * n, MPI_INT, 0, rank, MPI_COMM_WORLD, &st); // Принятия блока каждым процессом
		for (int k = 0; k < n; k++) { // Обход к-ых строчек
			int root = (int)(k / strInRow);
			if (root == rank) { // к-ая строчка находится в блоке с тем же рангом или нет?
				// Да, находится
				MPI_Bcast(&block[k % strInRow * n], n, MPI_INT, root, MPI_COMM_WORLD); // Отдадим k-ую строчку остальным процессам
				for (int i = 0; i < n; i++)
					kRow[i] = block[k % strInRow * n + i]; // Скопируем k-ую строчку
			}
			else {
				// Нет, не находится
				MPI_Bcast(&kRow[0], n, MPI_INT, root, MPI_COMM_WORLD); // Принимаем от другого процесса
			}
			for (int i = 0; i < strInRow; i++) // Проход строчек в блоке
				for (int j = 0; j < n; j++) // Проход столбцов в строчке
					if (block[i * n + k] < INT_MAX && kRow[j] < INT_MAX)
							block[i * n + j] = min(block[i * n + k] + kRow[j], block[i * n + j]);
		}
		MPI_Send(&block[0], strInRow * n, MPI_INT, 0, (size + rank) * strInRow, MPI_COMM_WORLD); // Отправка блока в 1 процесс
	}
	MPI_Finalize();
	return 0;
}
