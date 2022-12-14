#include <iostream>
#include <clocale>
#include <cstdio>
#include <ctime>
#include <fstream>
#include "mpi.h"

using namespace std;

int main(int *argc, char **argv) {
	int n, m, v;   // matrix dimensions and vector length
	double startTime, endTime; // time of sart and end of multiplication
	int rank, numProc, block;
	int* a, * b, * c, * buffer, * ans;
	FILE* inM = fopen("inM.txt", "rt"); // input matrix file
	FILE* inV = fopen("inV.txt", "rt"); // input vector file
	FILE* out = fopen("out.txt", "wt"); // output file

	MPI_Init(argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);// get current process number
	MPI_Comm_size(MPI_COMM_WORLD, &numProc);// get ammount of processes

	n = 5000;
	m = 5000;
	v = 5000;

	// Memory allocation for matrix and vectors
	a = (int*)malloc(sizeof(int) * n * m); // mathrix
	b = (int*)malloc(sizeof(int) * v); // vector
	c = (int*)malloc(sizeof(int) * v * 12); // result vector
	block = m / numProc; // block size
	buffer = (int*)malloc(sizeof(int) * block * n);
	ans = (int*)malloc(sizeof(int) * block);

	if (n != v) {
		cout << " Error: vector length must be equal to the number of matrix columns!" << endl;
		return -1;
	}

	// printf("My rank: %d, processes in work: %d\n", rank, numProc);

	if (rank==0) {
		startTime = MPI_Wtime();

		// reading matrix
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				fscanf(inM, "%d", &a[i*n + j]);
			}
			char skip =' ';
			while (skip != '\n') {
				fscanf(inM, "%c", &skip);
			}
			// fgets(NULL, n, inM);
		}
		fclose(inM);

		/*
		for (int i = 0; i < n*m; i++) {
			printf("matrix element: %d\n", a[i]);
		}
		*/

		// reading vector
		for (int i = 0; i < v; i++) {
			fscanf(inV, "%d", &b[i]);
		}
		fclose(inV);

		// sending vector to every process
		for (int i = 1; i < numProc; i++) {
			MPI_Send(b, v, MPI_INT, i, 0, MPI_COMM_WORLD);
		}

		// sending blocks of matrix to every process
		for (int l = 1; l < numProc; l++) {
			MPI_Send(a + (l - 1) * block * n, block * n, MPI_INT, l, 1, MPI_COMM_WORLD);
		}
		
		// getting process calculation
		for (int k = 1; k < numProc; k++) {
			MPI_Recv(ans, block, MPI_INT, k, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			// forwarding result to vector c
			for (int i = 0; i < block; i++) {
				c[(k - 1) * block + i] = ans[i];
			}
		}

		// calculate remaining data
		for (int i = (numProc - 1) * block; i < n; i++) {
			for (int j = 0; j < n; j++) {
				int temp = 0;
				for (int k = 0; k < n; k++)
					temp += a[i * n + k] * b[k + j];
				c[i + j] = temp;
			}
		}

		out = fopen("out.txt", "w");
		for (int i = 0; i < n; i++) {
			fprintf(out, "%d ", c[i]);
		}
		fclose(out);

		// calculation time stats
		endTime = MPI_Wtime();

		std::ofstream log_out;
		log_out.open("log.txt", std::ios::app);
		log_out << "data size: " << v << ", num of processes:\t" << numProc << ", elapsed time:\t" << endTime - startTime << endl;
		log_out.close();

		// printf("rank: %d, time: %lfs\n", rank, endTime - startTime);

		// Memory cleaning
		free(a);
		free(b);
		free(c);
		free(buffer);
		free(ans);
	} else {
		// read vector
		MPI_Recv(b, v, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		
		// read parts of matrix
		MPI_Recv(buffer, block * n, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		// calculate and send result to main process
		for (int i = 0; i < block; i++) {
			for (int j = 0; j < n; j++) {
				int temp = 0;
				for (int k = 0; k < n; k++)
					temp += buffer[i * n + k] * b[k];
				ans[i] = temp;
			}
		}
		// Отправить результат расчета в основной процесс
		MPI_Send(ans, block, MPI_INT, 0, 3, MPI_COMM_WORLD);
	}

	// end of MPI work
	MPI_Finalize();

	return 0;
}

