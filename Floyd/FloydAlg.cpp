#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <cmath>
#include "mpi.h"

using namespace std;

const int MAIN_PROC_RANK = 0;

const int matrix_size = 1000; //����������� �������

const int block_size = 200; // ����������� �����

const int BIG_VALUE = 10000; // ������� ����� - �������������

struct {
	int blockNum;
	int* block;
};

//�������� ������� ��������� ��� ����� � ������� ���������� ������ ���������, �.�. 1 � 2 � n, 2 � 1 � 3 � �.�.
void GetdAjacencyMatrixConnect(int**& matrix)
{
	matrix = new int* [matrix_size];
	for (int i = 0; i < matrix_size; i++)
	{
		matrix[i] = new int[matrix_size];
		for (int j = 0; j < matrix_size; j++)
		{
			if (i == j)
				matrix[i][j] = 0;
			else
				if (j == i + 1 || i == j + 1)
					matrix[i][j] = 1;
				else
					matrix[i][j] = BIG_VALUE;
		}
	}
}

void GetdAjacencyMatrixX(int**& matrix) // ������� ����������� �������
{
	matrix = new int* [matrix_size];
	int j;
	for (int i = 0; i < matrix_size; i++)
	{
		matrix[i] = new int[matrix_size];
		if (i < matrix_size / 2)
			if (i % 2 == 0)
			{
				j = matrix_size - i - 2;
			}
			else
			{
				j = matrix_size - i;
			}
		else if (i == matrix_size / 2)
		{
			j = i - 1;
		}
		else
			if (i % 2 == 0)
			{
				j = matrix_size - i;
			}
			else
			{
				j = matrix_size - i - 2;
			}

		matrix[i][j] = 1;
	}

	for (int i = 0; i < matrix_size; i++)
		for (int k = 0; k < matrix_size; k++)
		{
			if (i == k)
				matrix[i][k] = 0;
			else if (matrix[i][k] != 1)
				matrix[i][k] = BIG_VALUE;
		}
}

// ���� �� ������������� ����������� ������
void GetdAjacencyMatrixWheel(int**& matrix)
{
	matrix = new int* [matrix_size];
	for (int i = 0; i < matrix_size; i++)
	{
		matrix[i] = new int[matrix_size];
		for (int j = 0; j < matrix_size; j++)
		{
			if (i == j)
				matrix[i][j] = 0;
			else
				if (i == 0)
					matrix[i][j] = 1;
				else if (j == i + 1)
					matrix[i][j] = 1;
				else
					matrix[i][j] = BIG_VALUE;
		}
	}
}

// �������� ����� block1
// block2 � block3 - ���������������
void CalcBlock(int**& block1, int** block2, int** block3)
{
	int b;
	for (int k = 0; k < block_size; k++)
	{
		for (int i = 0; i < block_size; i++)
		{
			for (int j = 0; j < block_size; j++)
			{
				b = block2[i][k] + block3[k][j];
				if (block1[i][j] > b)
					block1[i][j] = b;

			}
		}
	}

}

void AlgFloyd(int****& blocks, int col_blocks)
{
	int myRank;
	int procCount;
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &procCount);

	int cBlocksCount = (col_blocks - 1) * 2; // ��������� ���������� ������ ���� � ��� ����, ����� ��������� �� �� ���������
	int uBlocksCount = (col_blocks - 1) * (col_blocks - 1); // � ������ ���� U

	int cBlocksLenght = cBlocksCount * block_size * block_size; // ������ ������ ��� ������ � ������ ��������
	int cBufLenght = procCount * cBlocksLenght; // ������ ������� ��� ��������� ���� ����������� ������ �
	int* crecvbuf;
	int* csendbuf;

	int uBlocksLenght = uBlocksCount * block_size * block_size; // ������ ������ ��� ������ U ����� ��������� 
	int uBufLenght = procCount * uBlocksLenght; // ������ ������� ��� ��������� ���� ����������� ������ �
	int* urecvbuf;
	int* usendbuf;


	for (int m = 0; m < col_blocks; m++)
	{
		// ��������� ���� D0
		CalcBlock(blocks[m][m], blocks[m][m], blocks[m][m]);


		csendbuf = new int[cBlocksLenght]; // ����� ��� ������ � ������ ��������

		for (int k = myRank, p = 0; k < cBlocksCount; k += procCount, p++) // ������� ����� ���� �
		{
			if (k < cBlocksCount / 2) // ���� ���� ���� �2
			{
				if (k < m) { // ���� ��������� ����� �� D0

					CalcBlock(blocks[m][k], blocks[m][m], blocks[m][k]);

					for (int t = 0; t < block_size; t++)
						for (int r = 0; r < block_size; r++)
							csendbuf[t * block_size + r + p * block_size * block_size] = blocks[m][k][t][r]; // ����� ���������� � �����
				}
				else { // ���� ��������� ������ �� D0

					CalcBlock(blocks[m][1 + k], blocks[m][m], blocks[m][1 + k]);

					for (int t = 0; t < block_size; t++)
						for (int r = 0; r < block_size; r++)
							csendbuf[t * block_size + r + p * block_size * block_size] = blocks[m][1 + k][t][r];
				}


			}
			else // ���� ���� ���� C1
			{
				int j = k % (cBlocksCount / 2);

				if (j < m) { // ���� ��������� ������ �� D0

					CalcBlock(blocks[j][m], blocks[j][m], blocks[m][m]);

					for (int t = 0; t < block_size; t++)
						for (int r = 0; r < block_size; r++)
							csendbuf[t * block_size + r + p * block_size * block_size] = blocks[j][m][t][r];
				}
				else { // ���� ��������� ����� �� D0

					CalcBlock(blocks[1 + j][m], blocks[1 + j][m], blocks[m][m]);

					for (int t = 0; t < block_size; t++)
						for (int r = 0; r < block_size; r++)
							csendbuf[t * block_size + r + p * block_size * block_size] = blocks[1 + j][m][t][r];
				}


			}

		}


		crecvbuf = new int[cBufLenght]; // ����� ��� ������ ������ � ��� ���� ���������

		MPI_Barrier(MPI_COMM_WORLD); // ����, ���� ��� ���������
		MPI_Allgather(csendbuf, cBlocksLenght, MPI_INT, crecvbuf, cBlocksLenght, MPI_INT, MPI_COMM_WORLD); // ���������� ������������ ������� �

		delete[] csendbuf;

		int* blockProc = new int[procCount] {0};
		for (int i = 0; i < cBlocksCount; i++) // ��������� ���������� ������������ ������ ������ ���������
		{
			int j = i % procCount;
			blockProc[j] += 1;
		}

		for (int k = 0; k < procCount; k++) {
			int disp = k * cBlocksLenght; //�������� �� �������-���������
			for (int p = 0; p < blockProc[k]; p++)
			{
				int cNum = p * procCount + k;
				if (cNum < cBlocksCount / 2) // ���� ���� ���� �2
				{
					if (cNum < m) { // ���� ���� ��������� ����� �� D0
						for (int t = 0; t < block_size; t++)
							for (int r = 0; r < block_size; r++) {
								blocks[m][cNum][t][r] = crecvbuf[t * block_size + r + p * block_size * block_size + disp];
							}
					}
					else { // ���� ���� ��������� ������ �� D0
						int ind = 1 + cNum;
						for (int t = 0; t < block_size; t++)
							for (int r = 0; r < block_size; r++) {
								blocks[m][ind][t][r] = crecvbuf[t * block_size + r + p * block_size * block_size + disp];
							}
					}
				}
				else // ���� ���� ���� C1
				{
					int j = cNum % (cBlocksCount / 2);

					if (j < m) { // ���� ��������� ������ �� D0

						for (int t = 0; t < block_size; t++)
							for (int r = 0; r < block_size; r++) {
								blocks[j][m][t][r] = crecvbuf[t * block_size + r + p * block_size * block_size + disp];
							}
					}
					else { // ���� ��������� ����� �� D0
						int ind = 1 + j;
						for (int t = 0; t < block_size; t++)
							for (int r = 0; r < block_size; r++) {
								blocks[ind][m][t][r] = crecvbuf[t * block_size + r + p * block_size * block_size + disp];
							}
					}


				}
			}

		}

		delete[] crecvbuf;

		usendbuf = new int[uBlocksLenght]; // ����� ��� ����������� ��������� ������ U

		for (int k = myRank, p = 0; k < uBlocksCount; k += procCount, p++) // c������ ����� U
		{
			int i = k / (col_blocks - 1);
			int j = k % (col_blocks - 1);

			i = (i < m) ? i : i + 1;
			j = (j < m) ? j : j + 1;

			CalcBlock(blocks[i][j], blocks[i][m], blocks[m][j]);

			for (int t = 0; t < block_size; t++)
				for (int r = 0; r < block_size; r++)
					usendbuf[t * block_size + r + p * block_size * block_size] = blocks[i][j][t][r];
		}

		urecvbuf = new int[uBufLenght];
		MPI_Barrier(MPI_COMM_WORLD); // ����, ���� ��� ���������
		MPI_Allgather(usendbuf, uBlocksLenght, MPI_INT, urecvbuf, uBlocksLenght, MPI_INT, MPI_COMM_WORLD); // ���������� ������������ ������� U


		delete[] blockProc;
		blockProc = new int[procCount] {0};
		for (int i = 0; i < uBlocksCount; i++)
		{
			int j = i % procCount;
			blockProc[j] += 1;
		}

		delete[] usendbuf;

		for (int k = 0; k < procCount; k++)
			for (int p = 0; p < blockProc[k]; p++)
			{
				int cNum = p * procCount + k;
				int i = cNum / (col_blocks - 1);
				int j = cNum % (col_blocks - 1);
				int disp = k * uBlocksLenght; //�������� �� �������-���������

				i = (i < m) ? i : i + 1;
				j = (j < m) ? j : j + 1;

				for (int t = 0; t < block_size; t++)
					for (int r = 0; r < block_size; r++) {
						int index = t * block_size + r + p * block_size * block_size + disp;
						blocks[i][j][t][r] = urecvbuf[index];
					}
			}

		delete[] urecvbuf;

		MPI_Barrier(MPI_COMM_WORLD); // ����, ���� ��� ���������

	}

}



int main(int argc, char** argv)
{
	int** adjMatr;
	//GetdAjacencyMatrixX(adjMatr);
	GetdAjacencyMatrixConnect(adjMatr); // ���������������� ������� D
	//GetdAjacencyMatrixWheel(adjMatr);

	int col_blocks = matrix_size / block_size;

	// ���������� ���������
	int procCount;

	// ���� �������� ��������
	int myRank;

	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &procCount);

	// ���� ��������� ������, ��� ������
	if (procCount > col_blocks* col_blocks) {
		if (myRank == MAIN_PROC_RANK) {
			cout << "Too many processes!" << endl;
		}

		MPI_Abort(MPI_COMM_WORLD, 0);
	}

	int**** blocks;
	blocks = new int*** [col_blocks];

	for (int p = 0; p < col_blocks; p++) // ��������� �����
	{
		blocks[p] = new int** [col_blocks];
		for (int k = 0; k < col_blocks; k++)
		{
			blocks[p][k] = new int* [block_size];
			for (int i = 0; i < block_size; i++)
			{
				blocks[p][k][i] = new int[block_size];
				for (int j = 0; j < block_size; j++)
				{
					blocks[p][k][i][j] = adjMatr[p * block_size + i][k * block_size + j];
				}
			}
		}
	}

	double startTime = MPI_Wtime();

	AlgFloyd(blocks, col_blocks);

	double endTime = MPI_Wtime();

	double elapsedTime = endTime - startTime;

	MPI_Finalize();

	if (myRank == MAIN_PROC_RANK) { // ������ ���������

		if (matrix_size < 101)
			for (int k = 0; k < col_blocks; k++) {
				for (int p = 0; p < block_size; p++)
				{
					for (int i = 0; i < col_blocks; i++) {
						for (int j = 0; j < block_size; j++)
						{
							cout << blocks[k][i][p][j] << " ";
						}

					}
					cout << endl;
				}
			}

		cout << endl << "Calculation time: " << elapsedTime << " s." << endl;

	}
}
