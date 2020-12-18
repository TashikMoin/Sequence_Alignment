//g++ SequenceAlignment.cpp -fopenmp -O3 -Wall -o SequenceAlignment
#include <string>
#include <omp.h>
#include <cstring>
#include <iostream>

using namespace std;

int Find_Minimum_Penalty(std::string x, std::string y, int Match_OR_Mismatch, int Gap_Penalty, int* First_Aligned_Gene, int* Second_Aligned_Gene);


int main()
{
	int Mismatch_Penalty;
	int Gap_Penalty;
	std::string First_Gene;
	std::string Second_Gene;
	std::cout << "Enter Mismatch Minimum_Penalty : ";
	std::cin >> Mismatch_Penalty;
	std::cout << "Enter Gap Minimum_Penalty : ";
	std::cin >> Gap_Penalty;
	std::cout << "Enter First Gene : ";
	std::cin >> First_Gene;
	std::cout << "Enter Second Gene : ";
	std::cin >> Second_Gene;
	std::cout << "Mismatch_Penalty = " << Mismatch_Penalty << std::endl;
	std::cout << "Gap_Penalty = " << Gap_Penalty << std::endl;

	int m = First_Gene.length(); 
	int n = Second_Gene.length(); 
	int l = m + n;
	int First_Aligned_Gene[l + 1], Second_Aligned_Gene[l + 1];

	int Minimum_Penalty = Find_Minimum_Penalty(First_Gene, Second_Gene, Mismatch_Penalty, Gap_Penalty, First_Aligned_Gene, Second_Aligned_Gene);
	int id = 1;
	int i;
	for (i = l; i >= 1; i--)
	{
		if ((char)Second_Aligned_Gene[i] == '_' && (char)First_Aligned_Gene[i] == '_')
		{
			id = i + 1;
			break;
		}
	}

	std::cout << "Minimum_Penalty in aligning the Genes = ";
	std::cout << Minimum_Penalty << std::endl;
	std::cout << "The Aligned Genes Are :" << std::endl;
	for (i = id; i <= l; i++)
	{
		std::cout << (char)First_Aligned_Gene[i];
	}
	std::cout << "\n";
	for (i = id; i <= l; i++)
	{
		std::cout << (char)Second_Aligned_Gene[i];
	}
	std::cout << "\n";

	return 0;
}


int Find_Minimum_Of_All_Three_Options(int a, int b, int c) 
{
	if (a <= b && a <= c) 
    {
		return a;
	}
	else if (b <= a && b <= c) 
    {
		return b;
	}
	else 
    {
		return c;
	}
}


int** Make_2D_Matrix(int Width, int Height)
{
	int** Matrix = new int* [Width];
	size_t Size = Width;
	Size *= Height;
	int* Matrix0 = new int[Size];
	if (!Matrix || !Matrix0)
	{
		std::cerr << "Find_Minimum_Penalty: new failed" << std::endl;
		exit(1);
	}
	Matrix[0] = Matrix0;
	for (int i = 1; i < Width; i++)
		Matrix[i] = Matrix[i - 1] + Height;

	return Matrix;
}

int Find_Minimum_Penalty(std::string x, std::string y, int Match_OR_Mismatch, int Gap_Penalty, int* First_Aligned_Gene, int* Second_Aligned_Gene)
{
	int i, j; 

	int m = x.length();
	int n = y.length(); 

	int** Matrix = Make_2D_Matrix(m + 1, n + 1);
	size_t Size = m + 1;
	Size *= n + 1;
	memset(Matrix[0], 0, Size);

	omp_set_num_threads(omp_get_max_threads());

	#pragma omp parallel for 
	for (i = 0; i <= m; i++)
	{
		Matrix[i][0] = i * Gap_Penalty;
	}

	#pragma omp parallel for 
	for (i = 0; i <= n; i++)
	{
		Matrix[0][i] = i * Gap_Penalty;
	}

	int End, Start;

	if (m == n)
	{
		Start = 1;
		End = n;
	}
	else if (m < n)
	{
		End = m; 
		Start = 0;
	}
	else
	{
		End = n;
		Start = m - n + 1;
	}


	#pragma omp parallel
	for (int k = 2; k <= End; k++)
	{
		#pragma omp for
		for (int j = k - 1; j > 0; j--)
		{
			
			int i = k - j;
			if (x[i - 1] == y[j - 1])
			{
				Matrix[i][j] = Matrix[i - 1][j - 1];
			}
			else
			{
				Matrix[i][j] = Find_Minimum_Of_All_Three_Options(Matrix[i - 1][j - 1] + Match_OR_Mismatch,
					Matrix[i - 1][j] + Gap_Penalty,
					Matrix[i][j - 1] + Gap_Penalty);
			}
		}
	}

	if (m < n)
	{
		#pragma omp parallel 
		for (int k = m + 1; k <= n; k++)
		{
			#pragma omp for 
			for (int i = 1; i <= m; i++)
			{
				int j = k - i;
				if (x[i - 1] == y[j - 1])
				{
					Matrix[i][j] = Matrix[i - 1][j - 1];
				}
				else
				{
					Matrix[i][j] = Find_Minimum_Of_All_Three_Options(Matrix[i - 1][j - 1] + Match_OR_Mismatch,
						Matrix[i - 1][j] + Gap_Penalty,
						Matrix[i][j - 1] + Gap_Penalty);
				}
			}
		}
	}
	else if (m > n)
	{
		#pragma omp parallel
		for (int k = n + 1; k <= m; k++)
		{
			#pragma omp for 
			for (int j = n; j > 0; j--)
			{

				int i = k - j;
				if (x[i - 1] == y[j - 1])
				{
					Matrix[i][j] = Matrix[i - 1][j - 1];
				}
				else
				{
					Matrix[i][j] = Find_Minimum_Of_All_Three_Options(Matrix[i - 1][j - 1] + Match_OR_Mismatch,
						Matrix[i - 1][j] + Gap_Penalty,
						Matrix[i][j - 1] + Gap_Penalty);
				}
			}
		}
	}

	if (Start == 0)
	{
		Start++;
	}

	#pragma omp parallel
	for (int k = Start; k <= m; k++)
	{
		#pragma omp for
		for (int i = k; i <= m; i++)
		{
			int j = k + n - i;
			if (i != 0 || j != 0)
			{
				if (x[i - 1] == y[j - 1])
				{
					Matrix[i][j] = Matrix[i - 1][j - 1];
				}
				else
				{
					Matrix[i][j] = Find_Minimum_Of_All_Three_Options(Matrix[i - 1][j - 1] + Match_OR_Mismatch,
						Matrix[i - 1][j] + Gap_Penalty,
						Matrix[i][j - 1] + Gap_Penalty);
				}
			}
		}

	}

	int l = n + m; 

	i = m; j = n;

	int Position_Of_First_Gene = l;
	int Position_Of_Second_Gene = l;

	while (!(i == 0 || j == 0))
	{
		if (x[i - 1] == y[j - 1])
		{
			First_Aligned_Gene[Position_Of_First_Gene--] = (int)x[i - 1];
			Second_Aligned_Gene[Position_Of_Second_Gene--] = (int)y[j - 1];
			i--; j--;
		}
		else if (Matrix[i - 1][j - 1] + Match_OR_Mismatch == Matrix[i][j])
		{
			First_Aligned_Gene[Position_Of_First_Gene--] = (int)x[i - 1];
			Second_Aligned_Gene[Position_Of_Second_Gene--] = (int)y[j - 1];
			i--; j--;
		}
		else if (Matrix[i - 1][j] + Gap_Penalty == Matrix[i][j])
		{
			First_Aligned_Gene[Position_Of_First_Gene--] = (int)x[i - 1];
			Second_Aligned_Gene[Position_Of_Second_Gene--] = (int)'_';
			i--;
		}
		else if (Matrix[i][j - 1] + Gap_Penalty == Matrix[i][j])
		{
			First_Aligned_Gene[Position_Of_First_Gene--] = (int)'_';
			Second_Aligned_Gene[Position_Of_Second_Gene--] = (int)y[j - 1];
			j--;
		}
	}
	while (Position_Of_First_Gene > 0)
	{
		if (i > 0) First_Aligned_Gene[Position_Of_First_Gene--] = (int)x[--i];
		else First_Aligned_Gene[Position_Of_First_Gene--] = (int)'_';
	}
	while (Position_Of_Second_Gene > 0)
	{
		if (j > 0) Second_Aligned_Gene[Position_Of_Second_Gene--] = (int)y[--j];
		else Second_Aligned_Gene[Position_Of_Second_Gene--] = (int)'_';
	}

	int ret = Matrix[m][n];

	delete[] Matrix[0];
	delete[] Matrix;

	return ret;
}

