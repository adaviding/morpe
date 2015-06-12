//	Almon David Ing (This code was taken from a random location on the internet and then altered).
//	Ctr. Perceptual Systems, University of Texas at Austin
//	Written July 6, 2009

//	The classic QuickSort algorithm
void QuickSort(int* idx, double* x, int left, int right)
{
	int i = left, j = right;
	int tmp;
	double pivot = x[idx[(left + right) / 2]];

	//	Partition
	while (i <= j)
	{
		while (x[idx[i]] < pivot) i++;
		while (x[idx[j]] > pivot) j--;
		if (i <= j)
		{
			  tmp = idx[i];
			  idx[i++] = idx[j];
			  idx[j--] = tmp;
		}
	}

	//	Recursion
	if (left < j)
		QuickSort(idx, x, left, j);
	if (i < right)
		QuickSort(idx, x, i, right);
}