#pragma once

namespace Morpe { namespace Static
{
	template <typename T> void Quicksort(T* x, int first, int last)
	{
		if(first<last)
		{
			T pivotvalue = x[first];
			int left = leftarg - 1;
			int right = rightarg + 1;

			for(;;)
			{
				while (x[--right] > pivotvalue);
				while (x[++left] < pivotvalue);

				if (left >= right) break;

				T temp = x[right];
				x[right] = x[left];
				x[left] = temp;
			}

			int pivot = right;
			Quicksort<T>(x, leftarg, pivot);
			Quicksort<T>(x, pivot + 1, rightarg);
		}
	};
	template <typename T> void QuicksortIndex(int* ind, T* x, int first, int last)
	{
		if(first<last)
		{
			T pivotvalue = x[ind[first]];
			int left = leftarg - 1;
			int right = rightarg + 1;

			for(;;)
			{
				while (x[ind[--right]] > pivotvalue);
				while (x[ind[++left]] < pivotvalue);

				if (left >= right) break;

				int temp = ind[right];
				ind[right] = ind[left];
				ind[left] = temp;
			}

			int pivot = right;
			QuicksortIndex<T>(ind, x, leftarg, pivot);
			QuicksortIndex<T>(ind, x, pivot + 1, rightarg);
		}
	};
}}