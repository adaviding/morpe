#pragma once
#include <vector>

namespace Morpe
{ 
	template<class T> class Tensor
	{
		//	Constructs an empty tensor
		public:		Tensor();
		//	Constructs a 1-D tensor of the specified size.
		public:		Tensor(int i);
		//	Constructs a 2-D tensor of the specified size on each dimension (i=rows, j=cols).
		public:		Tensor(int i, int j);
		//	Constructs a 3-D tensor of the specified size on each dimension (i=pages, j=rows, k=cols).
		public:		Tensor(int i, int j, int k);
		//	Constructs a 4-D tensor of the specified size on each dimension (i=?, j=pages, k=rows, l=cols).
		public:		Tensor(int i, int j, int k, int l);
		//	Constructs a N-D tensor of the specified size on each dimension (for an arbitrary number of dimensions).
		public:		Tensor(std::vector<int> dim);
		//	Constructs a new tensor by duplicating an existing tensor.
		public:		Tensor(Tensor<T> &copyFrom);
		public:		~Tensor();

		//	The dimensionality of the Tensor.
		public:		std::vector<int> Dim;
		//	The numerical values contained inside the Tensor.
		public:		T* Values;
		//	This vector is used by methods IndToSub and SubToInd.
		public:		std::vector<int> Stride;
		//	The number of values of type T that can be contained inside the Tensor.
		public:		int Size;

		//	The tensor will assume the same shape and the same values as the argument.
		public:		void CopyFrom(Tensor<T> &copyFrom);

		//	Converts the index to a subscript for this tensor.
		public:		std::vector<int> IndToSub(int i);
		//	Converts the subscript to an index for this tensor.
		public:		int SubToInd(std::vector<int> sub);
		//	Converts the subscript to an index for this tensor.
		public:		int SubToInd(int i, int j);
		//	Converts the subscript to an index for this tensor.
		public:		int SubToInd(int i, int j, int k);
		//	Converts the subscript to an index for this tensor.
		public:		int SubToInd(int i, int j, int k, int l);
		
		//	Returns true if the subscripts are within the tensor's range.  False otherwise.
		public:		bool Check(std::vector<int> sub); 
		//	Returns true if the index is within the tensor's range.  False otherwise.
		public:		bool Check(int i);
		//	Returns true if the subscripts are within the tensor's range.  False otherwise.
		public:		bool Check(int i, int j);
		//	Returns true if the subscripts are within the tensor's range.  False otherwise.
		public:		bool Check(int i, int j, int k);
		//	Returns true if the subscripts are within the tensor's range.  False otherwise.
		public:		bool Check(int i, int j, int k, int l);
			
		//	Gets the value for this tensor at the given subscript.  This operation may be unsafe.
		public:		T Get(std::vector<int> sub);
		//	Gets the value for this tensor at the given index.  This operation may be unsafe.
		public:		T Get(int i);
		//	Gets the value for this tensor at the given subscript.  This operation may be unsafe.
		public:		T Get(int i, int j);
		//	Gets the value for this tensor at the given subscript.  This operation may be unsafe.
		public:		T Get(int i, int j, int k);
		//	Gets the value for this tensor at the given subscript.  This operation may be unsafe.
		public:		T Get(int i, int j, int k, int l);

		//	Sets the value for this tensor at the given subscript.  This operation may be unsafe.
		public:		void Set(T value, std::vector<int> sub);
		//	Sets the value for this tensor at the given index.  This operation may be unsafe.
		public:		void Set(T value, int i);
		//	Sets the value for this tensor at the given subscript.  This operation may be unsafe.
		public:		void Set(T value, int i, int j);
		//	Sets the value for this tensor at the given subscript.  This operation may be unsafe.
		public:		void Set(T value, int i, int j, int k);
		//	Sets the value for this tensor at the given subscript.  This operation may be unsafe.
		public:		void Set(T value, int i, int j, int k, int l);

		//	Prepares the tensor to store information by allocating all necessary memory.
		protected:	void allocate();
	};

	template<class T> Tensor<T>::~Tensor()
	{
		if(this->Values!=NULL)
		{
			free(this->Values);
			this->Values = NULL;
			this->Size = 0;
		}
	}
	template<class T> Tensor<T>::Tensor()
	{
		this->Dim = std::vector<int>();
		this->Stride = std::vector<int>();
		this->Values = NULL;
	}
	template<class T> Tensor<T>::Tensor(int i)
	{
		this->Dim = std::vector<int> {i};
		this->allocate();
	}
	template<class T> Tensor<T>::Tensor(int i, int j)
	{
		//	The following line is broken for some reason:
		//this->Dim = std::vector<int>() {i, j};
		this->Dim = std::vector<int>(2);
		this->Dim.push_back(i);
		this->Dim.push_back(j);
		this->allocate();
	}
	template<class T> Tensor<T>::Tensor(int i, int j, int k)
	{
		this->Dim = std::vector<int>() {i, j, k};
		this->allocate();
	}
	template<class T> Tensor<T>::Tensor(int i, int j, int k, int l)
	{
		this->Dim = std::vector<int>() {i, j, k, l};
		this->allocate();
	}
	template<class T> Tensor<T>::Tensor(std::vector<int> dim)
	{
		this->Dim = std::vector<int>(dim);
		this->allocate();
	}
	template<class T> Tensor<T>::Tensor(Tensor<T> &copyFrom)
	{
		this->CopyFrom(&copyFrom);
	}

	template<class T> void CopyFrom(Tensor<T> &copyFrom)
	{
		this->Dim = std::vector<int>(copyFrom.Dim);
		this->allocate();
		int szBuff = sizeof(T)*this->Size;
		memcpy_s(this->Values, szBuff, copyFrom.Values, szBuff);
	}

	template<class T> void Tensor<T>::allocate()
	{
		if(this->Values!=NULL)
		{
			free(this->Values);
			this->Values = NULL;
			this->Size = 0;
		}

		this->Stride = std::vector<int>(this->Dim.size());
		int iDim = this->Dim.size();
		this->Size=1;
		while(--iDim > 0)
		{
			this->Stride[iDim] = this->Size;
			this->Size *= this->Dim[iDim];
		}
		this->Values = (int*)malloc(sizeof(T)*this->Size);
	}

	template<class T> bool Tensor<T>::Check(std::vector<int> sub)
	{
		int sz = sub.size();
		if(this->Dim.size() < sz)
			return false;
		for(int i=0; i<sz; i++)
			if(sz[i]<0 || sz[i]>=this->Dim[i])
				return false;
		return true;
	}
	template<class T> bool Tensor<T>::Check(int i)
	{
		if(this->Dim.size()>=1 && i>=0 && i<this->Size)
			return true;
		return false;
	}
	template<class T> bool Tensor<T>::Check(int i, int j)
	{
		if(this->Dim.size()==2 && i>=0 && j>=0 && i<this->Dim[0] && j<this->Dim[1])
			return true;
		return false;
	}
	template<class T> bool Tensor<T>::Check(int i, int j, int k)
	{
		if(Dim.size()==3 && i>=0 && j>=0 && k>=0 && i<Dim[0] && j<Dim[1] && k<Dim[2])
			return true;
		return false;
	}
	template<class T> bool Tensor<T>::Check(int i, int j, int k, int l)
	{
		if(Dim.size()==4 && i>=0 && j>=0 && k>=0 && l>=0 && i<Dim[0] && j<Dim[1] && k<Dim[2] && l<Dim[3])
			return true;
		return false;
	}

	template<class T> std::vector<int> Tensor<T>::IndToSub(int i)
	{
		int nDim = this->Dim.size();
		std::vector<int> output = std::vector<int>(nDim);
		for(int iDim=0; iDim<nDim; iDim++)
		{
			int sub = i/Stride[iDim];
			output[iDim] = sub;
			i %= Stride[iDim];
		}
	}
	template<class T> int Tensor<T>::SubToInd(std::vector<int> sub)
	{
		int output = 0;
		int n = sub.size();
		for(int i=0; i<n; i++)
			output += Stride[i]*sub[i];
	}
	template<class T> int Tensor<T>::SubToInd(int i, int j)
	{
		return Stride[0]*i + Stride[1]*j;
	}
	template<class T> int Tensor<T>::SubToInd(int i, int j, int k)
	{
		return Stride[0]*i + Stride[1]*j + Stride[2]*k;
	}
	template<class T> int Tensor<T>::SubToInd(int i, int j, int k, int l)
	{
		return Stride[0]*i + Stride[1]*j + Stride[2]*k + Stride[3]*l;
	}

	template<class T> T Tensor<T>::Get(std::vector<int> sub)
	{
		return this->Values[this->SubToInd(sub)];
	}
	template<class T> T Tensor<T>::Get(int i)
	{
		return this->Values[i];
	}
	template<class T> T Tensor<T>::Get(int i, int j)
	{
		return this->Values[this->SubToInd(i,j)];
	}
	template<class T> T Tensor<T>::Get(int i, int j, int k)
	{
		return this->Values[this->SubToInd(i,j,k)];
	}
	template<class T> T Tensor<T>::Get(int i, int j, int k, int l)
	{
		return this->Values[this->SubToInd(i,j,k,l)];
	}

	template<class T> void Tensor<T>::Set(T value, std::vector<int> sub)
	{
		this->Values[this->SubToInd(sub)] = value;
	}
	template<class T> void Tensor<T>::Set(T value, int i)
	{
		this->Values[i] = value; 
	}
	template<class T> void Tensor<T>::Set(T value, int i, int j)
	{
		this->Values[this->SubToInd(i,j)] = value; 
	}
	template<class T> void Tensor<T>::Set(T value, int i, int j, int k)
	{
		this->Values[this->SubToInd(i,j,k)] = value; 
	}
	template<class T> void Tensor<T>::Set(T value, int i, int j, int k, int l)
	{
		this->Values[this->SubToInd(i,j,k,l)] = value; 
	}
}