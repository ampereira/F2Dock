#include <cuda.h>


/* dynArray_cuda implements an array whose size can change when adding elements to the end of the array, mimmicking the C++ STL vector */

template<typename T>

class dynArray_cuda {

 private:
  T *elements;
  int size;
	int capacity;
	int initCapacity;

 public:
  __device__ dynArray_cuda() 
	{}

  __device__ dynArray_cuda(int capacity_) 
	{
    capacity = capacity_;
		initCapacity = capacity_;
    size = 0;
		elements = (T*) malloc(sizeof(T)*capacity); 
  }

	__device__ T get(int i)
	{
		return elements[i];
	}

	__device__ T* getAll()
	{
		return elements;
	}


	__device__ void set(int i, T val)
	{
		if (i < size)
			elements[i] = val;
	}

	__device__ int getSize()
	{
		return size;
	}

	__device__ bool empty()
	{
		if (size == 0)
			return true;
		else
			return false;
	}

	__device__ void clear()
	{
		size = 0;
		free(elements);
		elements = (T*) malloc(sizeof(T) * initCapacity);
	}

	__device__ void push_end(T element_)
	{
		if ( size < capacity )
		{
			elements[size] = element_;
			size++;
		}
		else
		{
			T *new_elements = (T*) malloc(sizeof(T)*capacity);  
			for (int i = 0; i < size; i++)
				new_elements[i] = elements[i];
			free(elements);

			capacity *= 2;
			elements = (T*) malloc(sizeof(T)*capacity);  

			for (int i = 0; i < size; i++)
				elements[i] = new_elements[i];
			free(new_elements);

			elements[size] = element_;
			size++;
		}
	}

	__device__ void freeElements()
	{
		size = 0;
		free(elements);
	}

};

