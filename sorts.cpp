#include <iostream>
#include <stdint.h>
#include <cstdlib>
#include <time.h>
#include <cstring>
#include <vector>
#include <cmath>
#include <algorithm>

#define ARRAY_SIZE 131072
float CLOCKS_PER_MS = (float)CLOCKS_PER_SEC / 1000;

using namespace std;

class bitArr
{
 private:
   uint32_t * arr;
   uint32_t size;
   uint32_t bSize;
   uint32_t arrPos;
   uint32_t bitPos;
   uint32_t i;
 public:
    bitArr(const uint32_t &);
    ~bitArr();
    uint32_t getSize();
    void setPos(uint32_t &);
    void unSetPos(uint32_t &);
    void togglePos(uint32_t &);
    bool getValue(uint32_t &);
};

bitArr::bitArr(const uint32_t & newSize) 
{  
  bSize = newSize;
  size = ceil((float)bSize/32);

  arr = (uint32_t *)malloc(size*sizeof(uint32_t));
  
  for(i = 0; i < size; i++)
  {
    arr[i] = 0;  
  }

  if (arr == NULL) 
  {
   exit(1);
  }
}

bitArr::~bitArr() 
{ 
  if(arr) 
  { 
   free(arr); 
   arr = NULL;
  } 
} 

void bitArr::setPos(uint32_t & pos)
{
  arrPos = floor(((float)pos/32));
  bitPos = pos % 32;

  arr[arrPos] |= 1 << bitPos;

  return;
}

void bitArr::unSetPos(uint32_t & pos)
{
  arrPos = floor(((float)pos/32));
  bitPos = pos % 32;

  arr[arrPos] &= ~(1 << bitPos);

  return;
}

uint32_t bitArr::getSize()
{
 return bSize;
}

void bitArr::togglePos(uint32_t & pos)
{
  arrPos = floor(((float)pos/32));
  bitPos = pos % 32;

  arr[arrPos] ^= 1 << bitPos;

  return;
}

bool bitArr::getValue(uint32_t & pos)
{
  arrPos = floor(((float)pos/32));

  bitPos = pos % 32;

  return arr[arrPos] & (1 << bitPos);
}

template <class T>
void ptrSwap(T * a, T * b)
{
 unsigned int tmp = *a;
 *a = *b;
 *b = tmp;
 
 return;
}

template <class T>
void populateWithRandom(T * arr, const uint32_t & size, const uint32_t & upperBound, const uint32_t & lowerBound)
{
  uint32_t i = 0;
  
  srand ( time(NULL) );
  for(i = 0; i < size; i++)
  {
   arr[i] = rand() % (upperBound - lowerBound) + 1;
  }

  return;
}

template <class T>
void printArr(const T * arr, const uint32_t & size)
{
  uint32_t i = 0;

  for(i = 0; i < size; i++)
  {
    cout << arr[i] << " ";
  }
  cout << "\n";

  return;
}

template <class T>
void quickSort(T * arr, uint32_t left, uint32_t right) 
{
  uint32_t i = left;
  uint32_t j = right;
  uint32_t pivot = arr[(left + right) / 2];
 
  /* partition */
  while (i <= j) 
  {
    while (arr[i] < pivot)
    {
      i++;
    }
    while (arr[j] > pivot)
    {
      j--;
    }
    
    if (i <= j) 
    {
     ptrSwap(&arr[i],&arr[j]);
     i++;
     j--;
    }
   }
 
   /* recursion */
   if (left < j)
   {
    quickSort(arr, left, j);
   }
   if (i < right)
   {
    quickSort(arr, i, right);
   }
}

template <class T>
void bubbleSort(T * arr, const uint32_t & size) 
{
 bool swapped = true;
 uint32_t j = 0;
 uint32_t i = 0;
 
 while (swapped) 
 {
  swapped = false;
  j++;
  for (i = 0; i < size - j; i++) 
  {
   if (arr[i] > arr[i + 1]) 
   {
    ptrSwap(&arr[i],&arr[i+1]);
    swapped = true;
   }
  }
 }
 return;
}

template <class T>
void insertionSort(T * arr, const uint32_t & size) 
{
 uint32_t i = 0;
 uint32_t j = 0;
 T tmp;

 for (i = 1; i < size; i++) 
 {
   j = i;
   while (j > 0 && arr[j - 1] > arr[j]) 
   {
    ptrSwap(&arr[j-1],&arr[j]);
    j--;
   }
 }

 return;
}

template <class T>
void selectionSort(T * arr, uint32_t size) 
{
  uint32_t i;
  uint32_t j;
  uint32_t minIndex;
  T tmp; 
   
  for (i = 0; i < size - 1; i++) 
  {
    minIndex = i;
    for (j = i + 1; j < size; j++)
    {
     if (arr[j] < arr[minIndex])
     {
      minIndex = j;
     }
    }
    if (minIndex != i) 
    {
      ptrSwap(&arr[minIndex],&arr[i]);
    }   
  }
  return;
}

template <class T>
void merge(T * arr, uint32_t start, uint32_t end)
{
    uint32_t mid = floor((start + end) / 2);
    uint32_t i = 0;
    uint32_t j = start;
    uint32_t k = mid + 1;
    uint32_t l;

    T temp[end-start+1];

    while ( j <= mid && k <= end )
    {
     if ( arr[j] < arr[k] )
     {
      temp[i++] = arr[j++];
     }
     else
     {
      temp[i++] = arr[k++];
     }
    }

    while ( j <= mid )
    {
     temp[i++] = arr[j++];
    }

    while ( k <= end )
    {
     temp[i++] = arr[k++];
    }

    for ( l = start; l <= end; l++ )
    {
     arr[l] = temp[l-start];
    }
}

template <class T>
void mergeSort(T * arr, uint32_t start, uint32_t end)
{
  uint32_t mid;

  if ( start < end )
  {
    mid = floor((start + end) / 2);
    mergeSort(arr, start, mid);
    mergeSort(arr, mid + 1, end);
    merge(arr, start, end);
  }
}

template <class T>
void copyArr(const T *arr1, T *arr2, const uint32_t & size)
{
  uint32_t i = 0;

  for(i = 0; i < size; i++)
  {
    arr2[i] = arr1[i];
  }

  return;
}

template <class T>
void bucketSort(T * arr, uint32_t size, uint32_t max)
{
  T buckets [max];
  uint32_t i;
  uint32_t j;
  uint32_t k;

  for(j = 0; j < max; j++)
  {
    buckets[j] = 0;
  }
  for(i = 0; i < size; i++)
  {
    buckets[arr[i]]++;
  }
  for(i = 0, j = 0; j < max; j++)
  {
    for(k = buckets[j]; k > 0; k--)
    {
      arr[i++] = j;
    }
  }

  return;
}

template <class T>
uint32_t bucketDedupSort(T * arr, uint32_t size, uint32_t max)
{
  bitArr buckets(max);
  uint32_t i;
  uint32_t j;
  uint32_t k;

  for(i = 0; i < size; i++)
  {
    buckets.setPos(arr[i]);
  }
  for(i = 0, j = 0; j < max; j++)
  {
    if(buckets.getValue(j))
    {
     arr[i++] = j;
    }
  }
  
  return i;
}

template <class T>
void radixSort(T * arr,uint32_t size)
{
  T tempArray[size];
  uint32_t i;
  uint32_t j;
  uint32_t k;
  uint32_t r = 8;
  uint32_t R = 1 << r;
  uint32_t p = ((sizeof(uint32_t)*8) + r - 1U)/r;
  uint32_t count[R];

  for(i = 0; i < p; i++)
  {
    for(j = 0; j < R; j++)
    {
     count[j] = 0;
    }
    for(k = 0; k < size; k++)
    {
     count[(arr[k] >> (r * i)) & (R - 1U)]++;
     tempArray[k] = arr[k];
    }
    uint32_t pos = 0;
    for(j = 0; j < R; j++)
    {
     const uint32_t tmp = count[j];
     count[j] = pos;
     pos += tmp;
    }
    for(k = 0; k < size; k++)
    {
     j = (tempArray[k] >> (r*i)) & (R -1U);
     arr[count[j]++] = tempArray[k];
    }
  }

  return;
}

int compare (const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}

template <class T>
void shellSort(T *arr,uint32_t size)
{
     uint32_t i;
     uint32_t temp;
     uint32_t flag = 1;
     uint32_t d = size;

     while( flag || (d > 1))
     {
       flag = 0;
       d = (d+1) / 2;
       for (i = 0; i < (size - d); i++)
       {
        if (arr[i + d] > arr[i])
        {
          ptrSwap(&arr[i+d],&arr[i]);
          flag = 1;
        }
       }
     }

     return;
}

int main ()
{
 uint32_t toSort[ARRAY_SIZE]; 
 uint32_t sorted[ARRAY_SIZE];
 uint32_t min;
 uint32_t max;
 uint32_t newSize;

 unsigned long clicksStart;
 unsigned long clicksEnd;

 clicksStart = clock();
 populateWithRandom(toSort, ARRAY_SIZE,ARRAY_SIZE,0);
 clicksEnd = clock();
 
 cout << "Size of Array : " << ARRAY_SIZE << "\n"; 
 copyArr (toSort, sorted, ARRAY_SIZE) ;
 cout << "random population time taken: " << (float)(clicksEnd - clicksStart)/CLOCKS_PER_MS << "ms\n";

 clicksStart = clock();
 bubbleSort(sorted,ARRAY_SIZE-1);
 clicksEnd = clock();

 cout << "bubble sort time taken: " << (float)(clicksEnd - clicksStart)/CLOCKS_PER_MS << "ms\n";

 copyArr (toSort, sorted, ARRAY_SIZE);

 clicksStart = clock();
 selectionSort(sorted,ARRAY_SIZE -1);
 clicksEnd = clock();

 cout << "selection sort time taken: " << (float)(clicksEnd - clicksStart)/CLOCKS_PER_MS << "ms\n";

 copyArr (toSort, sorted, ARRAY_SIZE);

 clicksStart = clock();
 insertionSort(sorted,ARRAY_SIZE-1);
 clicksEnd = clock();

 cout << "insertion sort time taken: " << (float)(clicksEnd - clicksStart)/CLOCKS_PER_MS << "ms\n";

 copyArr (toSort, sorted, ARRAY_SIZE);

 clicksStart = clock();
 shellSort(sorted,ARRAY_SIZE);
 clicksEnd = clock();
 cout << "shellSort time taken: " << (float)(clicksEnd - clicksStart)/CLOCKS_PER_MS << "ms\n";

 copyArr (toSort, sorted, ARRAY_SIZE);

 clicksStart = clock();
 mergeSort(sorted,0,ARRAY_SIZE -1);
 clicksEnd = clock();
 cout << "merge sort time taken: " << (float)(clicksEnd - clicksStart)/CLOCKS_PER_MS << "ms\n";
 
 copyArr (toSort, sorted, ARRAY_SIZE);

 clicksStart = clock();
 quickSort(sorted,0,ARRAY_SIZE-1);
 clicksEnd = clock();

 cout << "quick sort time taken: " << (float)(clicksEnd - clicksStart)/CLOCKS_PER_MS << "ms\n";

 max = sorted[ARRAY_SIZE-1];

 copyArr (toSort, sorted, ARRAY_SIZE);

 clicksStart = clock();
 bucketSort(sorted,ARRAY_SIZE,max+1);
 clicksEnd = clock();
 cout << "bucket sort time taken: " << (float)(clicksEnd - clicksStart)/CLOCKS_PER_MS << "ms\n";

 //printArr(sorted,ARRAY_SIZE);

 copyArr (toSort, sorted, ARRAY_SIZE);

 clicksStart = clock();
 radixSort(sorted,ARRAY_SIZE);
 clicksEnd = clock();
 cout << "radix sort time taken: " << (float)(clicksEnd - clicksStart)/CLOCKS_PER_MS << "ms\n";

 copyArr (toSort, sorted, ARRAY_SIZE);

 clicksStart = clock();
 qsort(sorted,ARRAY_SIZE,sizeof(uint32_t),compare);
 clicksEnd = clock();
 cout << "qsort time taken: " << (float)(clicksEnd - clicksStart)/CLOCKS_PER_MS << "ms\n";

 copyArr (toSort, sorted, ARRAY_SIZE);

 clicksStart = clock();
 sort(sorted,sorted+ARRAY_SIZE);
 clicksEnd = clock();
 cout << "sort time taken: " << (float)(clicksEnd - clicksStart)/CLOCKS_PER_MS << "ms\n";

 copyArr (toSort, sorted, ARRAY_SIZE);

 clicksStart = clock();
 newSize = bucketDedupSort(sorted,ARRAY_SIZE,max+1);
 clicksEnd = clock();
 cout << "bucketDedupSort time taken: " << (float)(clicksEnd - clicksStart)/CLOCKS_PER_MS << "ms\n";

 //printArr(sorted,newSize);

 return 0;
}
