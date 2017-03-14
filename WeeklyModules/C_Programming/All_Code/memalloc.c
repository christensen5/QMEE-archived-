#include <stdio.h>
#include <stdlib.h>

// malloc: memory allocation
// calloc: cleared memory allocation

int *create_array(int nelems)
{
	int *arrayptr;
	arrayptr = (int*)malloc(nelems * sizeof(int)); 
	//arrayptr = (int*)calloc(nelems, sizeof(int));	

	int i = 0;
	for (i = 0; i < nelems; ++i){
		printf("%i ", arrayptr[i]);
	}
	printf("\n\n");
	return arrayptr;
}


int main (void)
{
	int nelems = 220;
	//int array[100];
	//int array2[nelems];

	int *arrayptr = create_array(nelems);
	int i = 0;
	for (i = 0; i < nelems; ++i){
		printf("%i ", arrayptr[i]);
	}
	printf("\n\n");
	
	free(arrayptr);

	return 0;	
}
