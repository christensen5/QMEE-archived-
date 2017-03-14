#include <stdio.h>

void read_array_with_ptr(int *iptr, int nelems)
{
	int i = 0;
	for (i = 0; i < nelems; ++i){
		printf("%i ", *iptr);
		iptr = iptr + 1;
	}
	printf("\n");	
}

int main(void)
{
	int array[] = {72, 1, 22, 67, 98, 34};
	int nelems = 6;
	int *intptr;
	int *intptr2;
	char *charptr;

	intptr = array;
	//charptr = (char)array;
	read_array_with_ptr(array, nelems);
	return 0;
}
