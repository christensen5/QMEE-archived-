#include <stdio.h>

void print_int_array(int inputarray[], int nelems)
{
	int i;
	
	for (i = 0; i < nelems; ++i) {
		printf("%i ", inputarray[i]);
	}
	printf("\n");
}

int sum_array_values(int inputarray[], int nelems)
{
	int i = 0;
	int sum = 0;

	for (i = 0; i < nelems; ++i) {
		sum = sum + inputarray[i];
		// sum += inputarray[i];
	}
	
	return sum;
}

int main (void)
{
	int n_elems = 5;
	int intarray [n_elems];

	print_int_array(intarray, n_elems);
	printf("Sum of array: %i\n", 
			sum_array_values(
							intarray, 
							n_elems));

	return 0;

}
