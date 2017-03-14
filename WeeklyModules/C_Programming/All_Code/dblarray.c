#include <stdio.h>
#include "./printarray/p_array.h"

void double_array_values(float fparray[], int intarray[], int nelems)
{
    int i;

    if (fparray) {
        for (i = 0; i < nelems; ++i) {
            fparray[i] = fparray[i] * 2;
        }
    }
    else if (intarray) {
        for (i = 0; i < nelems; ++i) {
            intarray[i] = intarray[i] * 2;
        }
    }
}

int main (void)
{
	int nelems = 4;
	int intarray[] = {81, 8, 4, 30};
	float fparray[] = {2.30, 10.1, 10.0, 81.8};

	printf("Before pass to doubler:\n");
	p_int_array(intarray, nelems);
	printf("\n");
	p_float_array(fparray, nelems);
	printf("\n");	
		
	double_array_values(NULL, intarray, nelems);
	double_array_values(fparray, NULL, nelems);

	printf("\nAfter pass to doubler:\n");
	p_int_array(intarray, nelems);
	printf("\n");
	p_float_array(fparray, nelems);
	printf("\n");	

	return 0;

}
