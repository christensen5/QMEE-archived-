#include <stdio.h>

int* return_int_ptr(void)
{
	/* This is a dumb function don't do it */
	/* At any point in program exection, this could be overwritten*/
	int i = 42;
	int* iptr;
	
	iptr = &i;

	return iptr;
}

int main (void)
{
	int locali;
	int *localiptr;

	localiptr = return_int_ptr();
	locali = *return_int_ptr();
	
	printf("Result of return from fn: %i\n", *localiptr);
	printf("Result of return from fn: %i\n", locali);

	return 0;	
}
