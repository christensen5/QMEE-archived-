#include <stdio.h>

int main (void)
{

	int a = 5;
	int b = 2;
	int c = 0;
	int d = 0;

	// && AND (binary: two operands)
	if (a && b) {
		printf("%i and %i are non-zero\n", a, b);
	}
	
	if (a && c) {
		printf("%i and %i are non-zero\n", a, c);
	}
	else {
		printf("One of %i or %i is zero\n", a, c);
	}	

	// || OR  (binary: two operands)
	if (a || b) {
		printf("At least one of %i or %i is non-zero\n", a, b);
	}	

	if (a || c) {
		printf("At least one of %i or %i is non-zero\n", a, c);
	}
	else {
		printf("Both %i and %i are zero\n", a , c);
	}

	//  ! NOT (unary: one operand)


}
