#include <stdio.h>

int main (void)
{

	int x = 0;
	char a = 'a';

	x = x + 1;
	a = a + 1;

	printf("x is: %i; a is: %c\n", x, a);

	// ++ --

	// ++x or x++ -> x = x + 1
	// ++x:
	// x = x + 1;
	// return x;

	// x++
	// return x;
	// x = x + 1;
	
	// ++x -> return value of x; then increment x
	printf("x++ is: %i; a is: %c\n", x++, ++a);
	printf("x is: %i; a is: %c\n", x, ++a);
	
	return 0;
}
