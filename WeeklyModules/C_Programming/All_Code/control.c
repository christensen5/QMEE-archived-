#include <stdio.h>

int main (void)
{
	// Must distinguish = (assignment) from == ('is equal to')

	int i = 0;

	if (i != 1) {
		printf ("i is not 1\n");
	}
	else {
		printf("i is 1\n");
	}

	if ((i = 3)) {
		
		printf("Attempt to assign %i to i\n", i);
	}

	// a <= b;  a => b;

	int a = 2;
	int b = 3;

	if (a <= b) {
		printf("a is less or equal to b\n");
	}



	return 0;
}
