#include <stdio.h>

int set_bit(int tgt, int bitpos)
{
	if (!bitpos) {
		printf("Error: can't input 0 bitpos\n");
		return tgt;
	}	

	tgt = tgt | (1 << (bitpos - 1));

	return tgt;	
}

int main (void)
{

	int a = 0;
	printf("a is: %i\n", a);

	a = set_bit(a, 4);
	printf("a is now: %i\n", a);

	a = set_bit(a, 27);
	printf("a is now: %i\n", a);

	return 0;
}
