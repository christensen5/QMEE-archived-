/* A program to deliver the bit width of a byte on this
 * machine */
#include <stdio.h>

int count_bits(void)
{
	int i = 0;
	int shifter = 0;
	char twochars[] = {~0, 0};
	
	do {
		if (twochars[1] & (1 << i)) {
			break;
		}
		++i;	
	}while (i < 64);

	printf("The number of bits %i\n", i);
	
	return i;
	
}

int main (void)
{
	count_bits();
	return 0;
}
