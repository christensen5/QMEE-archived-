#include <stdio.h>

int main (void)
{
	
	int x = 4;
	int y = 7;
	int i_result = 0;

	float f_y = 7.0;	
	float f_result = 0.0;
	
	// + - / * %

	i_result = y / x;
	f_result = (float)y / x;

	printf("i_result: %i\n", i_result);
	printf("f_result: %f\n", f_result);

	return 0;	

}
