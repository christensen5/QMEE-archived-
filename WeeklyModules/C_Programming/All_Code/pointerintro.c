/* A program to illustrate elementary pointer operations */
#include <stdio.h>

int main (void)
{
	char c = 'M';
	char b = 'n';
	char *char_ptr = NULL;	
	
	printf("char_ptr: %p\n", char_ptr);	

	char_ptr = &c;

	printf("char_ptr after initialisation: %p\n", char_ptr);	

	printf("dereferencing char_ptr: %c\n", *char_ptr);

	*char_ptr = 'J';
		
	printf("dereferencing char_ptr again: %c\n", *char_ptr);

	return 0;
}
