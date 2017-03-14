#include <stdio.h>

int main (void)
{
	printf("This is an example of a string\n\n");
	
	int nelems = 10;
	char carray [] = {'n', 'o', 'w', ' ', 'a', ' ', 's', '\0'};
	char carray2 [nelems];
	char carray3 [] = "this is a \"string\"\n";

	printf("%s", carray3);
	printf("\n");
	printf("%s", carray);
	printf("\n");

	int i = 0;
	while(carray3[i]) {
		printf("%c", carray3[i]);
		++i;
	}
	printf("\n");	
	
	return 0;
		
}
