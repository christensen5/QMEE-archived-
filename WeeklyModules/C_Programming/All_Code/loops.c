#include <stdio.h>
#include <limits.h>

int main (void)
{	
	int i;

	int im_getting_bored = 50000;

	i = 0;
	while (i < 5) {
		i = i + 1; //++i;
	}	
	printf("after while looping i is: %i\n", i);


	printf("Running %i loops\n", INT_MAX);
	i = 0;
	do {
		++i;
		//printf("value of i: %i\n", i);
		
		continue;		
	
		if (i == im_getting_bored) {
			printf("Okay, I'm getting bored of this crap!\n");
			break;
		}
	} 
	while (i < INT_MAX);
	printf("after do-while looping i is: %i\n", i);

	for(i = 0; i < 5 ; ++i) {
		printf("for loops are nice\n");
	}

	return 0;
}
