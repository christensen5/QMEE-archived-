#include <stdio.h>
#include "safearray.h"

int main (void)
{
	int i = 0;
	int size = 10;
	isafe_t *newarray;

	newarray = create_safe_intarray(size);

	// Interesting stuff happens
	for (i = 0; i < size; ++i) {
		set_value_iarray(newarray, i, i + 1);
	}

	i = 0;
	int result;
	while(!get_value_iarray(newarray, i, &result)) {
		printf("Result: %i\n", result);
		++i;
	};

	destroy_safe_intarray(newarray);	
	
	return 0;
}
