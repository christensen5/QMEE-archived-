#include <stdlib.h>
#include "safearray.h"


isafe_t *create_safe_intarray(int size)
{	
	isafe_t *numbers;
	
	numbers = (isafe_t*)malloc(sizeof(isafe_t)); 

	numbers->nelems = 0;
	numbers->array_size = size;
	numbers->array = (int*)malloc(numbers->array_size * sizeof(int));

	return numbers;
}

void destroy_safe_intarray(isafe_t *oldarray)
{
	if (oldarray->array) {
		free(oldarray->array);
	}
	oldarray->array = NULL;	
}

int set_value_iarray(isafe_t *myarray, int index, int val)
{
	if (myarray) {
		if (index < myarray->array_size) {
			myarray->array[index] = val;
			return 0;
		}
		return 1;
	}
	return -1;
}

int get_value_iarray(isafe_t *myarray, int index, int* result)
{
	if (myarray) {
		if (index < myarray->array_size) {
			*result = myarray->array[index];
			return 0;
		}
		return 1;
	}
	return -1;
}
