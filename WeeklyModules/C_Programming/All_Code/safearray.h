/* My safe array library. A library for working safely with arrays*/
struct safe_array_int {
	int *array;
	int nelems;
	int array_size;
};

typedef struct safe_array_int isafe_t;

/* Function prototypes */
isafe_t *create_safe_intarray(int size);
void destroy_safe_intarray(isafe_t *oldarray);
int set_value_iarray(isafe_t *myarray, int index, int val);
int get_value_iarray(isafe_t *myarray, int index, int* result);
