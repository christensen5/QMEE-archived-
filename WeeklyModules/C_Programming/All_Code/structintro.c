#include <stdio.h>

struct site_data {
	float lat;
	float longit;
	int index;
	struct site_data *next;
};

// typedef <datatype> <alias>
typedef struct site_data site_data_t;


void print_struct_index(site_data_t s)
{
	printf("Site index: %i\n", s.index);
}

int main (void)
{
	site_data_t s1;
	site_data_t s2;

	s1.index = 1;
	s2.index = 2;

	print_struct_index(s1);
	
	return 0;

}
