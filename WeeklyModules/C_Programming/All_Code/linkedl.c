#include <stdio.h>

struct entry {
	int index;
	struct entry *next;
};

typedef struct entry entry_t;

void traverse_list(entry_t *n)
{
	//printf("Visiting entry %i\n", n->index);	
	if (n->next) {
		traverse_list(n->next);	
	}
	printf("Visiting entry %i\n", n->index);	
}

void delete_next_element(entry_t *n)
{
	entry_t *temp;

	temp = n->next->next;

	n->next->next = NULL;
	n->next = temp;
}

int main (void)
{
	int i = 0;
	entry_t n1, n2, n3;
	
	n1.index = 1; n2.index = 2; n3.index = 3;
	n1.next = &n2; n2.next = &n3; n3.next = NULL;

	// Access the index value of structure pointed to by n1.next
	i = n1.next->next->index;	

	traverse_list(&n1);	
	delete_next_element(&n1);
	printf("\n");
	traverse_list(&n1);
	
	return 0;
}
